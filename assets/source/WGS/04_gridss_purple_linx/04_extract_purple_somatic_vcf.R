#!/usr/bin/env Rscript

# module use /usr/local/package/modulefiles/
#module load R/4.0.2

working_dir <- getwd()

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "~/workspace/runs", study.dir )

source( file.path( run.dir, "config_files/common_config.R" ) )

suppressPackageStartupMessages({
	library(data.table)
	library(parallel)
	library(VariantAnnotation)})

germline.dir <- paste(run.dir, "/result/", gpl_prefix, sep="")
purple.dir <- paste(run.dir, "/result/", gpl_prefix, sep="")
ref.dir <- "/home/tjohnson/reference/HMF/38/dbs/ensembl_data_cache"

if ( length(list.dirs(paste(run.dir, "/result_summaries/", gpl_prefix, sep=""), recursive=FALSE))==0 ){
	system(paste("mkdir ", run.dir, "/result_summaries/", gpl_prefix, sep=""))
}

load( file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))

tumor.normal.pair.info <- tumor_normal_pair_info.GPL

tumor.sample.ids.ls <- tumor.normal.pair.info$sample.id.T
names(tumor.sample.ids.ls) <- tumor.sample.ids.ls
sample.ct <- length(tumor.sample.ids.ls)

threads <- min( c(18, sample.ct) )

chrom.map.ls <- as.integer(1:25)
names(chrom.map.ls) <- paste("chr", c(1:22, "X", "Y", "M"), sep="")

extract_SEC <- function( curr.SEC ){
	curr.SEC.ls <- unlist(curr.SEC)
	if (length(curr.SEC.ls)==6){
		data.table(
			gene = curr.SEC.ls[[1]],
			ENST = curr.SEC.ls[[2]],
			funcClass1 = curr.SEC.ls[[3]],
			funcClass2 = curr.SEC.ls[[4]],
			bp.descr = curr.SEC.ls[[5]],
			aa.descr = curr.SEC.ls[[6]])
	}else{
		data.table(
			gene = "",
			ENST = "",
			funcClass1 = "",
			funcClass2 = "",
			bp.descr = "",
			aa.descr = "")
	}
}

extract_SEW <- function( curr.SEW ){
	curr.SEW.ls <- unlist(curr.SEW)
	if (length(curr.SEW.ls)==5){
		data.table(
			SEW.gene = curr.SEW.ls[[1]],
			SEW.ENST = curr.SEW.ls[[2]],
			SEW.funcClass1 = curr.SEW.ls[[3]],
			SEW.funcClass2 = curr.SEW.ls[[4]],
			SEW.genes = curr.SEW.ls[[5]])
	}else{
		data.table(
			SEW.gene = "",
			SEW.ENST = "",
			SEW.funcClass1 = "",
			SEW.funcClass2 = "",
			SEW.genes = "")
	}
}

extract_LOF <- function( curr.LOF ){
	curr.LOF.ls <- unlist(curr.LOF)
	if (length(curr.LOF.ls)==5){
		data.table(
			LOF.gene = curr.LOF.ls[[1]],
			LOF.ENSG = curr.LOF.ls[[2]],
			LOF.ENST = curr.LOF.ls[[3]],
			LOF.gene.tx.ct = curr.LOF.ls[[4]],
			LOF.pct.affected.Tx = curr.LOF.ls[[5]])
	}else{
		data.table(
			LOF.gene = "",
			LOF.ENSG = "",
			LOF.ENST = "",
			LOF.gene.tx.ct = -99L,
			LOF.pct.affected.Tx = -99)
	}
}

extract_CLN <- function(curr.value){
	if ( curr.value=='character(0)' ){
		''
	}else{
		paste(curr.value, collapse=";")
	}
}


message( paste("Importing somatic variants for ", study.name, sep=""))

somatic.variant.info.ls <- mclapply(
	X = tumor.sample.ids.ls,
	FUN = function( curr.sample.id ){
##		curr.sample.id <- tumor.sample.ids.ls[[2]]
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		curr.normal.sample.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$sample.id.N
		curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
		curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
		
		message(paste("Extracting somatic vcf data for ", curr.sample.id, sep=""))
		file.gz <- paste(purple.dir, "/", curr.subject.id, "/purple/", curr.sample.id, "/", curr.sample.id, ".purple.somatic.vcf.gz", sep="")
		vcf <- readVcf(file.gz, "GRCh38")
	
		vcf.header.dt <- data.table(
			CHROM = as.character(seqnames(rowRanges(vcf))),
			ID = rownames( geno(vcf)$GT),
			start = as.integer(start(rowRanges(vcf))),
			end = as.integer(end(rowRanges(vcf))),
			REF = as.character(rowRanges(vcf)$REF),
			ALT = as.data.table(rowRanges(vcf)$ALT)$value,
			FILTER = as.character(rowRanges(vcf)$FILTER),
			QUAL = as.integer(rowRanges(vcf)$QUAL))

		vcf.header.dt[,chr:=chrom.map.ls[CHROM]]
		vcf.header.dt[,snp.id:=paste(chr, "-", formatC(start, format="d"), "-", REF, "-", ALT, sep="")]

		message(paste("Merging vcf header and variant data for ", curr.sample.id, sep=""))
		vcf.dt <- as.data.table(info(vcf))
		vcf.dt <- cbind(vcf.header.dt[,list(snp.id, CHROM, chr, ID, start, end, REF, ALT, FILTER, QUAL)], vcf.dt)
		vcf.dt[,subject.id:=curr.subject.id]
		vcf.dt[,tumor.sample.id:=curr.sample.id]
		vcf.dt[,subject.index:=curr.subject.index]
		vcf.dt[,tumor.index:=curr.tumor.index]

		if (curr.normal.sample.id=='none'){
			vcf.dt$GT.N <- as.character(NA)
		}else{
			vcf.dt$GT.N <- geno(vcf)$GT[,curr.normal.sample.id]
		}
		
		vcf.dt$GT.T <- geno(vcf)$GT[,curr.sample.id]
		
		vcf.dt[,AF:=unlist(AF)]
		if (nrow(vcf.dt[AF=='.',])>0){
			vcf.dt[which(AF=='.'), AF:=NA]	
		}
		vcf.dt[,AF:=as.numeric(AF)]
		
		vcf.dt[,AF_raw:=unlist(AF_raw)]
		if (nrow(vcf.dt[AF_raw=='.',])>0){
			vcf.dt[which(AF_raw=='.'), AF_raw:=NA]	
		}
		vcf.dt[,AF_raw:=as.numeric(AF_raw)]
		
		vcf.dt[,AF_eas:=unlist(AF_eas)]
		if (nrow(vcf.dt[AF_eas=='.',])>0){
			vcf.dt[which(AF_eas=='.'), AF_eas:=NA]	
		}
		vcf.dt[,AF_eas:=as.numeric(AF_eas)]
		
		if (curr.normal.sample.id=='none'){
			vcf.dt$AF.N <- as.numeric(NA)
		}else{
			vcf.dt$AF.N <- geno(vcf)$AF[,curr.normal.sample.id]
		}
		
		vcf.dt$AF.T <- geno(vcf)$AF[,curr.sample.id]
		
		AF.names.to.parse <- c('dbGaP_PopFreq', 'TOPMED', 'GnomAD', 'GnomAD_exomes', 'X1000Genomes', 'TOMMO', 'KOREAN', 'Korea1K')
		
		for (curr.AF.name in AF.names.to.parse){
			curr.AF.data <- unlist(copy(vcf.dt[[curr.AF.name]]))
			vcf.dt[[curr.AF.name]] <- curr.AF.data[seq(2,length(curr.AF.data),2)]
			vcf.dt[[curr.AF.name]] <- ifelse(is.na(vcf.dt[[curr.AF.name]]), -99, vcf.dt[[curr.AF.name]])
		}

		vcf.dt[,CLNSIG:=extract_CLN(CLNSIG), by=1:nrow(vcf.dt)]
		vcf.dt[,CLNSIGCONF:=extract_CLN(CLNSIGCONF), by=1:nrow(vcf.dt)]

		vcf.dt[,CLNSIG:=unlist(CLNSIG)]
		vcf.dt[,CLNSIGCONF:=unlist(CLNSIGCONF)]
		
		message(paste("Finished extracting somatic vcf data for ", curr.sample.id, sep=""))
		
		base.vcf.dt <- vcf.dt[,list(subject.id, tumor.sample.id, subject.index, tumor.index,
			snp.id, ID, CHROM, chr, start, end, REF, ALT, FILTER, QUAL, REPORTED,
			BIALLELIC, CLNSIG, CLNSIGCONF,
			LPS, LRS, GT.N, GT.T,
			AF.N, AF.T, PURPLE_AF, PURPLE_CN, PURPLE_GERMLINE, PURPLE_MACN, PURPLE_VCN, SUBCL, TIER,
			HOTSPOT, NEAR_HOTSPOT,
			MH, RC_MH, RC, RC_IDX, RC_LF, RC_RF,
			AF_TOT, AF_EAS, AF, AF_raw, AF_eas,
			RS, VC, dbGaP_PopFreq, TOPMED, GnomAD, GnomAD_exomes, X1000Genomes, TOMMO, KOREAN, Korea1K)]

		vcf.dt.SEC <- vcf.dt[,extract_SEC( SEC ), by = list(subject.id, tumor.sample.id, subject.index, tumor.index, snp.id, ID, CHROM, chr, start, end, REF, ALT)]
		vcf.dt.SEW <- vcf.dt[,extract_SEW( SEW ), by = list(subject.id, tumor.sample.id, subject.index, tumor.index, snp.id, ID, CHROM, chr, start, end, REF, ALT)]
		vcf.dt.LOF <- vcf.dt[,extract_LOF( LOF ), by = list(subject.id, tumor.sample.id, subject.index, tumor.index, snp.id, ID, CHROM, chr, start, end, REF, ALT)]

		vcf.dt <- merge(
			x = base.vcf.dt,
			y = vcf.dt.SEC,
			by = c('subject.id', 'tumor.sample.id', 'subject.index', 'tumor.index', 'snp.id', 'ID', 'CHROM', 'chr', 'start', 'end', 'REF', 'ALT'))
	
		vcf.dt <- merge(
			x = vcf.dt,
			y = vcf.dt.SEW,
			by = c('subject.id', 'tumor.sample.id', 'subject.index', 'tumor.index', 'snp.id', 'ID', 'CHROM', 'chr', 'start', 'end', 'REF', 'ALT'))
	
		vcf.dt <- merge(
			x = vcf.dt,
			y = vcf.dt.LOF,
			by = c('subject.id', 'tumor.sample.id', 'subject.index', 'tumor.index', 'snp.id', 'ID', 'CHROM', 'chr', 'start', 'end', 'REF', 'ALT'))

		vcf.dt[order(chr, start),]
	},
	mc.cores = threads)

somatic.variant.info.dt <- rbindlist(somatic.variant.info.ls)
somatic.driver.variant.info.dt <- somatic.variant.info.dt[REPORTED==TRUE,]

save(list=c('somatic.variant.info.dt'),
	file = paste(run.dir, "/result_summaries/", gpl_prefix, "/somatic.variant.info.Rdata", sep=""))

save(list=c('somatic.driver.variant.info.dt'),
	file = paste(run.dir, "/result_summaries/", gpl_prefix, "/somatic.driver.variant.info.Rdata", sep=""))

message( paste("Saved PURPLE processed somatic variant data for ", study.name, sep=""))
