#!/usr/bin/env Rscript

# module use /usr/local/package/modulefiles/
#module load R/4.0.2

args <- commandArgs(trailingOnly = TRUE)

run.dir <- args[[1]]
#run.dir <- "/home/tjohnson/workspace/runs/KakimiLab_WGS_HMF_20210426"

threads <- args[[2]]
curr.sample.index <- as.integer(args[[3]])

source( file.path( run.dir, "config_files/common_config.R" ) )

message( paste( 'Running extraction of PURPLE germline genomewide VCF for study ', study.name, ': sample index ', curr.sample.index, sep="") )

suppressPackageStartupMessages({
	library(data.table)
	library(parallel)
	library(VariantAnnotation)})

germline.dir <- paste(run.dir, "/result/", gpl_prefix, sep="")
purple.dir <- paste(run.dir, "/result/", gpl_prefix, sep="")
ref.dir <- "/home/tjohnson/reference/HMF_38/dbs/ensembl_data_cache"

if ( length(list.dirs(paste(run.dir, "/result_summaries/", gpl_prefix, sep=""), recursive=FALSE))==0 ){
	system(paste("mkdir ", run.dir, "/result_summaries/", gpl_prefix, sep=""))
}

load( file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))

tumor.normal.pair.info <- tumor_normal_pair_info.GPL[which(sample.id.N!='none',)]

tumor.sample.ids.ls <- tumor.normal.pair.info$sample.id.T
names(tumor.sample.ids.ls) <- tumor.sample.ids.ls
sample.ct <- length(tumor.sample.ids.ls)

curr.sample.id <- tumor.sample.ids.ls[[curr.sample.index]]

threads <- 12L

chroms.ls <- paste("chr", c(1:22, "X", "Y", "M"), sep="")

chrom.map.ls <- as.integer(1:25)
names(chrom.map.ls) <- chroms.ls

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


message( paste("Importing germline variants for ", study.name, sep=""))

curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
curr.normal.sample.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$sample.id.N
curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index

message(paste("Extracting germline vcf data for ", curr.sample.id, sep=""))
file.gz <- paste(purple.dir, "/", curr.subject.id, "/purple/", curr.sample.id, "/", curr.sample.id, ".purple.germline.vcf.gz", sep="")

germline.variant.info.dt.ls <- mclapply(
	X = chrom.map.ls,
	FUN = function( curr.chr ){
		curr.chrom <- chroms.ls[[curr.chr]]

		message(paste("Processing VCF for ", curr.chrom, " and ", curr.sample.id, sep=""))
		
		rng <- GRanges(seqnames=c(curr.chrom), IRanges(c(1), c(250000000)) )
		
		tab <- TabixFile(file.gz)
		vcf <- readVcf(tab, "GRCh38", param=rng)

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
		
		message(paste("Merging vcf header and variant data for ", curr.sample.id, " and ", curr.chrom, sep=""))
		vcf.dt <- as.data.table(info(vcf))

		AF.names.to.parse <- c('dbGaP_PopFreq', 'TOPMED', 'GnomAD', 'GnomAD_exomes', 'X1000Genomes', 'TOMMO', 'KOREAN', 'Korea1K')
		
		for (curr.AF.name in AF.names.to.parse){
			curr.AF.data <- unlist(copy(vcf.dt[[curr.AF.name]]))
			vcf.dt[[curr.AF.name]] <- curr.AF.data[seq(2,length(curr.AF.data),2)]
			vcf.dt[[curr.AF.name]] <- ifelse(is.na(vcf.dt[[curr.AF.name]]), -99, vcf.dt[[curr.AF.name]])
		}

		vcf.dt[,AF_TOT:=ifelse(is.na(AF_TOT), -99, AF_TOT)]
		vcf.dt[,AF_EAS:=ifelse(is.na(AF_EAS), -99, AF_EAS)]
		vcf.dt[,MAF_TOT:=ifelse(is.na(MAF_TOT), -99, MAF_TOT)]
		vcf.dt[,MAF_EAS:=ifelse(is.na(MAF_EAS), -99, MAF_EAS)]
		
		vcf.dt <- cbind(vcf.header.dt[,list(snp.id, CHROM, chr, ID, start, end, REF, ALT, FILTER, QUAL)], vcf.dt)
		
		rm(vcf.header.dt)
		gc()

		vcf.dt[,subject.id:=curr.subject.id]
		vcf.dt[,tumor.sample.id:=curr.sample.id]
		vcf.dt[,subject.index:=curr.subject.index]
		vcf.dt[,tumor.index:=curr.tumor.index]

		message(paste("Extracting genotype data for and CLNSIG annotations ", curr.sample.id, " and ", curr.chrom, sep=""))
		
		vcf.dt$GT.N <- geno(vcf)$GT[,curr.normal.sample.id]
		
		vcf.dt$AF.N <- geno(vcf)$AF[,curr.normal.sample.id]
		vcf.dt$AF.T <- geno(vcf)$AF[,curr.sample.id]
		
		vcf.dt$DP.N <- geno(vcf)$DP[,curr.normal.sample.id]
		vcf.dt$DP.T <- geno(vcf)$DP[,curr.sample.id]
		
		vcf.dt[,CLNSIG:=extract_CLN(CLNSIG), by=1:nrow(vcf.dt)]
		vcf.dt[,CLNSIGCONF:=extract_CLN(CLNSIGCONF), by=1:nrow(vcf.dt)]
		
		vcf.dt[,CLNSIG:=unlist(CLNSIG)]
		vcf.dt[,CLNSIGCONF:=unlist(CLNSIGCONF)]

		setkey(vcf.dt, chr, start, end, REF, ALT)

		rm(vcf)
		gc()

		curr.vcf.dt <- vcf.dt[,list(subject.id, tumor.sample.id, subject.index, tumor.index, snp.id, ID, CHROM, chr, start, end, REF, ALT, SEC, SEW, LOF)]

		message(paste("Extracting SEC annotations ", curr.sample.id, " and ", curr.chrom, sep=""))
		vcf.dt.SEC <- curr.vcf.dt[,extract_SEC( SEC ), by = list(subject.id, tumor.sample.id, subject.index, tumor.index, snp.id, ID, CHROM, chr, start, end, REF, ALT)]
		germline.variant.info.dt <- merge(
			x = vcf.dt[,list(subject.id, tumor.sample.id, subject.index, tumor.index,
				snp.id, ID, CHROM, chr, start, end, REF, ALT, FILTER, QUAL, REPORTED,
				BIALLELIC, CLNSIG, CLNSIGCONF, AF_TOT, AF_EAS, MAF_TOT, MAF_EAS,
				RS, VC, dbGaP_PopFreq, TOPMED, GnomAD, GnomAD_exomes, X1000Genomes, TOMMO, KOREAN, Korea1K,
				GT.N, AF.N, AF.T, DP.N, DP.T,
				LPS, LRS, PURPLE_AF, PURPLE_CN, PURPLE_MACN, PURPLE_VCN, TIER,
				HOTSPOT, NEAR_HOTSPOT, NMD, PATH)],
			y = vcf.dt.SEC,
			by = c('subject.id', 'tumor.sample.id', 'subject.index', 'tumor.index', 'snp.id', 'ID', 'CHROM', 'chr', 'start', 'end', 'REF', 'ALT'))

		rm(vcf.dt)
		rm(vcf.dt.SEC)
		gc()

		message(paste("Extracting SEW annotations ", curr.sample.id, " and ", curr.chrom, sep=""))
		vcf.dt.SEW <- curr.vcf.dt[,extract_SEW( SEW ), by = list(subject.id, tumor.sample.id, subject.index, tumor.index, snp.id, ID, CHROM, chr, start, end, REF, ALT)]
		germline.variant.info.dt <- merge(
			x = germline.variant.info.dt,
			y = vcf.dt.SEW,
			by = c('subject.id', 'tumor.sample.id', 'subject.index', 'tumor.index', 'snp.id', 'ID', 'CHROM', 'chr', 'start', 'end', 'REF', 'ALT'))

		rm(vcf.dt.SEW)
		gc()

		message(paste("Extracting LOG annotations ", curr.sample.id, " and ", curr.chrom, sep=""))
		vcf.dt.LOF <- curr.vcf.dt[,extract_LOF( LOF ), by = list(subject.id, tumor.sample.id, subject.index, tumor.index, snp.id, ID, CHROM, chr, start, end, REF, ALT)]
		germline.variant.info.dt <- merge(
			x = germline.variant.info.dt,
			y = vcf.dt.LOF,
			by = c('subject.id', 'tumor.sample.id', 'subject.index', 'tumor.index', 'snp.id', 'ID', 'CHROM', 'chr', 'start', 'end', 'REF', 'ALT'))
		
		rm(vcf.dt.LOF)
		gc()
		
		germline.variant.info.dt
	},
	mc.cores = threads)

germline.variant.info.dt <- rbindlist(germline.variant.info.dt.ls)
rm(germline.variant.info.dt.ls)
gc()

save(list=c('germline.variant.info.dt'),
	file = paste(purple.dir, "/", curr.subject.id, "/purple/", curr.sample.id, "/", curr.sample.id, ".purple.germline.vcf.Rdata", sep=""))

message( paste("Saved ", nrow(germline.variant.info.dt), " PURPLE processed genomewide germline variants for ", curr.sample.id, sep=""))

germline.variant.info.dt <- germline.variant.info.dt[REPORTED=="TRUE" | TIER=="PANEL" | funcClass2%in%c('MISSENSE', 'NONSENSE_OR_FRAMESHIFT', 'SPLICE', 'SYNONYMOUS') | SEW.funcClass2%in%c('MISSENSE', 'NONSENSE_OR_FRAMESHIFT', 'SPLICE', 'SYNONYMOUS'),]

save(list=c('germline.variant.info.dt'),
	file = paste(purple.dir, "/", curr.subject.id, "/purple/", curr.sample.id, "/", curr.sample.id, ".purple.germline.vcf.filtered.Rdata", sep=""))

message( paste("Saved ", nrow(germline.variant.info.dt), " filtered PURPLE processed genomewide germline variants for ", curr.sample.id, sep=""))
