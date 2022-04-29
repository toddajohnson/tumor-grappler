#!/usr/bin/env Rscript

#module use /usr/local/package/modulefiles/
#module load R/4.0.2

working_dir <- getwd()
gpl_prefix <- basename(working_dir)

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
ref.dir <- "/home/tjohnson/reference/HMF_38/dbs/ensembl_data_cache"

if ( length(list.dirs(paste(run.dir, "/result_summaries/", gpl_prefix, sep=""), recursive=FALSE))==0 ){
	system(paste("mkdir ", run.dir, "/result_summaries/", gpl_prefix, sep=""))
}

load( file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))

tumor.normal.pair.info <- tumor_normal_pair_info.GPL
tumor.sample.ids.ls <- tumor.normal.pair.info$sample.id.T
names(tumor.sample.ids.ls) <- tumor.sample.ids.ls
sample.ct <- length(tumor.sample.ids.ls)

threads <- min( c(18, sample.ct) )

message( paste("Read tumor-normal information with ", sample.ct, " samples for ", study.name, sep=""))

message( paste("Importing QC tables for ", study.name, sep=""))

purple.qc.ls <- mclapply(
	X = tumor.sample.ids.ls,
	FUN = function( curr.sample.id ){
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		curr.normal.sample.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$sample.id.N
		curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
		curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
		
		qc.file <- paste(purple.dir, "/", curr.subject.id, "/purple/", curr.sample.id, "/", curr.sample.id, ".purple.qc", sep="")
		curr.dt <- fread(qc.file, header = FALSE, col.names=c('Variable', 'Value'))
		curr.cols <- copy(colnames(curr.dt))
		curr.dt[,subject.id:=curr.subject.id]
		curr.dt[,tumor.sample.id:=curr.sample.id]
		curr.dt[,normal.sample.id:=curr.normal.sample.id]
		curr.dt[,subject.index:=curr.subject.index]
		curr.dt[,tumor.index:=curr.tumor.index]
		
		dcast(
			data = curr.dt,
			formula = subject.id + tumor.sample.id + normal.sample.id + subject.index + tumor.index ~ Variable,
			value.var = "Value")
	},
	mc.cores = threads)

purple.qc.dt <- rbindlist(purple.qc.ls)

message( paste("Merged QC tables for ", study.name, "\n", sep=""))

message( paste("Importing purity tables for ", study.name, sep=""))

purple.purity.ls <- mclapply(
	X = tumor.sample.ids.ls,
	FUN = function( curr.sample.id ){
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
		curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
		
		purity.file <- paste(purple.dir, "/", curr.subject.id, "/purple/", curr.sample.id, "/", curr.sample.id, ".purple.purity.tsv", sep="")
		curr.dt <- fread(purity.file, header = TRUE)
		curr.cols <- copy(colnames(curr.dt))
		curr.dt[,subject.id:=curr.subject.id]
		curr.dt[,tumor.sample.id:=curr.sample.id]
		curr.dt[,subject.index:=curr.subject.index]
		curr.dt[,tumor.index:=curr.tumor.index]
		curr.dt[,c('subject.id', 'tumor.sample.id', 'subject.index', 'tumor.index', curr.cols), with=FALSE]
	},
	mc.cores = threads)

purple.purity.dt <- rbindlist(purple.purity.ls)

message( paste("Merged purity tables for ", study.name, "\n", sep=""))

purple.qc.dt.failed.contamination <- purple.qc.dt[grep("FAIL_CONTAMINATION", QCStatus, fixed=TRUE),]

message( paste("Saving QC and purity tables for ", study.name, "\n", sep=""))

fwrite(purple.qc.dt,
		file = paste(run.dir, "/result_summaries/", gpl_prefix, "/purple.qc.summary.txt", sep=""),
		col.names = TRUE,
		sep = '\t')

fwrite(purple.purity.dt,
		file = paste(run.dir, "/result_summaries/", gpl_prefix, "/purple.purity.summary.txt", sep=""),
		col.names = TRUE,
		sep = '\t')

save(list=c('purple.qc.dt', 'purple.purity.dt'),
	file = paste(run.dir, "/result_summaries/", gpl_prefix, "/QC.info.Rdata", sep=""))


#load(file = paste(run.dir, "/result_summaries/", gpl_prefix, "/QC.info.Rdata", sep=""))

purple.qc.dt.failed.contamination <- purple.qc.dt[grep("FAIL_CONTAMINATION", QCStatus, fixed=TRUE),]

FAIL_CONTAMINATION.indexes <- grep("FAIL_CONTAMINATION", purple.qc.dt$QCStatus, fixed=TRUE)
FAIL_NO_TUMOR.indexes <- grep("FAIL_NO_TUMOR", purple.qc.dt$QCStatus, fixed=TRUE)

fwrite(purple.qc.dt[which(!1:nrow(purple.qc.dt)%in%c(FAIL_CONTAMINATION.indexes, FAIL_NO_TUMOR.indexes)),
	list(subject.id, tumor.sample.id, normal.sample.id)],
	file = paste(run.dir, "/result_summaries/", gpl_prefix, "/vep_annotation_sample_info.tsv", sep=""),
	col.names = TRUE,
	sep = '\t')

message( paste("Importing germline driver tables for ", study.name, sep=""))

tumor.sample.ids.for.drivers.ls <- tumor.sample.ids.ls[which(!tumor.sample.ids.ls%in%purple.qc.dt.failed.contamination$tumor.sample.id)]

germline.drivers.dt.ls <- lapply(
	X = tumor.sample.ids.for.drivers.ls,
	FUN = function( curr.sample.id ){
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
		curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
		curr.dt <- fread( paste(germline.dir, "/", curr.subject.id, "/purple/", curr.sample.id, "/", curr.sample.id, ".driver.catalog.germline.tsv", sep=""),
			header = TRUE)
		curr.dt[,subject.id:=curr.subject.id]
		curr.dt[,tumor.sample.id:=curr.sample.id]
		curr.dt[,subject.index:=curr.subject.index]
		curr.dt[,tumor.index:=curr.tumor.index]
		curr.dt
	})

germline.drivers.dt <- rbindlist(germline.drivers.dt.ls)


message( paste("Merged germline driver data for ", study.name, "\n", sep=""))


message( paste("Importing somatic driver tables for ", study.name, sep=""))

somatic.drivers.dt.ls <- lapply(
	X = tumor.sample.ids.for.drivers.ls,
	FUN = function( curr.sample.id ){
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
		curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
		curr.dt <- fread( paste(germline.dir, "/", curr.subject.id, "/purple/", curr.sample.id, "/", curr.sample.id, ".driver.catalog.somatic.tsv", sep=""),
		header = TRUE)
		curr.dt[,subject.id:=curr.subject.id]
		curr.dt[,tumor.sample.id:=curr.sample.id]
		curr.dt[,subject.index:=curr.subject.index]
		curr.dt[,tumor.index:=curr.tumor.index]
		curr.dt
	})

somatic.drivers.dt <- rbindlist(somatic.drivers.dt.ls)

message( paste("Merged somatic driver data for ", study.name, "\n", sep=""))

chrom.map.ls <- as.integer(1:25)
names(chrom.map.ls) <- paste("chr", c(1:22, "X", "Y", "M"), sep="")

samples.with.germline.drivers <- unique(germline.drivers.dt$tumor.sample.id)
names(samples.with.germline.drivers) <- samples.with.germline.drivers

parse_ALT <- function( curr.ALT ){
	paste(as.character(curr.ALT[[1]]), collapse=",")
}

message( paste("Loading somatic variants for ", study.name, sep=""))
load( file = paste(run.dir, "/result_summaries/", gpl_prefix, "/somatic.driver.variant.info.Rdata", sep=""))

message( paste("Loading germline variants for ", study.name, sep=""))
if (nrow(germline.drivers.dt)>0 ){
	load( file = paste(run.dir, "/result_summaries/", gpl_prefix, "/germline.variant.info.Rdata", sep=""))

	germline.driver.variant.info.dt <- germline.variant.info.dt[REPORTED==TRUE &
		!((AF_TOT>=maf.cutoff & AF_TOT<af.cutoff) | AF_TOT>=af.cutoff | 
		(AF_EAS>=maf.cutoff & AF_EAS<af.cutoff) | AF_EAS>=af.cutoff | 
		(TOMMO>=maf.cutoff & TOMMO<af.cutoff) | TOMMO>=af.cutoff |
		(KOREAN>=maf.cutoff & KOREAN<af.cutoff) | KOREAN>=af.cutoff |
		(dbGaP_PopFreq>=maf.cutoff & dbGaP_PopFreq<af.cutoff) | dbGaP_PopFreq>=af.cutoff),]

	germline.drivers.dt	<- germline.drivers.dt[gene%in%germline.driver.variant.info.dt$gene,]
	
	rm(germline.variant.info.dt)
	gc()
}else{
	message( paste("There were no germline drivers. Merged germline variants for ", study.name, " not exported.\n", sep=""))
	germline.driver.variant.info.dt <- NULL;
}


## Keeping the same sample list with samples removed that are contaminated
tumor.sample.ids.ls <- copy(tumor.sample.ids.for.drivers.ls)

message( paste("Importing GRIDSS SV variants for ", study.name, sep=""))

gridss.sv.variant.info.ls <- mclapply(
	X = tumor.sample.ids.ls[1],
	FUN = function( curr.sample.id ){
##			curr.sample.id <- tumor.sample.ids.ls[[1]]
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
		curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
		
		file.gz <- paste(purple.dir, "/", curr.subject.id, "/gripss/", curr.sample.id, ".gripss.somatic.filtered.vcf.gz", sep="")
		
		if( file.exists(file.gz)!=TRUE ){
			file.gz <- paste(purple.dir, "/", curr.subject.id, "/gripss/", curr.sample.id, ".gripss.filtered.vcf.gz", sep="")
		}
		vcf <- readVcf(file.gz, "GRCh38")

		vcf.header.dt <- as.data.table(as.data.frame(rowRanges(vcf)), keep.rownames = "ID")
		vcf.header.dt[,CHROM:=as.character(seqnames)]
		vcf.header.dt[,chr:=chrom.map.ls[CHROM]]
		vcf.header.dt[,chr:=ifelse(is.na(chr), 99L, chr)];
		vcf.header.dt[,ALT:=parse_ALT(curr.ALT=ALT), by=1:nrow(vcf.header.dt)]
		vcf.header.dt$ALT <- unlist(vcf.header.dt$ALT)
		vcf.header.dt[,subject.id:=curr.subject.id]
		vcf.header.dt[,tumor.sample.id:=curr.sample.id]
		vcf.header.dt[,subject.index:=curr.subject.index]
		vcf.header.dt[,tumor.index:=curr.tumor.index]
		
		vcf.dt <- as.data.table(info(vcf))

		cbind(vcf.header.dt[,list(subject.id, tumor.sample.id, chr, CHROM, POS=start, REF, ALT, ID, QUAL, FILTER)], vcf.dt)
	},
	mc.cores = threads)

gridss.sv.variant.info.dt <- rbindlist(gridss.sv.variant.info.ls)

message( paste("Merged GRIDSS SV variants for ", study.name, "\n", sep=""))


message( paste("Importing PURPLE SV variants for ", study.name, sep=""))

# Purple SV variant extraction
purple.sv.variant.info.ls <- mclapply(
	X = tumor.sample.ids.ls,
	FUN = function( curr.sample.id ){
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
		curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
		
		sv.destination.file <- tempfile()
		sv.file.gz <- paste(purple.dir, "/", curr.subject.id, "/purple/", curr.sample.id, "/", curr.sample.id, ".purple.sv.vcf.gz", sep="")
		sv.tabix.file <- TabixFile(sv.file.gz, yieldSize=100000)		
		sv.vcf <- readVcf(sv.file.gz, "GRCh38")
		
		sv.vcf.header.dt <- as.data.table(as.data.frame(rowRanges(sv.vcf)), keep.rownames = "ID")
		sv.vcf.header.dt[,CHROM:=as.character(seqnames)]
		sv.vcf.header.dt[,chr:=chrom.map.ls[CHROM]]
		sv.vcf.header.dt[,chr:=ifelse(is.na(chr), 99L, chr)];
		sv.vcf.header.dt[,ALT:=parse_ALT(curr.ALT=ALT), by=1:nrow(sv.vcf.header.dt)]
		sv.vcf.header.dt$ALT <- unlist(sv.vcf.header.dt$ALT)
		sv.vcf.header.dt[,subject.id:=curr.subject.id]
		sv.vcf.header.dt[,tumor.sample.id:=curr.sample.id]
		sv.vcf.header.dt[,subject.index:=curr.subject.index]
		sv.vcf.header.dt[,tumor.index:=curr.tumor.index]
		
		sv.vcf$info
		sv.vcf.dt <- as.data.table(info(sv.vcf))
		
		cbind(sv.vcf.header.dt[,list(subject.id, tumor.sample.id, chr, CHROM, POS=start, REF, ALT, ID, QUAL, FILTER)], sv.vcf.dt)
	},
	mc.cores = threads)

purple.sv.variant.info.dt <- rbindlist(purple.sv.variant.info.ls, fill=TRUE)

message( paste("Merged PURPLE SV variants for ", study.name, "\n", sep=""))


message( paste("Importing PURPLE CNV gene info for ", study.name, sep=""))

## Purple CNV and segment extraction
purple.cnv.gene.info.ls <- mclapply(
	X = tumor.sample.ids.ls,
	FUN = function( curr.sample.id ){
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
		curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
		
		curr.dt <- fread( paste(purple.dir, "/", curr.subject.id, "/purple/", curr.sample.id, "/", curr.sample.id, ".purple.cnv.gene.tsv", sep=""))
		curr.dt[,subject.id:=curr.subject.id]
		curr.dt[,tumor.sample.id:=curr.sample.id]
		curr.dt[,subject.index:=curr.subject.index]
		curr.dt[,tumor.index:=curr.tumor.index]
		
		curr.dt
	},
	mc.cores = threads)

purple.cnv.gene.info <- rbindlist(purple.cnv.gene.info.ls)

message( paste("Merged PURPLE CNV gene info for ", study.name, "\n", sep=""))


message( paste("Importing PURPLE CNV germline info for ", study.name, sep=""))

purple.cnv.germline.info.ls <- mclapply(
	X = tumor.sample.ids.ls,
	FUN = function( curr.sample.id ){
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
		curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
		
		curr.dt <- fread( paste(purple.dir, "/", curr.subject.id, "/purple/", curr.sample.id, "/", curr.sample.id, ".purple.cnv.germline.tsv", sep=""))
		curr.dt[,subject.id:=curr.subject.id]
		curr.dt[,tumor.sample.id:=curr.sample.id]
		curr.dt[,subject.index:=curr.subject.index]
		curr.dt[,tumor.index:=curr.tumor.index]
		
		curr.dt
	},
	mc.cores = threads)

purple.cnv.germline.info <- rbindlist(purple.cnv.germline.info.ls)

message( paste("Merged PURPLE CNV genrmline info for ", study.name, "\n", sep=""))


message( paste("Importing PURPLE CNV somatic info for ", study.name, sep=""))

purple.cnv.somatic.info.ls <- mclapply(
	X = tumor.sample.ids.ls,
	FUN = function( curr.sample.id ){
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
		curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
		
		curr.dt <- fread( paste(purple.dir, "/", curr.subject.id, "/purple/", curr.sample.id, "/", curr.sample.id, ".purple.cnv.somatic.tsv", sep=""))
		curr.dt[,subject.id:=curr.subject.id]
		curr.dt[,tumor.sample.id:=curr.sample.id]
		curr.dt[,subject.index:=curr.subject.index]
		curr.dt[,tumor.index:=curr.tumor.index]
		
		curr.dt
	},
	mc.cores = threads)

purple.cnv.somatic.info <- rbindlist(purple.cnv.somatic.info.ls)

message( paste("Merged PURPLE CNV somatic info for ", study.name, "\n", sep=""))


message( paste("Importing PURPLE segment info for ", study.name, sep=""))

purple.segment.info.ls <- mclapply(
	X = tumor.sample.ids.ls,
	FUN = function( curr.sample.id ){
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
		curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
		
		curr.dt <- fread( paste(purple.dir, "/", curr.subject.id, "/purple/", curr.sample.id, "/", curr.sample.id, ".purple.segment.tsv", sep=""))
		curr.dt[,subject.id:=curr.subject.id]
		curr.dt[,tumor.sample.id:=curr.sample.id]
		curr.dt[,subject.index:=curr.subject.index]
		curr.dt[,tumor.index:=curr.tumor.index]
		
		curr.dt
	},
	mc.cores = threads)

purple.segment.info <- rbindlist(purple.segment.info.ls)

message( paste("Merged PURPLE segment info for ", study.name, "\n", sep=""))


tumor.sample.ids.ls <- tumor.normal.pair.info[which(sample.id.N!='none'),]$sample.id.T

message( paste("Importing Linx breakend info for ", study.name, sep=""))

## Linx output extraction
linx.breakend.dt.ls <- lapply(
	X = tumor.sample.ids.ls,
	FUN = function( curr.sample.id ){
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
		curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
		
		curr.dt <- fread( paste(germline.dir, "/", curr.subject.id, "/linx/", curr.sample.id, "/", curr.sample.id, ".linx.breakend.tsv", sep=""))	
		curr.dt[,subject.id:=curr.subject.id]
		curr.dt[,tumor.sample.id:=curr.sample.id]
		curr.dt[,subject.index:=curr.subject.index]
		curr.dt[,tumor.index:=curr.tumor.index]
		
		curr.dt
	})
linx.breakend.dt <- rbindlist(linx.breakend.dt.ls)

message( paste("Merged Linx breakend info for ", study.name, "\n", sep=""))


message( paste("Importing Linx clusters for ", study.name, sep=""))

linx.clusters.dt.ls <- lapply(
	X = tumor.sample.ids.ls,
	FUN = function( curr.sample.id ){
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
		curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
		
		curr.dt <- fread( paste(germline.dir, "/", curr.subject.id, "/linx/", curr.sample.id, "/", curr.sample.id, ".linx.clusters.tsv", sep=""))	
		curr.dt[,subject.id:=curr.subject.id]
		curr.dt[,tumor.sample.id:=curr.sample.id]
		curr.dt[,subject.index:=curr.subject.index]
		curr.dt[,tumor.index:=curr.tumor.index]
		
		curr.dt
	})
linx.clusters.dt <- rbindlist(linx.clusters.dt.ls)

message( paste("Merged Linx clusters for ", study.name, "\n", sep=""))


message( paste("Importing Linx drivers for ", study.name, sep=""))

linx.driver.catalog.dt.ls <- lapply(
		X = tumor.sample.ids.ls,
		FUN = function( curr.sample.id ){
			curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
			curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
			curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
			curr.dt <- fread( paste(germline.dir, "/", curr.subject.id, "/linx/", curr.sample.id, "/", curr.sample.id, ".linx.driver.catalog.tsv", sep=""), header=TRUE)	
			
			if( nrow(curr.dt)>0 ){
				curr.dt[,subject.id:=curr.subject.id]
				curr.dt[,tumor.sample.id:=curr.sample.id]
				curr.dt[,subject.index:=curr.subject.index]
				curr.dt[,tumor.index:=curr.tumor.index]
			# may need to change when integrating current LINX 1.17 or 1.18
			# BHD had mix of 1.17 beta and released 1.17; columns slightly differed
				curr.dt[,list(subject.id, tumor.sample.id, subject.index, tumor.index,
					chromosome, chromosomeBand, gene, driver, category, likelihoodMethod, driverLikelihood, 
					missense, nonsense, splice, inframe, frameshift, biallelic, minCopyNumber, maxCopyNumber)]
			}else{
				NULL
			}
		})
linx.driver.catalog.dt <- rbindlist(linx.driver.catalog.dt.ls)

linx.drivers.dt.ls <- lapply(
	X = tumor.sample.ids.ls,
	FUN = function( curr.sample.id ){
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
		curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
		curr.dt <- fread( paste(germline.dir, "/", curr.subject.id, "/linx/", curr.sample.id, "/", curr.sample.id, ".linx.drivers.tsv", sep=""))	
		curr.dt[,subject.id:=curr.subject.id]
		curr.dt[,tumor.sample.id:=curr.sample.id]
		curr.dt[,subject.index:=curr.subject.index]
		curr.dt[,tumor.index:=curr.tumor.index]
		
		curr.dt
	})
linx.drivers.dt <- rbindlist(linx.drivers.dt.ls)

message( paste("Merged Linx drivers for ", study.name, "\n", sep=""))

if (nrow(germline.drivers.dt)>0 ){
	merged.drivers.dt <- rbind(
			germline.drivers.dt[,list(subject.id, tumor.sample.id, subject.index, tumor.index, chromosome, chromosomeBand, gene, driver, category, likelihoodMethod, driverLikelihood, missense, nonsense, splice, inframe, frameshift, biallelic, minCopyNumber, maxCopyNumber)],
			somatic.drivers.dt[,list(subject.id, tumor.sample.id, subject.index, tumor.index, chromosome, chromosomeBand, gene, driver, category, likelihoodMethod, driverLikelihood, missense, nonsense, splice, inframe, frameshift, biallelic, minCopyNumber, maxCopyNumber)],
			linx.driver.catalog.dt[,list(subject.id, tumor.sample.id, subject.index, tumor.index, chromosome, chromosomeBand, gene, driver, category, likelihoodMethod, driverLikelihood, missense, nonsense, splice, inframe, frameshift, biallelic, minCopyNumber, maxCopyNumber)])
}else{
	merged.drivers.dt <- rbind(
			somatic.drivers.dt[,list(subject.id, tumor.sample.id, subject.index, tumor.index, chromosome, chromosomeBand, gene, driver, category, likelihoodMethod, driverLikelihood, missense, nonsense, splice, inframe, frameshift, biallelic, minCopyNumber, maxCopyNumber)],
			linx.driver.catalog.dt[,list(subject.id, tumor.sample.id, subject.index, tumor.index, chromosome, chromosomeBand, gene, driver, category, likelihoodMethod, driverLikelihood, missense, nonsense, splice, inframe, frameshift, biallelic, minCopyNumber, maxCopyNumber)])
}

merged.drivers.dt <- unique(merged.drivers.dt[order(subject.index, tumor.index,  gene),])

message( paste("Saving merged driver data for ", study.name, "\n", sep=""))

save(list=c('germline.drivers.dt', 'somatic.drivers.dt', 'merged.drivers.dt'),
		file = paste(run.dir, "/result_summaries/", gpl_prefix, "/driver.catalog.ver.20220210.Rdata", sep=""))


message( paste("Importing Linx fusions for ", study.name, sep=""))

linx.fusion.dt.ls <- lapply(
	X = tumor.sample.ids.ls,
	FUN = function( curr.sample.id ){
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
		curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
		
		curr.dt <- fread( paste(germline.dir, "/", curr.subject.id, "/linx/", curr.sample.id, "/", curr.sample.id, ".linx.fusion.tsv", sep=""))	
		curr.dt[,subject.id:=curr.subject.id]
		curr.dt[,tumor.sample.id:=curr.sample.id]
		curr.dt[,subject.index:=curr.subject.index]
		curr.dt[,tumor.index:=curr.tumor.index]
		
		curr.dt
	})
linx.fusion.dt <- rbindlist(linx.fusion.dt.ls)

message( paste("Merged Linx fusions for ", study.name, "\n", sep=""))


message( paste("Importing Linx links for ", study.name, sep=""))

linx.links.dt.ls <- lapply(
	X = tumor.sample.ids.ls,
	FUN = function( curr.sample.id ){
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
		curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
		
		curr.dt <- fread( paste(germline.dir, "/", curr.subject.id, "/linx/", curr.sample.id, "/", curr.sample.id, ".linx.links.tsv", sep=""))	
		curr.dt[,subject.id:=curr.subject.id]
		curr.dt[,tumor.sample.id:=curr.sample.id]
		curr.dt[,subject.index:=curr.subject.index]
		curr.dt[,tumor.index:=curr.tumor.index]
		if('lowerBreakendId'%in%colnames(curr.dt)){
			setnames(curr.dt,
				old=c('lowerBreakendId', 'upperBreakendId'),
				new = c('lowerSvId', 'upperSvId'))
		}
		curr.dt
	})
linx.links.dt <- rbindlist(linx.links.dt.ls)

message( paste("Merged Linx links for ", study.name, "\n", sep=""))


message( paste("Importing Linx SVs for ", study.name, sep=""))

linx.svs.dt.ls <- lapply(
	X = tumor.sample.ids.ls,
	FUN = function( curr.sample.id ){
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
		curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
		
		curr.dt <- fread( paste(germline.dir, "/", curr.subject.id, "/linx/", curr.sample.id, "/", curr.sample.id, ".linx.svs.tsv", sep=""))	
		curr.dt[,subject.id:=curr.subject.id]
		curr.dt[,tumor.sample.id:=curr.sample.id]
		curr.dt[,subject.index:=curr.subject.index]
		curr.dt[,tumor.index:=curr.tumor.index]
		
		curr.dt[,list(subject.id, tumor.sample.id, subject.index, tumor.index,
			vcfId, svId, clusterId, clusterReason, fragileSiteStart, fragileSiteEnd, isFoldback,
			lineTypeStart, lineTypeEnd, junctionCopyNumberMin, junctionCopyNumberMax,
			geneStart, geneEnd,
			localTopologyIdStart, localTopologyIdEnd, localTopologyStart, localTopologyEnd, localTICountStart, localTICountEnd)]
	})
linx.svs.dt <- rbindlist(linx.svs.dt.ls)

message( paste("Merged Linx SVs for ", study.name, "\n", sep=""))


## message( paste("Importing Linx viral inserts for ", study.name, sep=""))
## 
## linx.viral_inserts.dt.ls <- lapply(
##     X = tumor.sample.ids.ls,
##     FUN = function( curr.sample.id ){
##         curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
##         curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
##         curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
## 
##         curr.dt <- fread( paste(germline.dir, "/", curr.subject.id, "/linx/", curr.sample.id, "/", curr.sample.id, ".linx.viral_inserts.tsv", sep=""))	
##         curr.dt[,subject.id:=curr.subject.id]
##         curr.dt[,tumor.sample.id:=curr.sample.id]
##         curr.dt[,subject.index:=curr.subject.index]
##         curr.dt[,tumor.index:=curr.tumor.index]
## 
##         curr.dt
##     })
## linx.viral_inserts.dt <- rbindlist(linx.viral_inserts.dt.ls)
## 
## message( paste("Merged Linx viral inserts for ", study.name, "\n", sep=""))


message( paste("Importing Linx vis-copy-number for ", study.name, sep=""))

linx.vis_copy_number.info.ls <- mclapply(
	X = tumor.sample.ids.ls,
	FUN = function( curr.sample.id ){
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
		curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
		
		cn.dt <- fread(paste(purple.dir, "/", curr.subject.id, "/linx/", curr.sample.id, "/", curr.sample.id, ".linx.vis_copy_number.tsv", sep=""))
		cn.dt[,subject.id:=curr.subject.id]
		cn.dt[,tumor.sample.id:=curr.sample.id]
		cn.dt[,subject.index:=curr.subject.index]
		cn.dt[,tumor.index:=curr.tumor.index]
		
		cn.dt
	},
	mc.cores = 8)

linx.vis_copy_number.info.dt <- rbindlist(linx.vis_copy_number.info.ls)

message( paste("Merged Linx vis-copy-number for ", study.name, "\n", sep=""))


message( paste("Importing Linx vis-fusion for ", study.name, sep=""))

linx.vis_fusion.info.ls <- mclapply(
	X = tumor.sample.ids.ls,
	FUN = function( curr.sample.id ){
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
		curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
		
		cn.dt <- fread(paste(purple.dir, "/", curr.subject.id, "/linx/", curr.sample.id, "/", curr.sample.id, ".linx.vis_fusion.tsv", sep=""))
		cn.dt[,subject.id:=curr.subject.id]
		cn.dt[,tumor.sample.id:=curr.sample.id]
		cn.dt[,subject.index:=curr.subject.index]
		cn.dt[,tumor.index:=curr.tumor.index]
		
		cn.dt
	},
	mc.cores = threads)

linx.vis_fusion.info.dt <- rbindlist(linx.vis_fusion.info.ls)

message( paste("Merged Linx vis-fusion for ", study.name, "\n", sep=""))


message( paste("Importing Linx vis-gene exon info for ", study.name, sep=""))

linx.vis_gene_exon.info.ls <- mclapply(
	X = tumor.sample.ids.ls,
	FUN = function( curr.sample.id ){
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
		curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
		
		cn.dt <- fread(paste(purple.dir, "/", curr.subject.id, "/linx/", curr.sample.id, "/", curr.sample.id, ".linx.vis_gene_exon.tsv", sep=""))
		cn.dt[,subject.id:=curr.subject.id]
		cn.dt[,tumor.sample.id:=curr.sample.id]
		cn.dt[,subject.index:=curr.subject.index]
		cn.dt[,tumor.index:=curr.tumor.index]
		
		cn.dt
	},
	mc.cores = threads)

linx.vis_gene_exon.info.dt <- rbindlist(linx.vis_gene_exon.info.ls)

message( paste("Merged Linx vis-gene exon info for ", study.name, "\n", sep=""))


message( paste("Importing Linx vis-protein domain info for ", study.name, sep=""))

linx.vis_protein_domain.info.ls <- mclapply(
	X = tumor.sample.ids.ls,
	FUN = function( curr.sample.id ){
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
		curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
		
		cn.dt <- fread(paste(purple.dir, "/", curr.subject.id, "/linx/", curr.sample.id, "/", curr.sample.id, ".linx.vis_protein_domain.tsv", sep=""))
		cn.dt[,subject.id:=curr.subject.id]
		cn.dt[,tumor.sample.id:=curr.sample.id]
		cn.dt[,subject.index:=curr.subject.index]
		cn.dt[,tumor.index:=curr.tumor.index]
		
		cn.dt
	},
	mc.cores = threads)

linx.vis_protein_domain.info.dt <- rbindlist(linx.vis_protein_domain.info.ls)

message( paste("Merged Linx vis-protein domain info for ", study.name, "\n", sep=""))


message( paste("Importing Linx vis-segments info for ", study.name, sep=""))

linx.vis_segments.info.ls <- mclapply(
	X = tumor.sample.ids.ls,
	FUN = function( curr.sample.id ){
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
		curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
		
		cn.dt <- fread(paste(purple.dir, "/", curr.subject.id, "/linx/", curr.sample.id, "/", curr.sample.id, ".linx.vis_segments.tsv", sep=""))
		cn.dt[,subject.id:=curr.subject.id]
		cn.dt[,tumor.sample.id:=curr.sample.id]
		cn.dt[,subject.index:=curr.subject.index]
		cn.dt[,tumor.index:=curr.tumor.index]
		
		cn.dt
	},
	mc.cores = threads)

linx.vis_segments.info.dt <- rbindlist(linx.vis_segments.info.ls)

message( paste("Merged Linx vis-segments info for ", study.name, "\n", sep=""))


message( paste("Importing Linx vis-SV data info for ", study.name, sep=""))

linx.vis_sv_data.info.ls <- mclapply(
	X = tumor.sample.ids.ls,
	FUN = function( curr.sample.id ){
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
		curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
		
		cn.dt <- fread(paste(purple.dir, "/", curr.subject.id, "/linx/", curr.sample.id, "/", curr.sample.id, ".linx.vis_sv_data.tsv", sep=""))
		cn.dt[,subject.id:=curr.subject.id]
		cn.dt[,tumor.sample.id:=curr.sample.id]
		cn.dt[,subject.index:=curr.subject.index]
		cn.dt[,tumor.index:=curr.tumor.index]
		
		cn.dt
	},
	mc.cores = threads)

linx.vis_sv_data.info.dt <- rbindlist(linx.vis_sv_data.info.ls)

if (import.isofox.rna.fusion.match==TRUE){
	linx.isofox.match.info.ls <- mclapply(
		X = tumor.sample.ids.ls,
		FUN = function( curr.sample.id ){
			curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
			curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
			curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
			
			curr.dt <- fread(paste(purple.dir, "/", curr.subject.id, "/linx/", curr.sample.id, "/LNX_RNA_FUSION_MATCH_ISOFOX.csv", sep=""))
			curr.cols <- copy(colnames(curr.dt))
			
			curr.dt[,subject.id:=curr.subject.id]
			curr.dt[,subject.index:=curr.subject.index]
			curr.dt[,tumor.index:=curr.tumor.index]
			
			curr.dt[,c('subject.id', 'subject.index', 'tumor.index', curr.cols), with=FALSE]
		},
		mc.cores = threads)
	linx.isofox.match.info.dt <- rbindlist(linx.isofox.match.info.ls)
}else{
	linx.isofox.match.info.dt <- NULL
}


linx.neo_epitope.info.ls <- mclapply(
	X = tumor.sample.ids.ls,
	FUN = function( curr.sample.id ){
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
		curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
		
		if (file.exists(paste(purple.dir, "/", curr.subject.id, "/linx/", curr.sample.id, "/", curr.sample.id, ".linx.neo_epitope.tsv", sep=""))){
			curr.dt <- fread(paste(purple.dir, "/", curr.subject.id, "/linx/", curr.sample.id, "/", curr.sample.id, ".linx.neo_epitope.tsv", sep=""), sep=tsv.sep)
			process.dt <- TRUE
		}else if (file.exists(paste(purple.dir, "/", curr.subject.id, "/linx/", curr.sample.id, "/", curr.sample.id, ".linx.neo_epitope.csv", sep=""))){
			curr.dt <- fread(paste(purple.dir, "/", curr.subject.id, "/linx/", curr.sample.id, "/", curr.sample.id, ".linx.neo_epitope.csv", sep=""), sep=csv.sep)
			process.dt <- TRUE
		}else{
			process.dt <- FALSE
		}
		
		if (process.dt==TRUE){
			curr.cols <- copy(colnames(curr.dt))
			
			curr.dt[,subject.id:=curr.subject.id]
			curr.dt[,tumor.sample.id:=curr.sample.id]
			curr.dt[,subject.index:=curr.subject.index]
			curr.dt[,tumor.index:=curr.tumor.index]
			
			curr.dt[,c('subject.id', 'subject.index', 'tumor.index', 'tumor.sample.id', curr.cols), with=FALSE]
		}else{
			NULL
		}
	},
	mc.cores = threads)

linx.neo_epitope.info.dt <- rbindlist(linx.neo_epitope.info.ls)


message( paste("Merged Linx vis-SV data info for ", study.name, "\n", sep=""))


save(list=c('purple.qc.dt', 'purple.purity.dt',
	'merged.drivers.dt', 'germline.driver.variant.info.dt', 'somatic.driver.variant.info.dt',
	'gridss.sv.variant.info.dt', 'purple.sv.variant.info.dt',
	'purple.cnv.gene.info', 'purple.cnv.germline.info', 'purple.cnv.somatic.info', 'purple.segment.info',
	'linx.breakend.dt', 'linx.clusters.dt', 'linx.drivers.dt', 'linx.fusion.dt', 'linx.links.dt', 'linx.svs.dt',
	'linx.driver.catalog.dt',
#	'linx.viral_inserts.dt',
	'linx.vis_copy_number.info.dt', 'linx.vis_fusion.info.dt', 'linx.vis_gene_exon.info.dt', 'linx.vis_protein_domain.info.dt', 'linx.vis_segments.info.dt', 'linx.vis_sv_data.info.dt',
	'linx.neo_epitope.info.dt', 'linx.isofox.match.info.dt'),
	file = paste(run.dir, "/result_summaries/", gpl_prefix, "/driver.catalog.germline.somatic.variant.SV.info.ver.20220210.Rdata", sep=""))

message( paste("Saved driver and variant data for ", study.name, sep=""))