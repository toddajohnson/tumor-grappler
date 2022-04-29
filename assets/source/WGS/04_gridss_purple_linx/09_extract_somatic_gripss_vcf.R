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

message( paste("Loading QC tables for ", study.name, sep=""))

load( file = paste(run.dir, "/result_summaries/", gpl_prefix, "/QC.info.Rdata", sep=""))

purple.qc.dt.failed.contamination <- purple.qc.dt[grep("FAIL_CONTAMINATION", QCStatus, fixed=TRUE),]

## Keeping the same sample list with samples removed that are contaminated
tumor.sample.ids.ls <- tumor.normal.pair.info[!sample.id.T%in%purple.qc.dt.failed.contamination$tumor.sample.id,]$sample.id.T
names(tumor.sample.ids.ls) <- tumor.sample.ids.ls
sample.ct <- length(tumor.sample.ids.ls)

threads <- min( c(18, sample.ct) )

chrom.map.ls <- as.integer(1:25)
names(chrom.map.ls) <- paste("chr", c(1:22, "X", "Y", "M"), sep="")

parse_ALT <- function( curr.ALT ){
	paste(as.character(curr.ALT[[1]]), collapse=",")
}

message( paste("Importing GRDISS SV variants for ", study.name, sep=""))

gridss.somatic.sv.variant.info.ls <- mclapply(
	X = tumor.sample.ids.ls,
	FUN = function( curr.sample.id ){
##			curr.sample.id <- tumor.sample.ids.ls[[9]]
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		curr.subject.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.index
		curr.tumor.index <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$tumor.index
		
		sv.destination.file <- tempfile()
		sv.file.gz <- paste(purple.dir, "/", curr.subject.id, "/gripss/", curr.sample.id, ".gripss.somatic.vcf.gz", sep="")
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
		
		sv.vcf.dt <- as.data.table(info(sv.vcf))

		cbind(sv.vcf.header.dt[,list(subject.id, tumor.sample.id, chr, CHROM, POS=start, REF, ALT, ID, QUAL, FILTER)], sv.vcf.dt)
	},
	mc.cores = threads)

gridss.somatic.sv.variant.info.dt <- rbindlist(gridss.somatic.sv.variant.info.ls)

message( paste("Merged GRIDSS SV variants for ", study.name, "\n", sep=""))

save(list=c('gridss.somatic.sv.variant.info.dt'),
	file = paste(run.dir, "/result_summaries/", gpl_prefix, "/somatic.variant.SV.info.Rdata", sep=""))

message( paste("Saved somatic gripss SV data for ", study.name, sep=""))