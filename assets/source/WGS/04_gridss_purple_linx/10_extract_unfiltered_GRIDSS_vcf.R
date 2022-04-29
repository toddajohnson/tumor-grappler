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
subject.ids.ls <- unique(tumor.normal.pair.info[!sample.id.T%in%purple.qc.dt.failed.contamination$tumor.sample.id,]$subject.id)
names(subject.ids.ls) <- subject.ids.ls
sample.ct <- length(subject.ids.ls)

threads <- min( c(18, sample.ct) )

chrom.map.ls <- as.integer(1:25)
names(chrom.map.ls) <- paste("chr", c(1:22, "X", "Y", "M"), sep="")

parse_ALT <- function( curr.ALT ){
	paste(as.character(curr.ALT[[1]]), collapse=",")
}

message( paste("Importing GRIDSS SV variants for ", study.name, sep=""))

gridss.unfiltered.sv.variant.info.ls <- mclapply(
	X = subject.ids.ls,
	FUN = function( curr.subject.id ){
		curr.tumor.sample.id <- tumor.normal.pair.info[subject.id==curr.subject.id,]$sample.id.T
		curr.normal.sample.id <- tumor.normal.pair.info[subject.id==curr.subject.id,]$sample.id.N

		sv.destination.file <- tempfile()
		sv.file.gz <- paste(purple.dir, "/", curr.subject.id, "/gridss/", curr.subject.id, ".gridss.unfiltered.vcf.gz", sep="")
		sv.tabix.file <- TabixFile(sv.file.gz, yieldSize=100000)		
		sv.vcf <- readVcf(sv.file.gz, "GRCh38")

		sv.vcf.header.dt <- as.data.table(as.data.frame(rowRanges(sv.vcf)), keep.rownames = "ID")
		sv.vcf.header.dt[,CHROM:=as.character(seqnames)]
		sv.vcf.header.dt[,chr:=chrom.map.ls[CHROM]]
		sv.vcf.header.dt[,chr:=ifelse(is.na(chr), 99L, chr)];
		sv.vcf.header.dt[,ALT:=parse_ALT(curr.ALT=ALT), by=1:nrow(sv.vcf.header.dt)]
		sv.vcf.header.dt$ALT <- unlist(sv.vcf.header.dt$ALT)
		sv.vcf.header.dt[,subject.id:=curr.subject.id]
		
		sv.vcf.BVF <- as.data.table(geno(sv.vcf)$BVF)
		setnames(sv.vcf.BVF, old=c(curr.normal.sample.id, curr.tumor.sample.id), new=c('BVF.N', 'BVF.T'))
		
		sv.vcf.VF <- as.data.table(geno(sv.vcf)$VF)
		setnames(sv.vcf.VF, old=c(curr.normal.sample.id, curr.tumor.sample.id), new=c('VF.N', 'VF.T'))

		sv.vcf.REFREAD <- as.data.table(geno(sv.vcf)$REF)
		setnames(sv.vcf.REFREAD, old=c(curr.normal.sample.id, curr.tumor.sample.id), new=c('REFREAD.N', 'REFREAD.T'))
		
		sv.vcf.REFPAIR <- as.data.table(geno(sv.vcf)$REFPAIR)
		setnames(sv.vcf.REFPAIR, old=c(curr.normal.sample.id, curr.tumor.sample.id), new=c('REFPAIR.N', 'REFPAIR.T'))
		
		sv.vcf.BQ <- as.data.table(geno(sv.vcf)$BQ)
		setnames(sv.vcf.BQ, old=c(curr.normal.sample.id, curr.tumor.sample.id), new=c('BQ.N', 'BQ.T'))
		
		sv.vcf.dt <- as.data.table(info(sv.vcf))
		setnames(sv.vcf.dt, old="REF", new="REFREAD")
		
		cbind(sv.vcf.header.dt[,list(subject.id, chr, CHROM, POS=start, REF, ALT, ID, QUAL, FILTER)], sv.vcf.dt, sv.vcf.BVF, sv.vcf.VF, sv.vcf.REFREAD, sv.vcf.REFPAIR, sv.vcf.BQ)
	},
	mc.cores = threads)

gridss.unfiltered.sv.variant.info.dt <- rbindlist(gridss.unfiltered.sv.variant.info.ls)

message( paste("Merged GRIDSS SV variants for ", study.name, "\n", sep=""))

save(list=c('gridss.unfiltered.sv.variant.info.dt'),
	file = paste(run.dir, "/result_summaries/", gpl_prefix, "/gridss.unfiltered.variant.SV.info.Rdata", sep=""))

message( paste("Saved unfiltered gridss SV data for ", study.name, sep=""))