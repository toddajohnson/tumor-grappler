#!/usr/bin/env Rscript

# module use /usr/local/package/modulefiles/
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
	library(parallel)})

purple.dir <- paste(run.dir, "/result/", gpl_prefix, sep="")

if ( length(list.dirs(paste(run.dir, "/result_summaries/", gpl_prefix, sep=""), recursive=FALSE))==0 ){
	system(paste("mkdir ", run.dir, "/result_summaries/", gpl_prefix, sep=""))
}

load( file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))

tumor.normal.pair.info <- tumor_normal_pair_info.GPL[which(sample.id.N!='none',)]

tumor.sample.ids.ls <- tumor.normal.pair.info$sample.id.T
names(tumor.sample.ids.ls) <- tumor.sample.ids.ls
sample.ct <- length(tumor.sample.ids.ls)

threads <- min( c(6, sample.ct) )

chrom.map.ls <- as.integer(1:25)
names(chrom.map.ls) <- paste("chr", c(1:22, "X", "Y", "M"), sep="")

message( paste("Merging germline variants for ", study.name, sep=""))

germline.variant.info.ls <- mclapply(
	X = tumor.sample.ids.ls,
	FUN = function( curr.sample.id ){
##		curr.sample.id <- tumor.sample.ids.ls[[5]]
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		message(paste("Extracting filtered germline variants for ", curr.sample.id, sep=""))
		
		load( file = paste(purple.dir, "/", curr.subject.id, "/purple/", curr.sample.id, "/", curr.sample.id, ".purple.germline.vcf.filtered.Rdata", sep=""))
		germline.variant.info.dt
	},
	mc.cores = threads)

germline.variant.info.dt <- rbindlist(germline.variant.info.ls)

save(list=c('germline.variant.info.dt'),
	file = paste(run.dir, "/result_summaries/", gpl_prefix, "/germline.variant.info.Rdata", sep=""))

message( paste("Saved filtered PURPLE processed genomewide germline variant data for ", study.name, sep=""))
