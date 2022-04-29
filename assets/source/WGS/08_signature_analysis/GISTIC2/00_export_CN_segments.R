#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

min.segment.length <- as.integer(args[1])

options(width=350)
options (repos = "http://cran.ism.ac.jp")

working_dir <- getwd()

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "~/workspace/runs", study.dir )

needed.packages <- c("data.table", "parallel")
install.packages(setdiff(needed.packages, rownames(installed.packages()))) 

suppressPackageStartupMessages({
	library(data.table)
	library(parallel)})

# these should be set in common_config but are here for backup
remove.copy.number.noise.samples <- FALSE
tumor.depth.cutoff <- 10

source( file.path(run.dir, "config_files/common_config.R") )

gistic.segment.dir <- file.path( run.dir, "result", "GISTIC2/CN_segments")

if ( file.exists(gistic.segment.dir , recursive=FALSE)==FALSE ){
	system(paste("mkdir -p ", gistic.segment.dir, sep=""))
}

#list=c('cn.segments', 'cn.segments.resegmented.dt'),
if(exists('version.date.for.CN.segments.file')==TRUE){
	load( file = file.path( run.dir, 'result/CN_segments', paste(study.name, '_merged_CN_gt_', min.segment.length, 'bp.', version.date.for.CN.segments.file, '.Rdata', sep="") ) )
}else{
	load( file = file.path( run.dir, 'result/CN_segments', paste(study.name, '_merged_CN_gt_', min.segment.length, 'bp.Rdata', sep="") ) )
}

CN.segments <-  cn.segments.resegmented.dt[,list(sample, chromosome, start=startpos, end=endpos, NumMarkers, Seg.CN=seg.mean.adj_diploid)]

base.GISTIC2.input.cols <-c('sample', 'chromosome', 'start', 'end', 'NumMarkers')

fwrite(
	x = CN.segments[,c(base.GISTIC2.input.cols, 'Seg.CN'), with=FALSE],
	file = file.path( run.dir, "result/GISTIC2/CN_segments", paste("merged_CN_segs.AllSamples.txt", sep="") ), sep="\t")

save(list=c('CN.segments'),
	file = file.path( run.dir, "result/GISTIC2/CN_segments/merged_CN_segs.GISTIC2.Rdata" ))
