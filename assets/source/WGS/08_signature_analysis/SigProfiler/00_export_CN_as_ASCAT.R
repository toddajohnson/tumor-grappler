#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

min.segment.length <- as.integer(args[1])

options(width=350)
working_dir <- getwd()

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "~/workspace/runs", study.dir )

suppressPackageStartupMessages({
			library(data.table)
			library(parallel)})

source( file.path(run.dir, "config_files/common_config.R") )

#ref.dir <- "/home/tjohnson/reference/HMF/38/dbs/ensembl_data_cache"

#HOME.dir <- "~/HGC_mounts/HGC/"
HOME.dir <- "/home/tjohnson"

#list=c('cn.segments', 'cn.segments.resegmented.dt'),
load( file = file.path( run.dir, 'result/CN_segments', paste(study.name, '_merged_CN_gt_', min.segment.length, 'bp.Rdata', sep="") ) )

output.dir <- paste(run.dir, '/result/SigProfiler/', study.name, '_resegmented_CN_gt_', min.segment.length, 'bp_signature_analysis', sep="")

if ( length(list.dirs(output.dir, recursive=FALSE))==0 ){
	system(paste("mkdir -p ", output.dir, sep=""))
}

print( paste("Writing resegmented CN as ASCAT file in ", output.dir, sep="") )

fwrite(
	x = cn.segments.resegmented.dt[,list(sample, chr=chromosome, startpos, endpos, nMajor, nMinor)],
	file = paste(output.dir, '/', study.name, '_resegmented_merged_CN_gt_', min.segment.length, 'bp.txt', sep=""),
	col.names=TRUE, sep='\t')

print(cn.segments.resegmented.dt[,list(sample, chromosome, startpos, endpos, nMajor, nMinor)])
