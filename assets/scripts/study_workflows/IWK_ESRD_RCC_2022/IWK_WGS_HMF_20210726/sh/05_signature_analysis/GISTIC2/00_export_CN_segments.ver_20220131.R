#!/usr/bin/env Rscript

# module use /usr/local/package/modulefiles/
# module load R/4.0.2
# . ~/.R-4.1.0_setup.sh

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

remove.copy.number.noise.samples <- FALSE
tumor.depth.cutoff <- 10

source( file.path(run.dir, "config_files/common_config.R") )

#ref.dir <- "/home/tjohnson/reference/HMF/38/dbs/ensembl_data_cache"

#HOME.dir <- "~/HGC_mounts/HGC/"
HOME.dir <- "/home/tjohnson"

gistic.segment.dir <- file.path( run.dir, "result", "GISTIC2/CN_segments")

if ( !file.exists(gistic.segment.dir , recursive=FALSE)==FALSE ){
	system(paste("mkdir -p ", gistic.segment.dir, sep=""))
}

#load( file = file.path( run.dir, "result", "maftools", "clinical_data.Rdata" ))
load( file = file.path(run.dir, 'result/sample_summaries/clinical_data_with_colors.Rdata'))

#histology.short <- c('clear_cell', 'chromophobe', 'papillary', 'clear_papillary', 'ACD_RCC', 'unclassified')
#names(histology.short) <- c('Clear cell', 'Chromophobe', 'Papillary', 'Papillary (Clear papillary)', 'ACD-RCC', 'Unclassified')

min.segment.length <- 0
#list=c('cn.segments', 'cn.segments.resegmented.dt'),
load( file = file.path( run.dir, 'result/CN_segments', paste(study.name, '_merged_CN_gt_', min.segment.length, 'bp.ver_20220131.Rdata', sep="") ) )

cn.segments.resegmented.dt <- merge(
	x = clinical.data[,list(sample=tumor.sample.id, Histology)],
	y = cn.segments.resegmented.dt,
	by = c('sample'))

cn.segments.resegmented.dt[,annotation:=ifelse(Histology=='ccRCC', 'ccRCC', 'non-ccRCC')]
CN.segments <-  cn.segments.resegmented.dt[,list(sample, annotation,
	chromosome, start=startpos, end=endpos, NumMarkers, Seg.CN=seg.mean.adj_diploid)]
#	chromosome, start=startpos, end=endpos, NumMarkers, Seg.CN=seg.mean.diploid)]

annotation.ls <- unique(CN.segments$annotation)
base.GISTIC2.input.cols <-c('sample', 'chromosome', 'start', 'end', 'NumMarkers')

fwrite(
	x = CN.segments[,c(base.GISTIC2.input.cols, 'Seg.CN'), with=FALSE],
	file = file.path( run.dir, "result/GISTIC2/CN_segments", paste("merged_CN_segs.ver_20220131.AllSamples.txt", sep="") ), sep="\t")

catch.output.ls <- lapply(
	X = annotation.ls,
	FUN = function(curr.annotation.type){
		fwrite(
			x = CN.segments[annotation==curr.annotation.type, c(base.GISTIC2.input.cols, 'Seg.CN'), with=FALSE],
			file = file.path( run.dir, "result/GISTIC2/CN_segments", paste("merged_CN_segs.ver_20220131.", curr.annotation.type, ".txt", sep="") ), sep="\t")
	})

save(list=c('annotation.ls', 'CN.segments'),
		file = file.path( run.dir, "result/GISTIC2/CN_segments/merged_CN_segs.GISTIC2.ver_20220131.Rdata" ))


cn.segments.resegmented.dt[,annotation:=ifelse(Histology%in%c('ccRCC', 'ACD-RCC', 'pRCC'), as.character(Histology), 'other-RCC')]
#CN.segments <-  cn.segments.resegmented.dt[,list(sample, annotation, chromosome, start=startpos, end=endpos, NumMarkers, Seg.CN=seg.mean.diploid)]
CN.segments <-  cn.segments.resegmented.dt[,list(sample, annotation, chromosome, start=startpos, end=endpos, NumMarkers, Seg.CN=seg.mean.adj_diploid)]

annotation.ls <- unique(CN.segments$annotation)
base.GISTIC2.input.cols <-c('sample', 'chromosome', 'start', 'end', 'NumMarkers')

catch.output.ls <- lapply(
	X = annotation.ls,
	FUN = function(curr.annotation.type){
		fwrite(
			x = CN.segments[annotation==curr.annotation.type, c(base.GISTIC2.input.cols, 'Seg.CN'), with=FALSE],
			file = file.path( run.dir, "result/GISTIC2/CN_segments", paste("merged_CN_segs.ver_20220131.", curr.annotation.type, ".txt", sep="") ), sep="\t")
	})

save(list=c('annotation.ls', 'CN.segments'),
		file = file.path( run.dir, "result/GISTIC2/CN_segments/merged_CN_segs.GISTIC2.more_groups.ver_20220131.Rdata" ))
