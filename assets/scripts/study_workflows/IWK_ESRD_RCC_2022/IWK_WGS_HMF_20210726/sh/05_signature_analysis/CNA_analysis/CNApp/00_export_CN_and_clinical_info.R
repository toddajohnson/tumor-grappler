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

gpl_prefix <- "GRIDSS-2.12.0"

suppressPackageStartupMessages({
			library(data.table)
			library(parallel)})

extra.samples.to.exclude <- c()
tumor.depth.cutoff <- 10;

source( file.path(run.dir, "config_files/common_config.R") )

#ref.dir <- "/home/tjohnson/reference/HMF/38/dbs/ensembl_data_cache"

#HOME.dir <- "~/HGC_mounts/HGC/"
HOME.dir <- "/home/tjohnson"

fig.dir <- file.path( run.dir, "result", "maftools", "figures")

dir.create(file.path( run.dir, "result", "maftools"))
if ( file.exists(fig.dir, recursive=FALSE)==FALSE ){
	dir.create(fig.dir)
}

isofox.dir <- file.path(rnaseq.run.dir, 'result', 'isofox')

#load( file = file.path( run.dir, "result", "maftools", "clinical_data.Rdata" ))
load( file = file.path(run.dir, 'result/sample_summaries/clinical_data_with_colors.Rdata'))
load(file = paste(isofox.dir, "/isofox.gene_expression.edgeR.normalized.Rdata", sep=""))

min.segment.length <- 0
#list=c('cn.segments', 'cn.segments.resegmented.dt'),
load( file = file.path( run.dir, 'result/CN_segments', paste(study.name, '_merged_CN_gt_', min.segment.length, 'bp.Rdata', sep="") ) )

cn.segments.resegmented.dt <- merge(
	x = clinical.data[,list(sample=tumor.sample.id, gender, purity, Histology, years.dialysis)],
	y = cn.segments.resegmented.dt,
	by = c('sample'))

cn.segments.resegmented.dt[,annotation:=ifelse(Histology=='ccRCC', 'ccRCC', 'non-ccRCC')]

cn.segments.resegmented.dt.in.RNA <- cn.segments.resegmented.dt[sample%in%names(cls8),]


#CN.segments <-  cn.segments.resegmented.dt[,list(ID=sample,
#	chr=chromosome, loc.start=startpos, loc.end=endpos,
#	seg.mean=seg.mean.ploidy, purity, BAF, annotation)]


annotations.dt <- clinical.data[,list(ID=tumor.sample.id, Histology, years.dialysis, DialysisPeriod=ifelse(years.dialysis>14.1, '02_Long', ifelse(years.dialysis<14.1, '01_Short', 'UNK')))]

annotations.dt <- clinical.data[,list(ID=tumor.sample.id, Histology, years.dialysis,
	DialysisPeriod=ifelse(years.dialysis>14.1, '02_Long', ifelse(years.dialysis<14.1, '01_Short', 'UNK')),
	DialysisPeriod.2=cut(years.dialysis, breaks=c(0,10,30), labels=c('<10y', '>=10y')),
	DialysisPeriod.3=cut(years.dialysis, breaks=3, labels=c('01_Short', '02_Mid', '03_Long')))]

annotations.dt[,Histology.adj:=ifelse(Histology%in%c('ACD-RCC', 'pRCC'), 'ACD-RCC/pRCC', as.character(Histology))]
annotations.dt[,DialysisPeriod_Histology:=paste(DialysisPeriod, ':', Histology, sep='')]
annotations.dt[,DialysisPeriod_Histology.adj:=paste(DialysisPeriod, ':', Histology.adj, sep='')]

annotations.dt[,Histology_DialysisPeriod:=paste(Histology, ':', DialysisPeriod, sep='')]
annotations.dt[,Histology.adj_DialysisPeriod.:=paste(Histology.adj, ':', DialysisPeriod, sep='')]

annotations.in.RNA.dt <- annotations.dt[ID%in%names(cls8),]
annotations.in.RNA.dt[,cls.3:=cls[ID]]
annotations.in.RNA.dt[,cls.8:=cls8[ID]]

fwrite(
	x = cn.segments.resegmented.dt[,list(ID=sample,
					chr=chromosome, loc.start=startpos, loc.end=endpos,
					seg.mean=seg.mean.ploidy, purity, BAF, annotation)],
	file = file.path( run.dir, "result", "CNApp/CN_segments", "merged_CN_segs_for_CNApp.txt" ), sep="\t")

fwrite(
		x = cn.segments.resegmented.dt[,list(ID=sample,
						chr=chromosome, loc.start=startpos, loc.end=endpos,
						seg.mean=seg.mean.diploid, purity, BAF, annotation)],
		file = file.path( run.dir, "result", "CNApp/CN_segments", "merged_CN_segs_for_CNApp.diploid.txt" ), sep="\t")

fwrite(
		x = cn.segments.resegmented.dt[,list(ID=sample,
						chr=chromosome, loc.start=startpos, loc.end=endpos,
						seg.mean=seg.mean.adj_diploid, BAF, annotation)],
		file = file.path( run.dir, "result", "CNApp/CN_segments", "merged_CN_segs_for_CNApp.adj_diploid.no_purity.txt" ), sep="\t")

fwrite(
		x = cn.segments.resegmented.dt[,list(ID=sample,
						chr=chromosome, loc.start=startpos, loc.end=endpos,
						seg.mean=seg.mean.ploidy, BAF, annotation)],
	file = file.path( run.dir, "result", "CNApp/CN_segments", "merged_CN_segs_for_CNApp.no_purity.txt" ), sep="\t")

fwrite(
		x = cn.segments.resegmented.dt[,list(ID=sample,
						chr=chromosome, loc.start=startpos, loc.end=endpos,
						seg.mean=seg.mean.diploid, BAF, annotation)],
		file = file.path( run.dir, "result", "CNApp/CN_segments", "merged_CN_segs_for_CNApp.diploid.no_purity.txt" ), sep="\t")

fwrite(
		x = cn.segments.resegmented.dt[annotation=='ccRCC',list(ID=sample,
						chr=chromosome, loc.start=startpos, loc.end=endpos,
						seg.mean=seg.mean.diploid, BAF, annotation)],
		file = file.path( run.dir, "result", "CNApp/CN_segments", "merged_CN_segs_for_CNApp.diploid.no_purity.ccRCC.txt" ), sep="\t")

fwrite(
		x = cn.segments.resegmented.dt[annotation=='non-ccRCC',list(ID=sample,
						chr=chromosome, loc.start=startpos, loc.end=endpos,
						seg.mean=seg.mean.diploid, BAF, annotation)],
		file = file.path( run.dir, "result", "CNApp/CN_segments", "merged_CN_segs_for_CNApp.diploid.no_purity.non_ccRCC.txt" ), sep="\t")

fwrite(
		x = cn.segments.resegmented.dt[Histology=='ACD-RCC',list(ID=sample,
						chr=chromosome, loc.start=startpos, loc.end=endpos,
						seg.mean=seg.mean.diploid, BAF, annotation)],
		file = file.path( run.dir, "result", "CNApp/CN_segments", "merged_CN_segs_for_CNApp.diploid.no_purity.ACD-RCC.txt" ), sep="\t")

fwrite(
		x = cn.segments.resegmented.dt[Histology=='pRCC',list(ID=sample,
						chr=chromosome, loc.start=startpos, loc.end=endpos,
						seg.mean=seg.mean.diploid, BAF, annotation)],
		file = file.path( run.dir, "result", "CNApp/CN_segments", "merged_CN_segs_for_CNApp.diploid.no_purity.pRCC.txt" ), sep="\t")

fwrite(
		x = cn.segments.resegmented.dt[Histology=='ccpRCC',list(ID=sample,
						chr=chromosome, loc.start=startpos, loc.end=endpos,
						seg.mean=seg.mean.diploid, BAF, annotation)],
		file = file.path( run.dir, "result", "CNApp/CN_segments", "merged_CN_segs_for_CNApp.diploid.no_purity.ccpRCC.txt" ), sep="\t")

fwrite(
		x = cn.segments.resegmented.dt[Histology=='chRCC',list(ID=sample,
						chr=chromosome, loc.start=startpos, loc.end=endpos,
						seg.mean=seg.mean.diploid, BAF, annotation)],
		file = file.path( run.dir, "result", "CNApp/CN_segments", "merged_CN_segs_for_CNApp.diploid.no_purity.chRCC.txt" ), sep="\t")


fwrite(
		x = cn.segments.resegmented.dt[!Histology%in%c('ccRCC', 'ACD-RCC', 'pRCC'),list(ID=sample,
						chr=chromosome, loc.start=startpos, loc.end=endpos,
						seg.mean=seg.mean.diploid, BAF, annotation)],
		file = file.path( run.dir, "result", "CNApp/CN_segments", "merged_CN_segs_for_CNApp.diploid.no_purity.other.txt" ), sep="\t")

fwrite(
	x = annotations.dt,
	file = file.path( run.dir, "result", "CNApp/CN_segments", "annotation.txt" ), sep="\t")


fwrite(
		x = cn.segments.resegmented.dt.in.RNA[,list(ID=sample,
						chr=chromosome, loc.start=startpos, loc.end=endpos,
						seg.mean=seg.mean.adj_diploid, BAF, annotation)],
		file = file.path( run.dir, "result", "CNApp/CN_segments", "merged_CN_segs_for_CNApp.adj_diploid.no_purity.in.RNA.txt" ), sep="\t")


fwrite(
		x = annotations.in.RNA.dt,
		file = file.path( run.dir, "result", "CNApp/CN_segments", "annotation.in.RNA.txt" ), sep="\t")


save(list=c('cn.segments.resegmented.dt', 'annotations.dt', 'cn.segments.resegmented.dt.in.RNA', 'annotations.in.RNA.dt'),
		file = file.path( run.dir, "result", "CNApp/CN_segments", "merged_CN_segs_for_CNApp.Rdata" ))
