#!/usr/bin/env Rscript

# module use /usr/local/package/modulefiles/
# module load R/4.0.2
# . ~/.R-4.1.0_setup.sh

options(width=350)
working_dir <- getwd()
gpl_prefix <- basename(working_dir)

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "~/workspace/runs", study.dir )
#wgs.run.dir <- '/Users/toddjohnson/OneDrive/Documents/RIKEN/Projects/IWK_ESRD_RCC/IWK_WGS_HMF_20210726'
#run.dir <- '/Users/toddjohnson/OneDrive/Documents/RIKEN/Projects/IWK_ESRD_RCC/IWK_RNAseq_20210805'
#wgs.run.dir <- '/Users/tajohnson/OneDrive/Documents/RIKEN/Projects/IWK_ESRD_RCC/IWK_WGS_HMF_20210726'
#run.dir <- '/Users/tajohnson/OneDrive/Documents/RIKEN/Projects/IWK_ESRD_RCC/IWK_RNAseq_20210805'

library(data.table)
library(parallel)
library(tidyverse)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(edgeR)
library(limma)


source( file.path(run.dir, "config_files/common_config.R") )

#ref.dir <- "/home/tjohnson/reference/HMF/38/dbs/ensembl_data_cache"

#HOME.dir <- "~/HGC_mounts/HGC/"
HOME.dir <- "/home/tjohnson"

isofox.dir <- file.path(run.dir, 'result', 'isofox')

load(file = paste(isofox.dir, "/isofox.gene_expression.edgeR.normalized.Rdata", sep=""))
load(file = paste(wgs.run.dir, "/result/sample_summaries/clinical_data_with_colors.Rdata", sep=""))
load(file = file.path(wgs.run.dir, 'result/variant_summaries/small_variant_summaries.Rdata'))

chr3.gene.drivers.dt <- merged.drivers.dt[gene%in%c('VHL', 'PBRM1', 'SETD2', 'BAP1'),
	list(subject.id, tumor.sample.id, Histology, purity, ploidy, years.dialysis, chromosome, chromosomeBand, gene)]
chr3.gene.drivers.dt <- chr3.gene.drivers.dt[order(tumor.sample.id, gene),]
chr3.gene.drivers.summary <- chr3.gene.drivers.dt[,list(chr.genes=paste(gene, collapse=';')), by=list(subject.id, tumor.sample.id, Histology, purity, ploidy, years.dialysis)][order(years.dialysis),]
## 	subject.id tumor.sample.id Histology purity ploidy years.dialysis       chr.genes
## 1:     IWK048        IWK048_T      pRCC   0.53   2.26    0.008333333             VHL
## 2:     IWK022        IWK022_T     ccRCC   0.50   1.90    0.294444444       PBRM1;VHL
## 3:     IWK010        IWK010_T     ccRCC   0.38   2.04    0.416666667 PBRM1;SETD2;VHL
## 4:     IWK033        IWK033_T     ccRCC   0.68   1.82    1.094444444       PBRM1;VHL
## 5:     IWK002        IWK002_T     ccRCC   0.58   1.90    2.216666667        BAP1;VHL
## 6:     IWK009        IWK009_T     ccRCC   0.63   1.88    3.083333333             VHL
## 7:     IWK025        IWK025_T     ccRCC   0.34   1.92    4.666666667           PBRM1
## 8:     IWK021        IWK021_T     ccRCC   0.47   1.90    5.502777778             VHL
## 9:     IWK005        IWK005_T     ccRCC   0.69   2.04    6.619444444           PBRM1
## 10:     IWK006        IWK006_T     ccRCC   0.73   1.28    6.780555556       PBRM1;VHL
## 11:     IWK020        IWK020_T     ccRCC   0.57   3.20    7.672222222             VHL
## 12:     IWK012        IWK012_T     ccRCC   0.74   2.00    8.166666667           PBRM1
## 13:     IWK024        IWK024_T     ccRCC   0.53   1.90    9.025000000       PBRM1;VHL
## 14:     IWK029        IWK029_T     ccRCC   0.46   1.94   22.030555556             VHL
## 15:     IWK015        IWK015_T     ccRCC   0.72   1.62   25.647222222             VHL

setkey(clinical.data, tumor.sample.id)

sample.ids <- colnames(mat.all.samples.natural.cpm)

sample.info.with.clinical.info.dt <- clinical.data[sample.ids,]
sample.info.with.clinical.info.dt[names(cls),cls.3:=cls]
sample.info.with.clinical.info.dt[names(cls8),cls.8:=cls8]

sample.info.with.clinical.info.dt[chr3.gene.drivers.summary$tumor.sample.id,ccRCC.class.status:=chr3.gene.drivers.summary$chr.genes]
sample.info.with.clinical.info.dt[,ccRCC.class.status:=ifelse(ccRCC.class.status=='BAP1;VHL','VHL', ifelse(ccRCC.class.status=='PBRM1;SETD2;VHL', 'PBRM1;VHL', ccRCC.class.status))]
sample.info.with.clinical.info.dt[,ccRCC.class.status:=ifelse(is.na(ccRCC.class.status), cls.3, ccRCC.class.status)]
sample.info.with.clinical.info.dt[,DialysisPeriod.status:=cut(years.dialysis, breaks=2, labels=c('Short', 'Long'), ordered=TRUE)]
sample.info.with.clinical.info.dt[,DialysisPeriod.status:=ifelse(cls.3%in%c('B', 'C'), as.character(DialysisPeriod.status), 'none')]

#ensembl.gene.file <- '/Volumes/ToddsDocs/references/HMF/38/Ensembl-Data-Cache/38/ensembl_gene_data.csv'
#ensembl.gene.file <- '/Users/tajohnson/Documents/references/HMF/38/Ensembl-Data-Cache/38/ensembl_gene_data.csv'
ensembl.gene.file <- '/home/tjohnson/reference/HMF/38/dbs/linx/ensembl_gene_data.csv'
ensembl.genes <- fread(ensembl.gene.file)
ensembl.genes[,chromband:=paste('chr', Chromosome, KaryotypeBand, sep="")]
setkey(ensembl.genes, 'GeneName')

mat.all.samples.natural.cpm.dt <- as.data.table(mat.all.samples.natural.cpm, keep.rownames=TRUE)
setnames(mat.all.samples.natural.cpm.dt, old='rn', new='NAME')
mat.all.samples.natural.cpm.dt[,Description:='NA']
setkey(mat.all.samples.natural.cpm.dt, NAME)

ensembl.genes.in.file <- ensembl.genes[mat.all.samples.natural.cpm.dt$NAME,]
#> nrow(ensembl.genes.in.file)
# 20682
#[1] 20229

ensembl.genes.in.file.no.loc <- ensembl.genes.in.file[is.na(Chromosome)]
ensembl.genes.in.file.loc <- ensembl.genes.in.file[!is.na(Chromosome)]

mat.all.samples.natural.cpm.dt.loc <- mat.all.samples.natural.cpm.dt[ensembl.genes.in.file.loc$GeneName,]

gmt.dt <- ensembl.genes.in.file.loc[,list(Geneset=paste(GeneName, collapse='\t'), description='Cytogenetic sub-band'), by=list(chromband)]
#gmt.file <- file.path(isofox.dir, "result/GSEA/isofox_gene_expression.gmt")
gmt.file <- file.path(run.dir, "result/GSEA/isofox_gene_expression.gmt")
fwrite(x=gmt.dt[,list(chromband, description, Geneset)], file=gmt.file, sep='\t', col.names=FALSE, append=FALSE, quote=FALSE)



gene.ct <- nrow(ensembl.genes.in.file.loc)
sample.ct <- nrow(sample.info.with.clinical.info.dt)

gct.dt <- mat.all.samples.natural.cpm.dt.loc[,c('NAME', 'Description', sample.info.with.clinical.info.dt$tumor.sample.id), with=FALSE]
#length(gct.dt$NAME)
#[1] 20682

#gct.file <- file.path(isofox.dir, "result/GSEA/isofox_gene_expression.gct")
gct.file <- file.path(run.dir, "result/GSEA/isofox_gene_expression.gct")

writeLines(c('#1.2', paste(c(gene.ct, sample.ct), collapse='\t')), gct.file)
fwrite(gct.dt, gct.file, sep='\t', col.names=TRUE, append=TRUE)

#ccRCC.class.cls <- file.path(isofox.dir, "result/GSEA/ccRCC_class.cls")
ccRCC.class.cls <- file.path(run.dir, "result/GSEA/ccRCC_class.cls")

writeLines( c(
	paste(c(sample.ct, as.character(length(unique(sample.info.with.clinical.info.dt$ccRCC.class.status))), '1'), collapse=' '),
	paste( c("#", unique(sample.info.with.clinical.info.dt$ccRCC.class.status)), collapse=' '),
	paste( c(sample.info.with.clinical.info.dt$ccRCC.class.status), collapse=' ')), ccRCC.class.cls)


chRCC.non.chRCC <- ifelse(sample.info.with.clinical.info.dt$Histology=='chRCC', 'chRCC', 'non-chRCC')
#chRCC.non.chRCC.cls <- file.path(isofox.dir, "result/GSEA/chRCC_nonchRCC.cls")
chRCC.non.chRCC.cls <- file.path(run.dir, "result/GSEA/chRCC_nonchRCC.cls")

writeLines( c(
				paste(c(sample.ct, '2', '1'), collapse=' '),
				paste( c("#", unique(as.character(chRCC.non.chRCC))), collapse=' '),
				paste( c(chRCC.non.chRCC), collapse=' ')), chRCC.non.chRCC.cls)


sample.info.with.clinical.info.dt.years.dialysis <- sample.info.with.clinical.info.dt[DialysisPeriod.status!='none',]
years.dialysis.2cuts <- as.character(sample.info.with.clinical.info.dt.years.dialysis$DialysisPeriod.status)
#years.dialysis.2cuts.cls <- file.path(isofox.dir, "result/GSEA/years.dialysis_2.cls")
years.dialysis.2cuts.cls <- file.path(run.dir, "result/GSEA/years.dialysis_2.cls")

years.dialysis.sample.ct <- nrow(sample.info.with.clinical.info.dt.years.dialysis)

writeLines( c(
				paste(c(years.dialysis.sample.ct, '2', '1'), collapse=' '),
				paste( c("#", unique(as.character(years.dialysis.2cuts))), collapse=' '),
				paste( c(years.dialysis.2cuts), collapse=' ')), years.dialysis.2cuts.cls)

years.dialysis.gct.dt <- mat.all.samples.natural.cpm.dt.loc[,c('NAME', 'Description', sample.info.with.clinical.info.dt.years.dialysis$tumor.sample.id), with=FALSE]
#years.dialysis.gct.file <- file.path(isofox.dir, "result/GSEA/isofox_gene_expression.years.dialysis.gct")
years.dialysis.gct.file <- file.path(run.dir, "result/GSEA/isofox_gene_expression.years.dialysis.gct")

writeLines(c('#1.2', paste(c(gene.ct, years.dialysis.sample.ct), collapse='\t')), years.dialysis.gct.file)
fwrite(years.dialysis.gct.dt, years.dialysis.gct.file, sep='\t', col.names=TRUE, append=TRUE)



sample.info.with.clinical.info.dt.cls3.B.C <- sample.info.with.clinical.info.dt[cls.3%in%c('B', 'C'),]
cls3.B.C.cls <- file.path(run.dir, "result/GSEA/cls3.B.C.cls")

cls3.B.C.sample.ct <- nrow(sample.info.with.clinical.info.dt.cls3.B.C)

writeLines( c(
				paste(c(cls3.B.C.sample.ct, '2', '1'), collapse=' '),
				paste( c("#", unique(sample.info.with.clinical.info.dt.cls3.B.C$cls.3)), collapse=' '),
				paste( c(sample.info.with.clinical.info.dt.cls3.B.C$cls.3), collapse=' ')), cls3.B.C.cls)

cls3.B.C.gct.dt <- mat.all.samples.natural.cpm.dt.loc[,c('NAME', 'Description', sample.info.with.clinical.info.dt.cls3.B.C$tumor.sample.id), with=FALSE]
cls3.B.C.gct.file <- file.path(run.dir, "result/GSEA/isofox_gene_expression.cls3.B.C.gct")

writeLines(c('#1.2', paste(c(gene.ct, cls3.B.C.sample.ct), collapse='\t')), cls3.B.C.gct.file)
fwrite(cls3.B.C.gct.dt, cls3.B.C.gct.file, sep='\t', col.names=TRUE, append=TRUE)



sample.info.with.clinical.info.dt.cls8.B.H <- sample.info.with.clinical.info.dt[cls.8%in%c('B', 'H'),]
cls8.B.H.cls <- file.path(run.dir, "result/GSEA/cls8.B.H.cls")

cls8.B.H.sample.ct <- nrow(sample.info.with.clinical.info.dt.cls8.B.H)

writeLines( c(
				paste(c(cls8.B.H.sample.ct, '2', '1'), collapse=' '),
				paste( c("#", unique(sample.info.with.clinical.info.dt.cls8.B.H$cls.8)), collapse=' '),
				paste( c(sample.info.with.clinical.info.dt.cls8.B.H$cls.8), collapse=' ')), cls8.B.H.cls)

cls8.B.H.gct.dt <- mat.all.samples.natural.cpm.dt.loc[,c('NAME', 'Description', sample.info.with.clinical.info.dt.cls8.B.H$tumor.sample.id), with=FALSE]
cls8.B.H.gct.file <- file.path(run.dir, "result/GSEA/isofox_gene_expression.cls8.B.H.gct")

writeLines(c('#1.2', paste(c(gene.ct, cls8.B.H.sample.ct), collapse='\t')), cls8.B.H.gct.file)
fwrite(cls8.B.H.gct.dt, cls8.B.H.gct.file, sep='\t', col.names=TRUE, append=TRUE)



sample.info.with.clinical.info.dt.ccRCC <- sample.info.with.clinical.info.dt[Histology=='ccRCC',]

ccRCC.gct.dt <- mat.all.samples.natural.cpm.dt.loc[,c('NAME', 'Description', sample.info.with.clinical.info.dt.ccRCC$sample.id), with=FALSE]
ccRCC.sample.ct <- nrow(sample.info.with.clinical.info.dt.ccRCC)
		
#ccRCC.gct.file <- file.path(isofox.dir, "result/GSEA/isofox_gene_expression.ccRCC.gct")
ccRCC.gct.file <- file.path(run.dir, "result/GSEA/isofox_gene_expression.ccRCC.gct")

writeLines(c('#1.2', paste(c(gene.ct, ccRCC.sample.ct), collapse='\t')), ccRCC.gct.file)
fwrite(ccRCC.gct.dt, ccRCC.gct.file, sep='\t', col.names=TRUE, append=TRUE)

years.dialysis.ccRCC <- cut(sample.info.with.clinical.info.dt.ccRCC$years.dialysis, breaks=2, labels=c('Short', 'Long'), ordered=TRUE)
#years.dialysis.ccRCC.cls <- file.path(isofox.dir, "result/GSEA/years.dialysis.ccRCC.cls")
years.dialysis.ccRCC.cls <- file.path(run.dir, "result/GSEA/years.dialysis.ccRCC.cls")


writeLines( c(
				paste(c(ccRCC.sample.ct, '2', '1'), collapse=' '),
				paste( c("#", unique(as.character(years.dialysis.ccRCC))), collapse=' '),
				paste( c(years.dialysis.ccRCC), collapse=' ')), years.dialysis.ccRCC.cls)




sample.info.with.clinical.info.dt.ccRCC_ACD.RCC <- sample.info.with.clinical.info.dt[Histology%in%c('ccRCC', 'ACD-RCC'),]
ccRCC_ACD.RCC.labels <- as.character(sample.info.with.clinical.info.dt.ccRCC_ACD.RCC$Histology)
ccRCC_ACD.RCC.cls <- file.path(run.dir, "result/GSEA/ccRCC_ACD.RCC.cls")

ccRCC_ACD.RCC.sample.ct <- nrow(sample.info.with.clinical.info.dt.ccRCC_ACD.RCC)

writeLines( c(
				paste(c(ccRCC_ACD.RCC.sample.ct, '2', '1'), collapse=' '),
				paste( c("#", unique(as.character(ccRCC_ACD.RCC.labels))), collapse=' '),
				paste( c(ccRCC_ACD.RCC.labels), collapse=' ')), ccRCC_ACD.RCC.cls)

ccRCC_ACD.RCC.gct.dt <- mat.all.samples.natural.cpm.dt.loc[,c('NAME', 'Description', sample.info.with.clinical.info.dt.ccRCC_ACD.RCC$tumor.sample.id), with=FALSE]
ccRCC_ACD.RCC.gct.file <- file.path(run.dir, "result/GSEA/isofox_gene_expression.ccRCC_ACD.RCC.gct")

writeLines(c('#1.2', paste(c(gene.ct, ccRCC_ACD.RCC.sample.ct), collapse='\t')), ccRCC_ACD.RCC.gct.file)
fwrite(ccRCC_ACD.RCC.gct.dt, ccRCC_ACD.RCC.gct.file, sep='\t', col.names=TRUE, append=TRUE)

#sample.info.with.clinical.info.dt.cls.C <- sample.info.with.clinical.info.dt[cls=='C',]
#cls.C.gct.dt <- mat.all.samples.natural.cpm.dt.loc[,c('NAME', 'Description', sample.info.with.clinical.info.dt.cls.C$tumor.sample.id), with=FALSE]
#cls.C.sample.ct <- nrow(sample.info.with.clinical.info.dt.cls.C)
#
#cls.C.gct.file <- file.path(run.dir, "result/GSEA/isofox_gene_expression.cls.C.gct")
#
#writeLines(c('#1.2', paste(c(gene.ct, cls.C.sample.ct), collapse='\t')), cls.C.gct.file)
#fwrite(cls.C.gct.dt, cls.C.gct.file, sep='\t', col.names=TRUE, append=TRUE)

#years.dialysis.cls.C <- as.character(sample.info.with.clinical.info.dt.cls.C$DialysisPeriod.status)
#years.dialysis.cls.C.cls <- file.path(run.dir, "result/GSEA/years.dialysis.cls.C.cls")
#
#writeLines( c(
#				paste(c(cls.C.sample.ct, '2', '1'), collapse=' '),
#				paste( c("#", unique(as.character(years.dialysis.cls.C))), collapse=' '),
#				paste( c(years.dialysis.cls.C), collapse=' ')), years.dialysis.cls.C.cls)
#
#
#sample.info.with.clinical.info.dt.cls.B <- sample.info.with.clinical.info.dt[cls=='B',]
#cls.B.gct.dt <- mat.all.samples.natural.cpm.dt.loc[,c('NAME', 'Description', sample.info.with.clinical.info.dt.cls.B$tumor.sample.id), with=FALSE]
#cls.B.sample.ct <- nrow(sample.info.with.clinical.info.dt.cls.B)
#
#cls.B.gct.file <- file.path(run.dir, "result/GSEA/isofox_gene_expression.cls.B.gct")
#
#writeLines(c('#1.2', paste(c(gene.ct, cls.B.sample.ct), collapse='\t')), cls.B.gct.file)
#fwrite(cls.B.gct.dt, cls.B.gct.file, sep='\t', col.names=TRUE, append=TRUE)
#
#years.dialysis.cls.B <- as.character(sample.info.with.clinical.info.dt.cls.B$DialysisPeriod.status)
#years.dialysis.cls.B.cls <- file.path(run.dir, "result/GSEA/years.dialysis.cls.B.cls")
#
#writeLines( c(
#				paste(c(cls.B.sample.ct, '2', '1'), collapse=' '),
#				paste( c("#", unique(as.character(years.dialysis.cls.B))), collapse=' '),
#				paste( c(years.dialysis.cls.B), collapse=' ')), years.dialysis.cls.B.cls)
