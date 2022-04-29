#!/usr/bin/env Rscript

# module use /usr/local/package/modulefiles/
# module load R/4.0.2
# . ~/.R-4.1.0_setup.sh


#AD.threshold <- 0.1
#GISTIC2.date <- "20220202_median"

options(width=220)
working_dir <- getwd()
gpl_prefix <- basename(working_dir)

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "~/workspace/runs", study.dir )

gpl_prefix <- "GRIDSS-2.12.0"

suppressPackageStartupMessages({
			library(data.table)
			library(parallel)})

library(stringr)
#library(R.matlab)

HGNC.dir <- "~/reference/HGNC"
load(file = file.path(HGNC.dir, 'hgnc_complete_set_20220131.Rdata'))

source( file.path(run.dir, "config_files/common_config.R") )

#ref.dir <- "/home/tjohnson/reference/HMF/38/dbs/ensembl_data_cache"

#HOME.dir <- "~/HGC_mounts/HGC/"
#HOME.dir <- "/home/tjohnson"

fig.dir <- file.path( run.dir, "result", "maftools", "figures")

#ref.dir <- "/home/tjohnson/reference/HMF/38/dbs/linx"
#load(file = paste(ref.dir, "/ensembl.gene.info.Rdata", sep=""))

##hg38.genes <- readMat( file.path('~/tools/GISTIC2/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat') )

print("Loading cytoband info")
cytoband.dir <- '~/tools/snpEff/data/GRCh38.99'
cytoband.dt <- fread(file.path(cytoband.dir, 'cytoBand.txt.gz'), col.names=c('chrom', 'start', 'end', 'cytoband', 'stain'))

chrom.map.ls <- c(1:22, 'X', 'Y', 'M')
names(chrom.map.ls) <- paste('chr', c(1:22, 'X', 'Y', 'M'), sep='')

cytoband.dt[,chr:=chrom.map.ls[chrom]]
cytoband.dt <- cytoband.dt[!is.na(chr),]
cytoband.dt[,chromband:=paste(chr, cytoband, sep='')]
cytoband.dt[,chromarm:=chrom]
cytoband.dt[grep('p', chromband),chromarm:=paste(chromarm, 'p', sep='')]
cytoband.dt[grep('q', chromband),chromarm:=paste(chromarm, 'q', sep='')]


#load( file = file.path( run.dir, "result", "maftools", "clinical_data.Rdata" ))
load( file = file.path(run.dir, 'result/sample_summaries/clinical_data_with_colors.Rdata'))

clinical.data <- merge(
		x = clinical.data,
		y = clinical.data[,list(subtype.ct=.N), by=list(Histology)],
		by = 'Histology')

chrom.map.ls <- c(as.character(1:22), 'X', 'Y')
names(chrom.map.ls) <- as.character(1:24)

annotation.ls <- unique(as.character(clinical.data$Histology))
annotation.ls <- c('AllSamples', 'ccRCC', 'non-ccRCC', 'ACD-RCC', 'pRCC', 'ccpRCC', 'chRCC')

cna_freqs.dt.ls <- lapply(
	X = annotation.ls,
	FUN = function(curr.annotation.type){
			curr.dir <- file.path( run.dir, "result/CNApp/figures/AllSamples_adj_diploid")
			curr.file <- system( paste('ls ', curr.dir, '/CNAs_freq_', curr.annotation.type, '.tsv', sep=""), intern=TRUE)
			curr.dt <- fread( file = curr.file, sep="\t")
			curr.cols <- copy(colnames(curr.dt))
			curr.dt[,annotation:=curr.annotation.type]
			curr.dt[,c('annotation', curr.cols), with=FALSE]
		})
cna_freqs.dt <- rbindlist(cna_freqs.dt.ls)
cna_freqs.dt[,chrom:=strsplit(Region, split=':', fixed=TRUE)[[1]][[1]], by=1:nrow(cna_freqs.dt)]
cna_freqs.dt[,bp_range:=strsplit(Region, split=':', fixed=TRUE)[[1]][[2]], by=1:nrow(cna_freqs.dt)]
cna_freqs.dt[,start:=as.integer(strsplit(bp_range, split='-', fixed=TRUE)[[1]][[1]]), by=1:nrow(cna_freqs.dt)]
cna_freqs.dt[,end:=as.integer(strsplit(bp_range, split='-', fixed=TRUE)[[1]][[2]]), by=1:nrow(cna_freqs.dt)]


## cna_freqs.range.info <- unique(cna_freqs.dt[,list(chrom, start, end)])
## 
## annotate_range_with_cytoband <- function( curr.chrom, curr.start, curr.end ){
##     curr.cytoband.dt <- cytoband.dt[chrom==curr.chrom & start<curr.end & end>curr.start,]
##     data.table(
##         cytobands = paste(curr.cytoband.dt$cytoband, collapse=','),
##         chromarm = unique(curr.cytoband.dt$chromarm))
## }
## cna_freqs.range.info.with.cytobands <- cna_freqs.range.info[,annotate_range_with_cytoband(chrom, start, end), by=list(chrom, start, end)]

##cna_freqs.dt[,chromband:=strsplit(regions, split='_', fixed=TRUE)[[1]][[2]], by=1:nrow(cna_freqs.dt)]


save(list=c('cna_freqs.dt'),
	file = file.path( run.dir, "result/CNApp/merged_recurrent_cna_freqs.Rdata" ))
