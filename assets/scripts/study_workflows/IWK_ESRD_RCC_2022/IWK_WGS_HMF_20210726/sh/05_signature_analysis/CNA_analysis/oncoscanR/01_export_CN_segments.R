#!/usr/bin/env Rscript

# module use /usr/local/package/modulefiles/
# module load R/4.0.2
# . ~/.R-4.1.0_setup.sh

options(width=350)
#options(width=210)
working_dir <- getwd()

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "~/workspace/runs", study.dir )

#run.dir <- '/Users/tajohnson/Library/CloudStorage/OneDrive-Personal/Documents/RIKEN/Projects/IWK_ESRD_RCC/IWK_WGS_HMF_20210726'

# install oncoscanR from https://github.com/yannchristinat/oncoscanR-public
suppressPackageStartupMessages({
			library(data.table)
			library(parallel)})
library(oncoscanR)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)

remove.copy.number.noise.samples <- FALSE
tumor.depth.cutoff <- 10

source( file.path(run.dir, "config_files/common_config.R") )

#run.dir <- '/Users/tajohnson/Library/CloudStorage/OneDrive-Personal/Documents/RIKEN/Projects/IWK_ESRD_RCC/IWK_WGS_HMF_20210726'

#ref.dir <- "/home/tjohnson/reference/HMF/38/dbs/ensembl_data_cache"

HOME.dir <- "/home/tjohnson"

oncoscan.segment.dir <- file.path( run.dir, "result", "oncoscanR/CN_segments")

if ( file.exists(oncoscan.segment.dir , recursive=FALSE)==FALSE ){
	system(paste("mkdir -p ", oncoscan.segment.dir, sep=""))
}

load( file = file.path(run.dir, 'result/sample_summaries/clinical_data_with_colors.Rdata'))

chrom.map.ls <- as.integer(c(1:24))
names(chrom.map.ls) <- paste('chr', c(1:22, 'X', 'Y'), sep="")

min.segment.length <- 0
#list=c('cn.segments', 'cn.segments.resegmented.dt'),
load( file = file.path( run.dir, 'result/CN_segments', paste(study.name, '_merged_CN_gt_', min.segment.length, 'bp.Rdata', sep="") ) )

cn.segments.resegmented.dt[,chrom:=chrom.map.ls[chromosome]]

cn.segments.resegmented.dt <- merge(
	x = clinical.data[,list(sample=tumor.sample.id, Sex=ifelse(gender=='MALE', 'M', 'F'), Histology)],
	y = cn.segments.resegmented.dt)


cn.segments.resegmented.dt[,total.chrom.length:=sum(seg.length), by=list(sample, chromosome)]

cn.segments.resegmented.dt[,full.loc:=paste(chromosome, ':', startpos, '-', endpos, sep='')]
cn.segments.resegmented.dt[,seg.length.propn:=seg.length/total.chrom.length]
cn.segments.resegmented.dt[,type:='none']
cn.segments.resegmented.dt[chrom%in%1:22 | (chrom==23 & Sex=='F'),
	type:=ifelse(BAF>0.8 & nTotal==2, 'LOH', ifelse(nTotal<2, 'Loss', ifelse(nTotal>2, 'Gain', 'Neutral')))]
cn.segments.resegmented.dt[chrom==24 | (chrom==23 & Sex=='M'),type:=ifelse(nTotal<1, 'Loss', ifelse(nTotal>1, 'Gain', 'Neutral'))]

clinical.data <- clinical.data[tumor.sample.id%in%cn.segments.resegmented.dt$sample,]

sample.ids.ls <- clinical.data$sample
names(sample.ids.ls) <- sample.ids.ls

catch.results.ls <- lapply(
	X = sample.ids.ls,
	FUN = function(curr.sample.id){
		
		cn.segments <- cn.segments.resegmented.dt[sample==curr.sample.id,
			list(nTotal, type, NumMarkers, full.loc)]
		
		setnames(cn.segments,
			old=c('nTotal', 'NumMarkers', 'type', 'full.loc'),
			new=c('CN State', 'Marker Count', 'Type', 'Full Location'))
		
		
		fwrite(
				cn.segments[which(Type!='Neutral'),],
				file = file.path(oncoscan.segment.dir, paste('CN_segments.', curr.sample.id, '.txt', sep='')),
				sep = '\t',
				col.names = TRUE)
	})

clinical.data[,Sex:=ifelse(gender=='MALE', 'M', 'F')]

chrom.arms.dt <- data.table(
		chrom = c(c(1:24), c(1:24)),
		chromosome = c(c(1:22, 'X', 'Y'), c(1:22, 'X', 'Y')),
		chrom.arm = c(paste(c(1:22, 'X', 'Y'), 'p', sep=""), paste(c(1:22, 'X', 'Y'), 'q', sep="")))

chrom.arms.dt <- chrom.arms.dt[!chrom.arm%in%paste(c('13', '14', '15', '22'), 'p', sep=""),]
chrom.arms.dt <- chrom.arms.dt[order(chrom, chrom.arm),]
chrom.arms.dt[,chrom.arm.f:=factor(chrom.arm, levels=chrom.arms.dt$chrom.arm, ordered=TRUE)]


CN.oncoscan.ls <- mclapply(
	X = sample.ids.ls,
	FUN = function(curr.sample.id){
#		curr.sample.id <- sample.ids.ls[[1]]
		curr.file <- file.path(oncoscan.segment.dir, paste('CN_segments.', curr.sample.id, '.txt', sep=''))
		curr.gender <- clinical.data[tumor.sample.id==curr.sample.id,]$Sex

		oncoscan.out <- workflow_oncoscan.run(chas.fn = curr.file, gender = curr.gender)
	},
	mc.cores = 8) 

CN.oncoscan.dt.ls <- lapply(
	X = sample.ids.ls,
	FUN = function(curr.sample.id){
		print(curr.sample.id)
		oncoscan.out <- CN.oncoscan.ls[[curr.sample.id]]
		curr.AMP <- oncoscan.out$armlevel$AMP
		curr.LOSS <- oncoscan.out$armlevel$LOSS
		curr.LOH <- oncoscan.out$armlevel$LOH
		curr.GAIN <- oncoscan.out$armlevel$GAIN
		
		curr.dt <- data.table(
			chrom.arm = c(
				curr.AMP,
				curr.LOSS,
				curr.LOH,
				curr.GAIN
			),
			Type = c(
				rep('AMP', length(curr.AMP)),
				rep('LOSS', length(curr.LOSS)),
				rep('LOH', length(curr.LOH)),
				rep('GAIN', length(curr.GAIN))
			))
		if(nrow(curr.dt)>0){
			curr.dt[,tumor.sample.id:=curr.sample.id]
			curr.dt[,list(tumor.sample.id, chrom.arm, Type)]
		}else{
			curr.dt <- data.table(
				tumor.sample.id = rep(curr.sample.id, length(chrom.arms.dt$chrom.arm)),
				chrom.arm = chrom.arms.dt$chrom.arm,
				Type = rep('Neutral', length(chrom.arms.dt$chrom.arm)))
		}
		curr.dt[,list(tumor.sample.id, chrom.arm, Type)]
	}) 
CN.oncoscan.dt <- rbindlist(CN.oncoscan.dt.ls)

CN.oncoscan.dt <- merge(
	x = clinical.data[,list(tumor.sample.id, Sex, Histology, years.dialysis)],
	y = CN.oncoscan.dt)

CN.oncoscan.dt[,Type:=factor(Type, levels=c('LOSS', 'LOH', 'Neutral', 'GAIN', 'AMP'), ordered=TRUE)]

CN.oncoscan.dt <- merge(
	x = chrom.arms.dt,
	y = CN.oncoscan.dt,
	by = c('chrom.arm'))

save(list=c('clinical.data', 'sample.ids.ls', 'CN.oncoscan.ls', 'CN.oncoscan.dt'),
	file = file.path( run.dir, "result/oncoscanR/oncoscan_arm_analysis.Rdata" ))
