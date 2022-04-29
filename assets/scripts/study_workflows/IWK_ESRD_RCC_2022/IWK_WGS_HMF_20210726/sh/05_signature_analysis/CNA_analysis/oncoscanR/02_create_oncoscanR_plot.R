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

library(data.table)
library(parallel)
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

chrom.arms.dt <- data.table(
		chrom = c(c(1:24), c(1:24)),
		chromosome = c(c(1:22, 'X', 'Y'), c(1:22, 'X', 'Y')),
		chrom.arm = c(paste(c(1:22, 'X', 'Y'), 'p', sep=""), paste(c(1:22, 'X', 'Y'), 'q', sep="")))

chrom.arms.dt <- chrom.arms.dt[!chrom.arm%in%paste(c('13', '14', '15', '22'), 'p', sep=""),]
chrom.arms.dt <- chrom.arms.dt[order(chrom, chrom.arm),]
chrom.arms.dt[,chrom.arm.f:=factor(chrom.arm, levels=chrom.arms.dt$chrom.arm, ordered=TRUE)]

load( file = file.path( run.dir, "result/oncoscanR/oncoscan_arm_analysis.Rdata" ))

clinical.data[,Sex:=ifelse(gender=='MALE', 'M', 'F')]

#CN.oncoscan.dt[,Type:=ifelse(Type=='AMP', 'GAIN', Type)]

sample.summary <- clinical.data[,list(male.ct=sum(ifelse(gender=='MALE', 1, 0)), total.ct=.N), by=list(Histology)]
All.sample.summary <- clinical.data[,list(All.male.ct=sum(ifelse(gender=='MALE', 1, 0)), All.total.ct=.N)]


CN.oncoscan.dt <- merge(
	x = CN.oncoscan.dt,
	y = sample.summary,
	by = c('Histology'))

CN.oncoscan.dt[,All.male.ct:=All.sample.summary$All.male.ct]
CN.oncoscan.dt[,All.total.ct:=All.sample.summary$All.total.ct]

CN.oncoscan.dt[,sample.ct:=ifelse(Sex=='M' & chromosome%in%c('X', 'Y'), male.ct, total.ct)]
CN.oncoscan.dt[,All.sample.ct:=ifelse(Sex=='M' & chromosome%in%c('X', 'Y'), All.male.ct, All.total.ct)]

CN.gain.summary.all <- CN.oncoscan.dt[,list(gain.ct=sum(ifelse(Type%in%c('GAIN', 'AMP'), 1, 0))), by=list(sample.ct=All.sample.ct, chromosome, chrom, chrom.arm)]
CN.loss.summary.all <- CN.oncoscan.dt[,list(loss.ct=sum(ifelse(Type%in%c('LOSS', 'LOH'), 1, 0))), by=list(sample.ct=All.sample.ct, chromosome, chrom, chrom.arm)]

CN.gain.summary <- CN.oncoscan.dt[,list(gain.ct=sum(ifelse(Type%in%c('GAIN', 'AMP'), 1, 0))), by=list(Histology, sample.ct, chromosome, chrom, chrom.arm)]
CN.loss.summary <- CN.oncoscan.dt[,list(loss.ct=sum(ifelse(Type%in%c('LOSS', 'LOH'), 1, 0))), by=list(Histology, sample.ct, chromosome, chrom, chrom.arm)]

CN.gain.summary.all[,gain:=gain.ct/sample.ct]
CN.loss.summary.all[,loss:=loss.ct/sample.ct]

CN.gain.summary[,gain:=gain.ct/sample.ct]
CN.loss.summary[,loss:=loss.ct/sample.ct]

CN.gain.summary[,max.gain:=max(gain), by=list(chrom.arm)]
CN.gain.summary[,max.gain.ct:=max(gain.ct), by=list(chrom.arm)]

CN.loss.summary[,max.loss:=max(loss), by=list(chrom.arm)]
CN.loss.summary[,max.loss.ct:=max(loss.ct), by=list(chrom.arm)]

CN.gain.summary <- merge(
	x = CN.gain.summary.all[,list(chromosome, chrom, chrom.arm, all.gain=gain, all.gain.ct=gain.ct)],
	y = CN.gain.summary)

CN.loss.summary <- merge(
		x = CN.loss.summary.all[,list(chromosome, chrom, chrom.arm, all.loss=loss, all.loss.ct=loss.ct)],
		y = CN.loss.summary)
CN.gain.summary.filt <- CN.gain.summary[max.gain>=0.2 | loss>0.2,]

CN.gain.summary.wide <- dcast(
	data = CN.gain.summary[all.gain>0.2 | max.gain>0.2,],
	formula = chromosome + chrom + chrom.arm + all.gain ~ Histology,
	value.var = c('gain'),
	fun.aggregate = max,
	fill = 0)

CN.loss.summary.wide <- dcast(
		data = CN.loss.summary[all.loss>0.2 | max.loss>0.2,],
		formula = chromosome + chrom + chrom.arm + all.loss ~ Histology,
		value.var = c('loss'),
		fun.aggregate = max,
		fill = 0)

setnames(CN.gain.summary.wide, old='ACD-RCC', new='ACD.RCC')
setnames(CN.loss.summary.wide, old='ACD-RCC', new='ACD.RCC')

CN.filtered.summary <- rbind(
	CN.gain.summary.wide[,list(CNA='Gain', chromosome, chrom.arm, All=all.gain, ccRCC, ACD.RCC, pRCC, ccpRCC, chRCC)][order(-All),],
	CN.loss.summary.wide[,list(CNA='Loss', chromosome, chrom.arm, All=all.loss, ccRCC, ACD.RCC, pRCC, ccpRCC, chRCC)][order(-All),])

library(openxlsx)

wb <- createWorkbook()
addWorksheet(wb, sheetName = "Copy number summary", gridLines = FALSE)
writeDataTable(wb, sheet = "Copy number summary",
	x = CN.filtered.summary, colNames = TRUE, rowNames = FALSE)

saveWorkbook(wb, file.path(run.dir, 'result/oncoscanR/CNA_frequencies.xlsx'), overwrite = TRUE) 

save(list=c('CN.filtered.summary', 'CN.gain.summary', 'CN.loss.summary', 'CN.oncoscan.dt'),
	file = file.path(run.dir, 'result/oncoscanR/CNA_summary.Rdata'))

wide.all.samples <- as.data.frame(dcast(
	data = CN.oncoscan.dt[,list(chrom.arm.f, tumor.sample.id, Type)],
	formula = chrom.arm.f ~ tumor.sample.id,
	fill = factor('Neutral', levels=c('LOSS', 'LOH', 'Neutral', 'GAIN'), ordered=TRUE),
#	fill = factor('Neutral', levels=c('LOSS', 'LOH', 'Neutral', 'GAIN', 'AMP'), ordered=TRUE),
	value.var = 'Type'))

mat.all.samples <- wide.all.samples[,2:ncol(wide.all.samples)]
rownames(mat.all.samples) <- wide.all.samples$chrom.arm.f

mat.all.samples <- data.matrix(mat.all.samples)

clinical.data <- clinical.data[order(Histology, tumor.sample.id),]

# prepare data for histology annotation
his <- clinical.data[,list(id=tumor.sample.id, CancerType=Histology)]
his[,ccRCC:=ifelse(CancerType=='ccRCC', TRUE, FALSE)]
his[,ACD.RCC:=ifelse(CancerType=='ACD-RCC', TRUE, FALSE)]
his[,pRCC:=ifelse(CancerType=='pRCC', TRUE, FALSE)]
his[,chRCC:=ifelse(CancerType=='chRCC', TRUE, FALSE)]
his[,ccpRCC:=ifelse(CancerType=='ccpRCC', TRUE, FALSE)]

his <- as.data.frame(his)
rownames(his) <- his$id

bw <- c("TRUE" = "black", "FALSE" = "gray80")
# prepare annotation tracks
CN.col <- list('Loss' = 'blue', 'LOH' = 'purple',
	'Neutral' = 'grey', 'Gain' = 'red')
#	'Neutral' = 'grey', 'Gain' = 'orange', 'Amp' = 'red')

col <-  colorRamp2(c(1:length(CN.col)), unlist(CN.col))
#col <-  list(c(1, 2, 3, 4, 5) = c('blue', 'purple', 'grey', 'red', 'red'))

his.col <- colorRamp2(names(Histology.colors.ls), Histology.colors.ls)

top_annotation <- HeatmapAnnotation(
		df = data.frame(
			Histology = clinical.data$Histology),
		col = list(
			Histology = Histology.colors.ls),
		annotation_name_gp = gpar(fontsize = 8),
		show_legend = TRUE,
		annotation_legend_param = list(
				Histology = list( 
						title_gp = gpar(fontsize = 10), 
						labels_gp = gpar(fontsize = 8), 
						nrow=1)),
		simple_anno_size = unit(0.25, "cm"))

h <- Heatmap(
		mat.all.samples[,clinical.data$tumor.sample.id],
		name = "CNAs",
		heatmap_legend_param = list(
			title = 'CNA',
			title_gp = gpar(fontsize=9),
			labels = c('Loss', 'LOH', 'Neutral', 'Gain'),
#			labels = c('Loss', 'LOH', 'Neutral', 'Gain', 'Amp'),
			labels_gp = gpar(fontsize=8), color_bar = "discrete"),
		col = col,
		show_row_names=T,
		row_names_side = 'left',
#		cluster_columns = cluster_columns,
		cluster_rows = FALSE,
		cluster_columns = FALSE,
		top_annotation = top_annotation,
		column_names_gp = gpar(fontsize = 7),
		row_names_gp = gpar(fontsize = 7),
		height = unit(9.0, "cm")
)

h.p <- draw(
		h,
		annotation_legend_side = "top",
		heatmap_legend_side = 'right'
#		annotation_legend_list = lgd_list,
#		padding = unit(c(10, 1, 2, 8), "mm")
)

if ( file.exists( file.path(run.dir, "result/oncoscanR/figures") , recursive=FALSE)==FALSE ){
	system(paste("mkdir -p ", file.path(run.dir, "result/oncoscanR/figures"), sep=""))
}

ggexport( h.p, width=4.5, height=4.7, filename=file.path(run.dir, "result/oncoscanR/figures/CN_armlevel_heatmap.pdf"))

pdf( file.path(run.dir, "result/oncoscanR/figures/CN_armlevel_heatmap.pdf"),
	width=7, height=6.6)
draw(
		h,
		annotation_legend_side = "top",
		heatmap_legend_side = 'right'
#		annotation_legend_list = lgd_list,
#		padding = unit(c(1, 1, 2, 8), "mm")
)
dev.off()
