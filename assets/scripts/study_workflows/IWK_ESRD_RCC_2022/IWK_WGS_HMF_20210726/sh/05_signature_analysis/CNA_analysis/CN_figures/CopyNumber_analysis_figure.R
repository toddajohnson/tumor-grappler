#!/usr/bin/env Rscript

# module use /usr/local/package/modulefiles/
# module load R/4.0.2
# . ~/.R-4.1.0_setup.sh

options(width=210)
working_dir <- getwd()

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "~/workspace/runs", study.dir )

HOME.dir <- '/home/tjohnson'

gpl_prefix <- "GRIDSS-2.12.0"

extra.samples.to.exclude <- c()
tumor.depth.cutoff <- 10;

source( file.path(run.dir, "config_files/common_config.R") )

library(data.table)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(magick)

## Extract data to process for making CNA figure
load( file = file.path(macbook.run.dir , 'result/sample_summaries/clinical_data_with_colors.Rdata'))
load(file = file.path(macbook.run.dir , 'result/CN_segments/CNA_frequency_analysis/processed_copyNumberSummaries.Rdata'))

## Count segments that cover a bin by 80% or more
seg.bin.propn.cutoff <- 0.8

no_cores <- 4
cl <- new_cluster(no_cores)

processedCopyNumbers.dt <- as.data.table(processedCopyNumbers[,1:(ncol(processedCopyNumbers)-1)])

curr.genome.level <- 'chrom.arm'
bins <- hg38_coordinate.lists.gr.ls[[curr.genome.level]]
ol = as.matrix(findOverlaps(bins, processedCopyNumbers$region, type = "any"))
olCopyNumbers = cbind(
		bin = ol[, 1],
		chromosome=seqnames(bins[ol[,1]]),
		name = elementMetadata(bins[ol[,1]])$name,
		bin.start=start(bins[ol[,1]]), bin.end=end(bins[ol[,1]]),
		processedCopyNumbers[ol[, 2], c("sampleId", "cancerType", "start", "end", "loh", "absDel", "relDel", "gain", 'amp2_0')])

olCopyNumbers.dt <- as.data.table(olCopyNumbers)
olCopyNumbers.dt[,bin.width:=bin.end-bin.start+1]
olCopyNumbers.dt[,start.in.bin:=ifelse(start<bin.start, bin.start, start)]
olCopyNumbers.dt[,end.in.bin:=ifelse(end>bin.end, bin.end, end)]
olCopyNumbers.dt[,seg.width:=end.in.bin-start.in.bin+1]

olCopyNumbers.summary <- olCopyNumbers.dt[,list(seg.total=sum(seg.width)),
	by=list(bin, chromosome, name, bin.start, bin.end, bin.width, sampleId, cancerType, loh, absDel, relDel, gain, amp2_0)]
olCopyNumbers.summary[,seg.bin.propn:=seg.total/bin.width]

binSampleSummary = olCopyNumbers.summary[seg.bin.propn>seg.bin.propn.cutoff,
				list(chr=chromosome, name, region.start=bin.start, region.end=bin.end, sampleId, cancerType, loh, absDel, relDel, gain, amp2_0)] %>% 
		group_by(chr, name, region.start, region.end, sampleId, cancerType) %>% partition(cluster = cl) %>%
		summarise(cancerType = dplyr::first(cancerType),
				LOH = any(loh), LOSS = any(absDel) | any(relDel), GAIN = any(gain), AMP = any(amp2_0))  %>%
		collect() %>%
		as_tibble()

sampleChromarmCopyNumberSummary.dt <- as.data.table(binSampleSummary)
setnames(sampleChromarmCopyNumberSummary.dt, old='name', new=curr.genome.level)

## Setup chromosome arm info to merge with CN summary
chrom.map.ls <- as.integer(c(1:24))
names(chrom.map.ls) <- paste('chr', c(1:22, 'X', 'Y'), sep="")

chrom.arms.dt <- data.table(
		chrom = c(c(1:24), c(1:24)),
		chromosome = c(c(1:22, 'X', 'Y'), c(1:22, 'X', 'Y')),
		chrom.arm = c(paste(c(1:22, 'X', 'Y'), 'p', sep=""), paste(c(1:22, 'X', 'Y'), 'q', sep="")))

chrom.arms.dt <- chrom.arms.dt[!chrom.arm%in%paste(c('13', '14', '15', '22'), 'p', sep=""),]
chrom.arms.dt <- chrom.arms.dt[order(chrom, chrom.arm),]
chrom.arms.dt[,chrom.arm.f:=factor(chrom.arm, levels=chrom.arms.dt$chrom.arm, ordered=TRUE)]

sampleChromarmCopyNumberSummary.dt <- merge(
	x = chrom.arms.dt[,list(chromosome, chrom.arm, chrom.arm.f)],
	y = sampleChromarmCopyNumberSummary.dt[,
		list(chr, chrom.arm, region.start, region.end, tumor.sample.id=sampleId, cancerType,
			LOSS, LOH, GAIN, AMP,
			Type = ifelse(LOSS==TRUE, 'LOSS', ifelse(LOH==TRUE, 'LOH',
				ifelse(GAIN==TRUE, 'GAIN', ifelse(AMP==TRUE, 'AMP', 'Neutral')))))],
	by = c('chrom.arm'))

sampleChromarmCopyNumberSummary.dt[,Type:=factor(Type, levels=c('LOSS', 'LOH', 'Neutral', 'GAIN'), ordered=TRUE)]

# Setup up CN summary into matrix
wide.all.samples <- as.data.frame(dcast(
				data = sampleChromarmCopyNumberSummary.dt[,list(chrom.arm.f, tumor.sample.id, Type)],
				formula = chrom.arm.f ~ tumor.sample.id,
				fill = factor('Neutral', levels=c('LOSS', 'LOH', 'Neutral', 'GAIN'), ordered=TRUE),
				value.var = 'Type'))

mat.all.samples <- wide.all.samples[,2:ncol(wide.all.samples)]
rownames(mat.all.samples) <- wide.all.samples$chrom.arm.f

mat.all.samples <- data.matrix(mat.all.samples)

clinical.data <- clinical.data[order(Histology, tumor.sample.id),]

# prepare annotation tracks
CN.col <- list('Loss' = 'blue', 'LOH' = 'purple',
		'Neutral' = 'grey', 'Gain' = 'red')

CN.col.ramp <-  colorRamp2(c(1:length(CN.col)), unlist(CN.col))

his.col <- colorRamp2(names(Histology.colors.ls), Histology.colors.ls)

CNA.chrom.arm.heatmap.top_annotation <- HeatmapAnnotation(
		df = data.frame(
				Histology = clinical.data[order(Histology, tumor.sample.id),]$Histology),
		col = list(
				Histology = Histology.colors.ls),
		annotation_name_gp = gpar(fontsize = 7),
		show_legend = FALSE,
		annotation_legend_param = list(
				Histology = list( 
						title = '',
						title_gp = gpar(fontsize = 7), 
						labels_gp = gpar(fontsize = 7), 
						nrow=1)),
		simple_anno_size = unit(0.25, "cm"))

CNA.chrom.arm.heatmap.h <- Heatmap(
		mat.all.samples[,clinical.data[order(Histology, tumor.sample.id),]$tumor.sample.id],
		name = "CNAs",
		heatmap_legend_param = list(
				title = 'CNA',
				title_gp = gpar(fontsize=8),
				labels = c('Loss', 'LOH', 'Neutral', 'Gain'),
				labels_gp = gpar(fontsize=7),
				color_bar = "discrete"),
		col = CN.col.ramp,
		show_row_names=T,
		row_names_side = 'left',
		cluster_rows = FALSE,
		cluster_columns = FALSE,
		top_annotation = CNA.chrom.arm.heatmap.top_annotation,
		column_names_gp = gpar(fontsize = 6),
		row_names_gp = gpar(fontsize = 6),
		height = unit(8.2, "cm")
)

CNA.chrom.arm.heatmap.h.p <- draw(
		CNA.chrom.arm.heatmap.h,
		annotation_legend_side = "top",
		align_annotation_legend = 'heatmap_center',
		heatmap_legend_side = 'right',
)

ggexport( h.p, width=3.8, height=4.0, filename=file.path(run.dir, "result/CN_segments/CNA_frequency_analysis/CN_armlevel_heatmap.pdf"))


unsupp_corr.dt <- fread( file.path(run.dir, 'result/CNApp/figures/AllSamples/Unsup_corr_mat_Cytobands_by_Histology_2022-02-02 17_23_32.tsv'))

setkey(clinical.data, 'tumor.sample.id')

unsupp_corr_mat <- as.matrix(unsupp_corr.dt[,2:ncol(unsupp_corr.dt)])
rownames(unsupp_corr_mat) <- unsupp_corr.dt$V1
corr.top_annotation <- HeatmapAnnotation(
		df = data.frame(
				Histology = clinical.data[colnames(unsupp_corr_mat),]$Histology),
		col = list(
				Histology = Histology.colors.ls),
		annotation_name_gp = gpar(fontsize = 8),
		show_legend = TRUE,
		annotation_legend_param = list(
				Histology = list(
						title_gp = gpar(fontsize = 8), 
						labels_gp = gpar(fontsize = 6), 
						nrow=1)),
		simple_anno_size = unit(0.25, "cm"))

corr.h <- Heatmap(
		unsupp_corr_mat[nrow(unsupp_corr_mat):1,],
		name = "CNA corr.",
		heatmap_legend_param = list(
				title = "Pearson's r",
				title_gp = gpar(fontsize=8),
				labels_gp = gpar(fontsize=6)),
		show_row_names=T,
		row_names_side = 'left',
		cluster_rows = FALSE,
		cluster_columns = FALSE,
		top_annotation = corr.top_annotation,
		column_names_gp = gpar(fontsize = 6),
		row_names_gp = gpar(fontsize = 6),
		height = unit(7.2, "cm")
)

corr.h.p <- draw(
		corr.h,
		annotation_legend_side = "top",
)
ggexport( corr.h.p, width=4.0, height=4.0, filename=file.path(run.dir, "result/CN_segments/CNA_frequency_analysis/CN_correlation_heatmap.pdf"))

# Produce Figure 2 from panels
library(magick)

macbook.run.dir <- '/Users/tajohnson/HGC_mounts/HGC_tjohnson/workspace/runs/IWK_WGS_HMF_20210726'

panelA.file <- file.path(macbook.run.dir, "result/CN_segments/CNA_frequency_analysis/CN_armlevel_heatmap.pdf")
panelB.file <- file.path(macbook.run.dir, "result/CN_segments/CNA_frequency_analysis/CN_correlation_heatmap.pdf")
panelC.file <- file.path(macbook.run.dir, 'result/CNA_circos/plots/All.png')
panelD.file <- file.path(macbook.run.dir, 'result/CNA_circos/plots/ACD-RCC.png')

Figure1.file <- file.path(macbook.run.dir, "result/CN_segments/CNA_frequency_analysis/Figure_2_CNA.pdf")

panelA <- image_read_pdf(panelA.file) %>%
		image_annotate(panelA, "A", size = 10, color = "black", gravity = 'northwest', location="+0+0")

panelB <- image_read_pdf(panelB.file) %>%
		image_annotate(panelB, "B", size = 10, color = "black", gravity = 'northwest', location="+0+0")

panelC <- image_read(panelC.file) %>%
		image_annotate(panelC, "C", size = 120, color = "black", gravity = 'northwest', location="+50+0") %>%
		image_annotate("All", size = 120, color = "black", gravity = 'northwest', location="+250+150")

panelD <- image_read(panelD.file) %>%
		image_annotate(panelD, "ACD-RCC", size = 120, color = "black", gravity = 'northwest', location="+0+150")

pdf( width=8, height=8, file=Figure1.file)
#quartz(width=7.5, height=7.5)
layout(matrix(1:4, nrow=2, byrow=TRUE))
par(mai=c(0, 0.1, 0.1, 0))
plot(panelA)
par(mai=c(0, 0, 0.1, 0))
plot(panelB)
par(mai=c(0, 0.1, 0, 0))
plot(panelC)
plot(panelD)
dev.off()

rm(list=c('panelA', 'panelB', 'panelC', 'panelD'))
# Produce Figure S1 from panels
load( file = file.path(macbook.run.dir , 'result/sample_summaries/clinical_data_with_colors.Rdata'))

FigureS1.file <- file.path(macbook.run.dir, "result/CN_segments/CNA_frequency_analysis/Figure_S1_CNA.pdf")

Histology.ls <- c('All', names(Histology.colors.ls))
names(Histology.ls) <- Histology.ls

FigureS2.panels.ls <- lapply(
	X = Histology.ls,
	FUN = function(curr.histology){
		curr.file <- file.path(macbook.run.dir, paste('result/CNA_circos/plots/', curr.histology, '.png', sep=""))
		image_read(curr.file) %>%
			image_annotate(curr.histology, size = 120, color = "black", gravity = 'northwest', location="+0+150")
	})

panel.ct <- length(Histology.ls)
if(panel.ct%%2 == 1){
	plotmat <- matrix(1:(panel.ct+1), ncol=2, byrow=TRUE)
}else{
	plotmat <- matrix(1:(panel.ct), ncol=2, byrow=TRUE)
}

pdf( width=8.5, height=11, file=FigureS1.file)
#quartz(width=8.5, height=11)
layout(plotmat)
par(mai=c(0, 0, 0, 0))
lapply(
	X = Histology.ls,
	FUN = function(curr.histology){
		plot(FigureS2.panels.ls[[curr.histology]])
	})

dev.off()