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


library(data.table)
library(parallel)
library(tidyverse)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(edgeR)
library(limma)
library(ggpubr)

source( file.path(run.dir, "config_files/common_config.R") )

#ref.dir <- "/home/tjohnson/reference/HMF/38/dbs/ensembl_data_cache"

#HOME.dir <- "~/HGC_mounts/HGC/"
HOME.dir <- "/home/tjohnson"

isofox.dir <- file.path(run.dir, 'result', 'isofox')

load(file = paste(isofox.dir, "/isofox.gene_expression.edgeR.normalized.Rdata", sep=""))

load(file=file.path(isofox.dir, 'imported.tx.gene.expression.Rdata'))
load(file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))

load(file = paste(wgs.run.dir, "/result/sample_summaries/clinical_data_with_colors.Rdata", sep=""))

txi <- copy(imported.gene.expr.data)

adjTPM <- txi$abundance
# adjTPM should sum to close to 1 million
#colSums(adjTPM)
##	IWK001_T IWK002_T IWK005_T IWK006_T IWK008_T IWK009_T IWK010_T IWK012_T IWK013_T IWK015_T IWK016_T IWK017_T IWK018_T IWK019_T IWK020_T IWK021_T IWK022_T IWK024_T IWK025_T IWK026_T IWK028_T IWK029_T IWK039_T IWK040_T IWK048_T IWK049_T 
##	1006125  1004020  1008519  1059575  1049016  1003132  1000226  1000241  1026875  1002011  1012465  1000152  1069968  1003368  1001236  1001430  1000358  1009295  1003912  1017387  1029253  1001414  1079054  1009931  1012630  1022590 

cts <- txi$counts
normMat <- txi$length

setkey(clinical.data, tumor.sample.id)

sample.ids <- names(cls)

sample.info.with.clinical.info.dt <- clinical.data[sample.ids,]

sample.info.with.clinical.info.dt[names(cls),cls3:=cls]
sample.info.with.clinical.info.dt[names(cls8),cls8:=cls8]

sample.info <- sample.info.with.clinical.info.dt[cls3%in%c('B', 'C'),]
clusters.B.C.sample.ct <- nrow(sample.info)
sample.info[,exprClass:=ifelse(cls3=='B', 'ccRCC', ifelse(cls3=='C', 'non-ccRCC', 'UNK'))]
samples.to.keep <- sample.info$tumor.sample.id

## > nrow(cts)
# 39357
## [1] 37640
## > ncol(cts)
## [1] 26

cols.to.keep <- colnames(cts)
cols.to.keep <- cols.to.keep[which(cols.to.keep%in%samples.to.keep)]

cts <- cts[,cols.to.keep]
normMat <- normMat[,cols.to.keep]


# Obtaining per-observation scaling factors for length, adjusted to avoid
# changing the magnitude of the counts.
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat

# Computing effective library sizes from scaled counts, to account for
# composition biases between samples.
eff.lib <- calcNormFactors(normCts) * colSums(normCts)

# Combining effective library sizes with the length factors, and calculating
# offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

# Creating a DGEList object for use in edgeR.

groups <- sample.info$exprClass
names(groups) <- sample.info$tumor.sample.id
groups <- groups[colnames(cts)]
groups.f <- factor(x = groups, levels=c('non-ccRCC', 'ccRCC'), ordered = TRUE)

y <- DGEList(cts, group=groups.f)

y <- scaleOffset(y, normMat)
# filtering

keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]

pdf(file=file.path(run.dir, 'result/DGE/ccRCC_non_ccRCC/MDS.pdf'))
plotMDS(y)
dev.off()

#
# TMM normalization and log transformation
#
y <- estimateDisp(y)

pdf(file=file.path(run.dir, 'result/DGE/ccRCC_non_ccRCC/BCV.pdf'))
plotBCV(y)
dev.off()

fit <- glmFit(y)
lrt <- glmLRT(fit)
topTags(lrt)
## > topTags(lrt)

## Coefficient:  y$samples$group.L 
##         logFC     logCPM        LR       PValue          FDR
#HIF1A-AS3  8.956064  5.4543828 259.68214 2.012811e-58 3.484176e-54
#VEGFA      3.527579 11.6377215 149.73424 1.981746e-34 1.715201e-30
#SYT14     -6.239955  0.9476457 148.61441 3.481974e-34 2.009099e-30
#CCND1      2.226044  9.0551956 139.56530 3.313403e-32 1.433875e-28
#PCSK6      3.015586  5.9818297 124.35032 7.060939e-29 2.444497e-25
#NDUFA4L2   4.769912 10.1462471 113.51969 1.660243e-26 4.789801e-23
#EHD2       1.821105  7.7257455 106.78307 4.966645e-25 1.228180e-21
#CP         5.643000  7.8692890 104.50537 1.567679e-24 3.392066e-21
#SLC47A2   -5.288879  7.3949468  97.96533 4.257731e-23 8.189036e-20
#PAPPA-AS1 -4.524459  0.6435428  92.72429 6.011126e-22 1.011253e-18


o <- order(lrt$table$PValue)
cpm(y)[o[1:10],]
## 
#IWK001_T    IWK002_T    IWK005_T   IWK008_T     IWK009_T   IWK010_T    IWK012_T    IWK013_T     IWK015_T   IWK016_T     IWK017_T  IWK018_T    IWK019_T     IWK020_T     IWK021_T    IWK022_T
#HIF1A-AS3   0.0000000   49.551638   97.782189   0.000000   69.4825046   88.37873   66.336635    0.000000 5.488378e+01   0.000000 1.563938e+01  0.000000   0.0000000 1.129642e+01 6.386854e+01   85.176387
#VEGFA     166.2938025 5421.939155 3494.941531 215.431437 3597.9641317 8160.85140 3268.140397   83.751982 8.169195e+03 195.123610 2.186147e+03 58.036187  99.0346251 3.122734e+03 1.289871e+04 4791.269352
#SYT14       7.9670376    0.000000    0.000000   3.764529    0.0000000    0.00000    0.000000   15.567404 0.000000e+00   1.652797 4.261412e-02  1.102587   2.3633631 0.000000e+00 0.000000e+00    0.000000
#CCND1      83.0020460 1141.431272 1289.452155  70.622561 1043.6525918  658.22966  723.407925  116.674278 6.689758e+02  87.161657 6.596666e+02 38.139498  86.7321872 4.963830e+02 8.251098e+02  911.114185
#PCSK6       3.0419598  130.994059   71.615124  11.544555   91.8301079   41.57240   43.397799    6.044963 6.905145e+01   3.367964 2.914806e+01  4.510585   2.9461101 1.905886e+02 1.717277e+02  160.951799
#NDUFA4L2    9.7053003 2295.405127 3642.160951   2.459492 1356.1286869 1764.03554 2669.679981    6.109962 1.175885e+03   0.623697 9.374681e+02 41.497378  14.4391771 2.066321e+03 1.598208e+03 1531.649533
#EHD2       46.7882389  297.197341  301.039968  74.939221  364.4415918  293.64869  385.477372   21.547367 2.850446e+02  76.059850 2.358265e+02 67.207709  51.8321134 1.913358e+02 2.826756e+02  375.364996
#CP          0.3621381  493.660307  321.460726   2.258717  293.9344155   89.47537  137.156121    2.177487 1.770321e+02   3.523888 5.671940e+01  1.603763   0.8417458 2.517738e+02 6.708688e+02  437.163108
#SLC47A2    27.8846315    9.786589    1.804625 677.866159    0.3903511    0.99694    3.719811 1443.023568 2.042187e+00  26.538308 6.051205e+00 70.415236 204.0585950 1.230738e+00 1.992778e-01    4.257046
#PAPPA-AS1   0.7242761    0.000000    0.000000   7.930607    0.0000000    0.00000    0.000000    4.159974 6.381835e-02   2.245309 0.000000e+00  1.503528   5.0180997 4.395493e-02 3.487362e-01    0.000000
#IWK024_T     IWK025_T    IWK026_T     IWK028_T     IWK029_T     IWK048_T
#HIF1A-AS3   38.315310   92.7211586   0.0000000   0.05491907   48.8914957 1.741094e+02
#VEGFA     2110.923595 3078.7835369 189.9364388 245.76285257 4567.9953449 4.773824e+03
#SYT14        0.000000    0.0000000   3.6229487   4.94271659    0.0000000 3.787456e-02
#CCND1      723.887003  746.8991361  84.5019232 140.37315110  642.8445418 4.412765e+02
#PCSK6       95.098076  115.0624512   4.0254985   4.94271659   57.5761693 9.184582e+01
#NDUFA4L2  1244.035052 1428.1647323  22.3415168  35.09328778 1642.0108762 1.488129e+03
#EHD2       293.986991  438.8673658  38.8796066  29.38170416  231.6627742 2.653492e+02
#CP          42.083417  251.5073412   0.4360957   0.27459537  504.3901161 1.388444e+03
#SLC47A2      4.439651    0.5753123 515.3309030 826.47713260    0.6433092 9.847387e-01
#PAPPA-AS1    0.000000    0.0000000   6.4072518   4.33860678    0.1072182 3.787456e-02

## cpm(y)['HIF1A-AS3',][order(cpm(y)['HIF1A-AS3',])]
#IWK001_T     IWK008_T     IWK013_T     IWK016_T     IWK018_T     IWK019_T     IWK026_T     IWK028_T     IWK020_T     IWK017_T     IWK024_T     IWK029_T     IWK002_T     IWK015_T     IWK021_T     IWK012_T 
#0.00000000   0.00000000   0.00000000   0.00000000   0.00000000   0.00000000   0.00000000   0.05491907  11.29641709  15.63938267  38.31530959  48.89149570  49.55163846  54.88377839  63.86854390  66.33663548 
#IWK009_T     IWK022_T     IWK010_T     IWK025_T     IWK005_T     IWK048_T 
#69.48250459  85.17638746  88.37873376  92.72115856  97.78218869 174.10937043 

mean(cpm(y)['HIF1A',][which(cpm(y)['HIF1A-AS3',]==0)])
#[1] 123.3643
mean(cpm(y)['HIF1A',][which(cpm(y)['HIF1A-AS3',]>0)])
#[1] 62.14457
#> cpm(y)['HIF1A',][which(cpm(y)['HIF1A-AS3',]==0)]
#IWK001_T  IWK008_T  IWK013_T  IWK016_T  IWK018_T  IWK019_T  IWK026_T 
#32.15533 121.02523 211.37023 284.42258  76.92739  61.47837  86.93637

summary(decideTests(lrt))
## y$samples$group.L
#Down                2246
#NotSig             13054
#Up                  2010


pdf(file=file.path(run.dir, 'result/DGE/ccRCC_non_ccRCC/MD.pdf'))
plotMD(lrt)
abline(h=c(-1, 1), col="blue")
dev.off()

library(EnhancedVolcano)

# Import edgeR data
## Coefficient:  y$samples$group.L 
## logFC     logCPM        LR       PValue          FDR

volcanodata <- as.data.table(topTags(lrt, n=nrow(lrt$table)), keep.rownames=TRUE)
setnames(volcanodata,
	old=paste('table.', c('rn', 'logFC', 'logCPM', 'LR', 'PValue', 'FDR'), sep=""),
	new = c("Gene", "logFC",  "logCPM", "LR", "pvalue", "FDR"))

pdf(file=file.path(run.dir, 'result/DGE/ccRCC_non_ccRCC/ccRCC_nonccRCC_Volcano.pdf'))
EnhancedVolcano(
		volcanodata, 
		x = 'logFC',
		y = 'pvalue',
		lab = volcanodata$Gene,
#		AdjustedCutoff = 0.05, 
#		LabellingCutoff = 0.05, 
#		FCCutoff = 2.0,
		title = "ccRCC vs. non-ccRCC edgeR")
#		col = c("grey30", "forestgreen", "royalblue", "red2"))

dev.off()

curr.p <- EnhancedVolcano(
				volcanodata, 
				x = 'logFC',
				y = 'pvalue',
				lab = volcanodata$Gene,
				title = NULL,
				subtitle = NULL,
				axisLabSize = 8,
				captionLabSize = 8,
				caption = NULL,
				legendPosition = 'none',
#				pCutoff = 0.001,
#				FCcutoff = 1,
#		legendLabSize = 8,
#		legendIconSize = 3,
				borderWidth = 0.5,
				pointSize = 0.8,
				labSize = 2) +
		theme(plot.margin = margin(0.1, 0.1, 0.05, 0.05, "cm")) +
		theme( axis.ticks = element_line(size = 0.5),
				axis.line = element_line(size = 0.5)) +
		theme( panel.grid.minor = element_line(size = 0.5),
				panel.grid.major = element_line(size = 0.5)) +
		theme(axis.title.y = element_text(margin = margin(t = 0, r = 0.6, b = 0, l = 0))) +
		theme(axis.title.x = element_text(margin = margin(t = 1, r = 0, b = 0, l = 0))) +
		theme(axis.text.y = element_text(margin = margin(t = 0, r = 1, b = 0, l = 0))) +
		theme(axis.text.x = element_text(margin = margin(t = 1, r = 0, b = 0, l = 0)))

ggexport( curr.p, width=2.2, height=2.2, filename=file.path(run.dir, 'result/DGE/ccRCC_non_ccRCC/ccRCC_nonccRCC_Volcano.for.figure.2.2in.pdf'))
ggexport( curr.p, width=2.0, height=2.0, filename=file.path(run.dir, 'result/DGE/ccRCC_non_ccRCC/ccRCC_nonccRCC_Volcano.for.figure.2.0.in.pdf'))


pdf(file=file.path(run.dir, 'result/DGE/ccRCC_non_ccRCC/ccRCC_nonccRCC_Volcano_byFDR.pdf'))
EnhancedVolcano(
		volcanodata, 
		x = 'logFC',
		y = 'FDR',
		lab = volcanodata$Gene,
#		AdjustedCutoff = 0.05, 
#		LabellingCutoff = 0.05, 
		FCcutoff = 2.0,
		title = "ccRCC vs. non-ccRCC edgeR")
#		col = c("grey30", "forestgreen", "royalblue", "red2"))
dev.off()


save( list=c('volcanodata', 'sample.info', 'lrt', 'y'),
		file=file.path(run.dir, 'result/DGE/ccRCC_non_ccRCC/ccRCC_nonccRCC_DE.Rdata'))

fwrite(
		volcanodata[,list(Gene, logFC, logCPM, LR, pvalue, FDR)],
		file = file.path(run.dir, 'result/DGE/ccRCC_non_ccRCC/ccRCC_nonccRCC_DGE_list.tsv'),
		sep = '\t')

extracted.cpm <- cpm(y, log=FALSE)[c('HIF1A', 'HIF1A-AS3', 'VHL', 'EPAS1', 'CDH1', 'SLC22A7', 'RPL37'),]
extracted.log.cpm <- cpm(y, log=TRUE)[c('HIF1A', 'HIF1A-AS3', 'VHL', 'EPAS1', 'CDH1', 'SLC22A7', 'RPL37'),]

volcanodata[Gene%in%c('HIF1A', 'HIF1A-AS3', 'VHL', 'EPAS1'),]
## > volcanodata[Gene%in%c('HIF1A', 'HIF1A-AS2', 'VHL', 'EPAS1'),]
## Gene      logFC   logCPM         LR       pvalue          FDR adjust.method        comparison test
## 1: HIF1A-AS3  8.9560637 5.454383 259.682136 2.012811e-58 3.484176e-54            BH y$samples$group.L  glm
## 2:     EPAS1  1.3994533 8.895595  35.240398 2.914154e-09 2.025863e-07            BH y$samples$group.L  glm
## 3:       VHL -0.5664421 5.016849  11.533075 6.836885e-04 5.441896e-03            BH y$samples$group.L  glm
## 4:     HIF1A -0.6208408 6.326212   5.414238 1.997317e-02 7.201323e-02            BH y$samples$group.L  glm

extracted.cpm.dt <- as.data.table(extracted.cpm, keep.rownames=TRUE)
setnames(extracted.cpm.dt, old='rn', new='gene')

extracted.cpm.dt.long <- melt(
		data = extracted.cpm.dt,
		id.vars = c('gene'),
		measure.vars = colnames(extracted.cpm.dt[,-c('gene')]),
		variable.name = 'tumor.sample.id',
		value.name = 'ct')

extracted.log.cpm.dt <- as.data.table(extracted.log.cpm, keep.rownames=TRUE)
setnames(extracted.log.cpm.dt, old='rn', new='gene')

extracted.log.cpm.dt.long <- melt(
		data = extracted.log.cpm.dt,
		id.vars = c('gene'),
		measure.vars = colnames(extracted.log.cpm.dt[,-c('gene')]),
		variable.name = 'tumor.sample.id',
		value.name = 'log2.ct')

extracted.cpm.dt.long <- merge(
		x = data.table(
				tumor.sample.id = rownames(y$samples),
				exprClass = y$samples$group),
		y = extracted.cpm.dt.long,
		by = c('tumor.sample.id'))

extracted.cpm.dt.long <- merge(
		x = extracted.cpm.dt.long,
		y = extracted.log.cpm.dt.long,
		by = c('tumor.sample.id', 'gene'))

extracted.cpm.dt.by.gene <- dcast(
		extracted.cpm.dt.long,
		formula = exprClass + tumor.sample.id  ~ gene,
		value.var = c('ct', 'log2.ct'))

extracted.cpm.dt.by.gene[,Histology:=exprClass]

exprClass.colors <- c('brown2', 'aquamarine4')
names(exprClass.colors) <- c("ccRCC", "non-ccRCC")

extracted.cpm.dt.long[,exprClass.color:=exprClass.colors[exprClass]]

my_comparisons1 <- list(
		c("ccRCC", "non-ccRCC"))

curr.alpha <- 1.0
jitter.size <- 1.0

VHL.boxplot.p <- ggboxplot(
				extracted.cpm.dt.long[gene=='VHL',list(tumor.sample.id, ct, exprClass, exprClass.color)],
				x = "exprClass", y = "ct", 
#				outlier.shape = NA,
#				shape = 21,
#				color = "exprClass",
#				fill = "exprClass",
#				size = 0.5,
#				ylim = c(-2, 200),
#				title = 'A',
				add = "jitter",
				add.params = list(shape = 16, size = jitter.size, color = "exprClass"),
				ylab = "VHL CPM", xlab = "RCC class") +
		scale_color_manual(values = alpha(exprClass.colors, curr.alpha)) +
#		stat_compare_means(comparisons = my_comparisons1, method="wilcox.test", vjust=0) +
		rremove('legend') +
		rotate_x_text(45) +
		theme(plot.margin = margin(0.3, 0.1, 0.1, 0.1, "cm")) +
		theme(
				title = element_text(size=10, face='bold', color='black'),
				legend.key.size = unit(0.6,"line"),
				legend.text = element_text(size=8),
				axis.title = element_text(size=8), axis.text = element_text(size=8))

HIF1A_AS3.boxplot.p <- ggboxplot(
				extracted.cpm.dt.long[gene=='HIF1A-AS3',list(tumor.sample.id, ct, exprClass, exprClass.color)],
				x = "exprClass", y = "ct", 
#				outlier.shape = NA,
#				shape = 21,
#				color = "exprClass",
#				fill = "exprClass",
#				size = 0.5,
#				ylim = c(-2, 200),
#				title = 'A',
				add = "jitter",
				add.params = list(shape = 16, size = jitter.size, color = "exprClass"),
				ylab = "HIF1A-AS3 CPM", xlab = "RCC class") +
		scale_color_manual(values = alpha(exprClass.colors, curr.alpha)) +
#		stat_compare_means(comparisons = my_comparisons1, method="wilcox.test", vjust=0) +
		rremove('legend') +
		rotate_x_text(45) +
		theme(plot.margin = margin(0.3, 0.1, 0.1, 0.1, "cm")) +
		theme(
				title = element_text(size=10, face='bold', color='black'),
				legend.key.size = unit(0.6,"line"),
				legend.text = element_text(size=8),
				axis.title = element_text(size=8), axis.text = element_text(size=8))

HIF1A.boxplot.p <- ggboxplot(
				extracted.cpm.dt.long[gene=='HIF1A',list(tumor.sample.id, ct, exprClass)],
				x = "exprClass", y = "ct", 
#				outlier.shape = NA,
#				title = 'B',
#				color = 'black',
#				fill = "exprClass",
				add = "jitter",
				add.params = list(shape = 16, size = jitter.size, color = "exprClass"),
				ylab = "HIF1A CPM", xlab = "RCC class") +
		scale_color_manual(values = alpha(exprClass.colors, curr.alpha)) +
#		stat_compare_means(comparisons = my_comparisons1, method="t.test", vjust=0.4)+
		rremove('legend') +
		rotate_x_text(45) +
		theme(plot.margin = margin(0.3, 0.1, 0.1, 0.1, "cm")) +
		theme(
				title = element_text(size=10, face='bold', color='black'),
				legend.key.size = unit(0.6,"line"),
				legend.text = element_text(size=8),
				axis.title = element_text(size=8), axis.text = element_text(size=8))


EPAS1.boxplot.p <- ggboxplot(
				extracted.cpm.dt.long[gene=='EPAS1',list(tumor.sample.id, ct, exprClass)],
				x = "exprClass", y = "ct", 
#				outlier.shape = NA,
#				title = 'B',
#				color = 'black',
#				fill = "exprClass",
				add = "jitter",
				add.params = list(shape = 16, size = jitter.size, color = "exprClass"),
				ylab = "EPAS1 CPM", xlab = "RCC class") +
		scale_color_manual(values = alpha(exprClass.colors, curr.alpha)) +
#		stat_compare_means(comparisons = my_comparisons1, method="t.test", vjust=0.4)+
		rremove('legend') +
		rotate_x_text(45) +
		theme(plot.margin = margin(0.3, 0.1, 0.1, 0.1, "cm")) +
		theme(
				title = element_text(size=10, face='bold', color='black'),
				legend.key.size = unit(0.6,"line"),
				legend.text = element_text(size=8),
				axis.title = element_text(size=8), axis.text = element_text(size=8))

VHL.HIF1A_AS3.scatter.p <- ggscatter(
				data = extracted.cpm.dt.by.gene[order(-exprClass),],
				x = 'ct_HIF1A-AS3',
				y  = "ct_VHL",
				shape = 16,
#				color = 'black',
				color = "Histology",
#				title = 'C',
				add = "reg.line", 
				add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
				xlab = "HIF1A-AS3 CPM",
				ylab = "VHL CPM") +
		stat_cor(method = "pearson", label.x = 15, label.y = 45, size=2.5) +
		rremove("legend") +
		scale_color_manual(values = alpha(exprClass.colors, curr.alpha)) +
		theme(title = element_text(size=10, face='bold', color='black'),
				legend.key.size = unit(0.6,"line"),
				legend.text = element_text(size=8),
				axis.title = element_text(size=8), axis.text = element_text(size=8))


HIF1A.EPAS1.scatter.p <- ggscatter(
				data = extracted.cpm.dt.by.gene[order(-exprClass),],
				x = 'ct_HIF1A',
				y  = "ct_EPAS1",
				shape = 16,
#				color = 'black',
				color = "Histology",
#				title = 'E',
				add = "reg.line", 
				add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
				xlab = "HIF1A CPM",
				ylab = "EPAS1 CPM") +
		stat_cor(method = "pearson", label.x = 50, label.y = 550, size=2.5) +
		rremove("legend") +
		scale_color_manual(values = alpha(exprClass.colors, curr.alpha)) +
		theme(title = element_text(size=10, face='bold', color='black'),
				legend.key.size = unit(0.6,"line"),
				legend.text = element_text(size=8),
				axis.title = element_text(size=8), axis.text = element_text(size=8))

EPAS1.HIF1A_AS3.scatter.p <- ggscatter(
				data = extracted.cpm.dt.by.gene[order(-exprClass),],
				x = 'ct_HIF1A-AS3',
				y  = "ct_EPAS1",
				shape = 16,
#				color = 'black',
				color = "Histology",
#				title = 'D',
				add = "reg.line", 
				add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
				xlab = "HIF1A-AS3 CPM",
				ylab = "EPAS1 CPM") +
		stat_cor(method = "pearson", label.x = 35, label.y = 250, size=2.5) +
		rremove("legend") +
		scale_color_manual(values = alpha(exprClass.colors, curr.alpha)) +
		theme(title = element_text(size=10, face='bold', color='black'),
				legend.key.size = unit(0.6,"line"),
				legend.text = element_text(size=8),
				axis.title = element_text(size=8), axis.text = element_text(size=8))

HIF1A.HIF1A_AS3.scatter.p <- ggscatter(
				data = extracted.cpm.dt.by.gene[order(-exprClass),],
				x = 'ct_HIF1A-AS3',
				y  = "ct_HIF1A",
				shape = 16,
#				color = 'black',
				color = "Histology",
#				title = 'C',
				add = "reg.line", 
				add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
				xlab = "HIF1A-AS3 CPM",
				ylab = "HIF1A CPM") +
		stat_cor(method = "pearson", label.x = 25, label.y = 130, size=2.5) +
		rremove("legend") +
		scale_color_manual(values = alpha(exprClass.colors, curr.alpha)) +
		theme(title = element_text(size=10, face='bold', color='black'),
				legend.key.size = unit(0.6,"line"),
				legend.text = element_text(size=8),
				axis.title = element_text(size=8), axis.text = element_text(size=8))


ccRCC.HIF.fig.titles.A.to.H <- ggarrange(
		ggarrange(
				HIF1A_AS3.boxplot.p + ggtitle('A'),
				HIF1A.boxplot.p + ggtitle('B'),
				EPAS1.boxplot.p + ggtitle('C'),
				VHL.boxplot.p + ggtitle('D'),
				ncol=4, nrow=1),
		ggarrange(
				HIF1A.HIF1A_AS3.scatter.p + ggtitle('E'),
				HIF1A.EPAS1.scatter.p + ggtitle('F'),
				EPAS1.HIF1A_AS3.scatter.p + ggtitle('G'),
				VHL.HIF1A_AS3.scatter.p + ggtitle('H'),
				ncol = 4, nrow=1),
		nrow = 2, heights = c(1, 0.8))

ggexport(ccRCC.HIF.fig.titles.A.to.H, width=7, height=5,
	filename=file.path(run.dir, "result/DGE/ccRCC_non_ccRCC/ccRCC.nonccRCC_DGE.figure.pdf"))

ccRCC.HIF.fig.no.titles <- ggarrange(
		ggarrange(
				HIF1A_AS3.boxplot.p,
				HIF1A.boxplot.p,
				EPAS1.boxplot.p,
				VHL.boxplot.p ,
				ncol=4, nrow=1),
		ggarrange(
				HIF1A.HIF1A_AS3.scatter.p,
				HIF1A.EPAS1.scatter.p,
				EPAS1.HIF1A_AS3.scatter.p,
				VHL.HIF1A_AS3.scatter.p,
				ncol = 4, nrow=1),
		nrow = 2, heights = c(1, 0.8))

ggexport(ccRCC.HIF.fig.no.titles, width=7, height=5,
filename=file.path(run.dir, "result/DGE/ccRCC_non_ccRCC/ccRCC.nonccRCC_DGE.no.titles.figure.pdf"))


## log2ct plots

VHL.boxplot.p <- ggboxplot(
				extracted.cpm.dt.long[gene=='VHL',list(tumor.sample.id, ct, log2.ct, exprClass, exprClass.color)],
				x = "exprClass", y = "log2.ct", 
#				outlier.shape = NA,
#				shape = 21,
#				color = "exprClass",
#				fill = "exprClass",
#				size = 0.5,
#				ylim = c(-2, 200),
#				title = 'A',
				add = "jitter",
				add.params = list(shape = 16, size = jitter.size, color = "exprClass"),
				ylab = "VHL log2(cpm)", xlab = "RCC class") +
		scale_color_manual(values = alpha(exprClass.colors, curr.alpha)) +
#		stat_compare_means(comparisons = my_comparisons1, method="wilcox.test", vjust=0) +
		rremove('legend') +
		rotate_x_text(45) +
		theme(plot.margin = margin(0.3, 0.1, 0.1, 0.1, "cm")) +
		theme(
				title = element_text(size=10, face='bold', color='black'),
				legend.key.size = unit(0.6,"line"),
				legend.text = element_text(size=8),
				axis.title = element_text(size=8), axis.text = element_text(size=8))

HIF1A_AS3.boxplot.p <- ggboxplot(
				extracted.cpm.dt.long[gene=='HIF1A-AS3',list(tumor.sample.id, ct, log2.ct, exprClass, exprClass.color)],
				x = "exprClass", y = "log2.ct", 
#				outlier.shape = NA,
#				shape = 21,
#				color = "exprClass",
#				fill = "exprClass",
#				size = 0.5,
#				ylim = c(-2, 200),
#				title = 'A',
				add = "jitter",
				add.params = list(shape = 16, size = jitter.size, color = "exprClass"),
				ylab = "HIF1A-AS3 log2(cpm)", xlab = "RCC class") +
		scale_color_manual(values = alpha(exprClass.colors, curr.alpha)) +
#		stat_compare_means(comparisons = my_comparisons1, method="wilcox.test", vjust=0) +
		rremove('legend') +
		rotate_x_text(45) +
		theme(plot.margin = margin(0.3, 0.1, 0.1, 0.1, "cm")) +
		theme(
				title = element_text(size=10, face='bold', color='black'),
				legend.key.size = unit(0.6,"line"),
				legend.text = element_text(size=8),
				axis.title = element_text(size=8), axis.text = element_text(size=8))

HIF1A.boxplot.p <- ggboxplot(
				extracted.cpm.dt.long[gene=='HIF1A',list(tumor.sample.id, ct, log2.ct, exprClass)],
				x = "exprClass", y = "log2.ct", 
#				outlier.shape = NA,
#				title = 'B',
#				color = 'black',
#				fill = "exprClass",
				add = "jitter",
				add.params = list(shape = 16, size = jitter.size, color = "exprClass"),
				ylab = "HIF1A log2(cpm)", xlab = "RCC class") +
		scale_color_manual(values = alpha(exprClass.colors, curr.alpha)) +
#		stat_compare_means(comparisons = my_comparisons1, method="t.test", vjust=0.4)+
		rremove('legend') +
		rotate_x_text(45) +
		theme(plot.margin = margin(0.3, 0.1, 0.1, 0.1, "cm")) +
		theme(
				title = element_text(size=10, face='bold', color='black'),
				legend.key.size = unit(0.6,"line"),
				legend.text = element_text(size=8),
				axis.title = element_text(size=8), axis.text = element_text(size=8))


EPAS1.boxplot.p <- ggboxplot(
				extracted.cpm.dt.long[gene=='EPAS1',list(tumor.sample.id, ct, log2.ct, exprClass)],
				x = "exprClass", y = "log2.ct", 
#				outlier.shape = NA,
#				title = 'B',
#				color = 'black',
#				fill = "exprClass",
				add = "jitter",
				add.params = list(shape = 16, size = jitter.size, color = "exprClass"),
				ylab = "EPAS1 log2(cpm)", xlab = "RCC class") +
		scale_color_manual(values = alpha(exprClass.colors, curr.alpha)) +
#		stat_compare_means(comparisons = my_comparisons1, method="t.test", vjust=0.4)+
		rremove('legend') +
		rotate_x_text(45) +
		theme(plot.margin = margin(0.3, 0.1, 0.1, 0.1, "cm")) +
		theme(
				title = element_text(size=10, face='bold', color='black'),
				legend.key.size = unit(0.6,"line"),
				legend.text = element_text(size=8),
				axis.title = element_text(size=8), axis.text = element_text(size=8))

VHL.HIF1A_AS3.scatter.p <- ggscatter(
				data = extracted.cpm.dt.by.gene[order(-exprClass),],
				x = 'log2.ct_HIF1A-AS3',
				y  = "log2.ct_VHL",
				shape = 16,
				color = "Histology",
				add = "reg.line", 
				add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
				xlab = "HIF1A-AS3 log2(cpm)",
				ylab = "VHL log2(cpm)") +
		rremove("legend") +
		scale_color_manual(values = alpha(exprClass.colors, curr.alpha)) +
		theme(title = element_text(size=10, face='bold', color='black'),
				axis.title = element_text(size=8), axis.text = element_text(size=8))

curr.x.range <- ggplot_build(VHL.HIF1A_AS3.scatter.p)$panel$ranges[[1]]$x.range
curr.y.range <- ggplot_build(VHL.HIF1A_AS3.scatter.p)$panel$ranges[[1]]$y.range

		VHL.HIF1A_AS3.scatter.p +
		stat_cor(method = "pearson",
			label.x = 0.1*(curr.x.range[[2]]-curr.x.range[[1]]) + curr.x.range[[1]],
			label.y = 0.9*(curr.y.range[[2]]-curr.y.range[[1]]) + curr.y.range[[1]],
			size=2.5)

HIF1A.EPAS1.scatter.p <- ggscatter(
				data = extracted.cpm.dt.by.gene[order(-exprClass),],
				x = 'log2.ct_HIF1A',
				y  = "log2.ct_EPAS1",
				shape = 16,
				color = "Histology",
				add = "reg.line", 
				add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
				xlab = "HIF1A log2(cpm)",
				ylab = "EPAS1 log2(cpm)") +
		rremove("legend") +
		scale_color_manual(values = alpha(exprClass.colors, curr.alpha)) +
		theme(title = element_text(size=10, face='bold', color='black'),
				axis.title = element_text(size=8), axis.text = element_text(size=8))

		curr.x.range <- ggplot_build(HIF1A.EPAS1.scatter.p)$panel$ranges[[1]]$x.range
		curr.y.range <- ggplot_build(HIF1A.EPAS1.scatter.p)$panel$ranges[[1]]$y.range
		
		HIF1A.EPAS1.scatter.p +
				stat_cor(method = "pearson",
						label.x = 0.1*(curr.x.range[[2]]-curr.x.range[[1]]) + curr.x.range[[1]],
						label.y = 0.9*(curr.y.range[[2]]-curr.y.range[[1]]) + curr.y.range[[1]],
						size=2.5)

EPAS1.HIF1A_AS3.scatter.p <- ggscatter(
				data = extracted.cpm.dt.by.gene[order(-exprClass),],
				x = 'log2.ct_HIF1A-AS3',
				y  = "log2.ct_EPAS1",
				shape = 16,
				color = "Histology",
				add = "reg.line", 
				add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
				xlab = "HIF1A-AS3 log2(cpm)",
				ylab = "EPAS1 log2(cpm)") +
		rremove("legend") +
		scale_color_manual(values = alpha(exprClass.colors, curr.alpha)) +
		theme(title = element_text(size=10, face='bold', color='black'),
				axis.title = element_text(size=8), axis.text = element_text(size=8))
		
		curr.x.range <- ggplot_build(EPAS1.HIF1A_AS3.scatter.p)$panel$ranges[[1]]$x.range
		curr.y.range <- ggplot_build(EPAS1.HIF1A_AS3.scatter.p)$panel$ranges[[1]]$y.range
		
		EPAS1.HIF1A_AS3.scatter.p +
				stat_cor(method = "pearson",
						label.x = 0.1*(curr.x.range[[2]]-curr.x.range[[1]]) + curr.x.range[[1]],
						label.y = 0.9*(curr.y.range[[2]]-curr.y.range[[1]]) + curr.y.range[[1]],
						size=2.5)

HIF1A.HIF1A_AS3.scatter.p <- ggscatter(
				data = extracted.cpm.dt.by.gene[order(-exprClass),],
				x = 'log2.ct_HIF1A-AS3',
				y  = "log2.ct_HIF1A",
				shape = 16,
				color = "Histology",
				add = "reg.line", 
				add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
				xlab = "HIF1A-AS3 log2(cpm)",
				ylab = "HIF1A log2(cpm)") +
		rremove("legend") +
		scale_color_manual(values = alpha(exprClass.colors, curr.alpha)) +
		theme(title = element_text(size=10, face='bold', color='black'),
				axis.title = element_text(size=8), axis.text = element_text(size=8))

		curr.x.range <- ggplot_build(HIF1A.HIF1A_AS3.scatter.p)$panel$ranges[[1]]$x.range
		curr.y.range <- ggplot_build(HIF1A.HIF1A_AS3.scatter.p)$panel$ranges[[1]]$y.range
		
		HIF1A.HIF1A_AS3.scatter.p +
				stat_cor(method = "pearson",
						label.x = 0.1*(curr.x.range[[2]]-curr.x.range[[1]]) + curr.x.range[[1]],
						label.y = 0.9*(curr.y.range[[2]]-curr.y.range[[1]]) + curr.y.range[[1]],
						size=2.5)

ccRCC.HIF.fig.titles.A.to.H <- ggarrange(
		ggarrange(
				HIF1A_AS3.boxplot.p + ggtitle('A'),
				HIF1A.boxplot.p + ggtitle('B'),
				EPAS1.boxplot.p + ggtitle('C'),
				VHL.boxplot.p + ggtitle('D'),
				ncol=4, nrow=1),
		ggarrange(
				HIF1A.HIF1A_AS3.scatter.p + ggtitle('E'),
				HIF1A.EPAS1.scatter.p + ggtitle('F'),
				EPAS1.HIF1A_AS3.scatter.p + ggtitle('G'),
				VHL.HIF1A_AS3.scatter.p + ggtitle('H'),
				ncol = 4, nrow=1),
		nrow = 2, heights = c(1, 0.8))

ggexport(ccRCC.HIF.fig.titles.A.to.H, width=7, height=5,
		filename=file.path(run.dir, "result/DGE/ccRCC_non_ccRCC/ccRCC.nonccRCC_DGE.log2ct.figure.pdf"))

ccRCC.HIF.fig.no.titles <- ggarrange(
		ggarrange(
				HIF1A_AS3.boxplot.p,
				HIF1A.boxplot.p,
				EPAS1.boxplot.p,
				VHL.boxplot.p ,
				ncol=4, nrow=1),
		ggarrange(
				HIF1A.HIF1A_AS3.scatter.p,
				HIF1A.EPAS1.scatter.p,
				EPAS1.HIF1A_AS3.scatter.p,
				VHL.HIF1A_AS3.scatter.p,
				ncol = 4, nrow=1),
		nrow = 2, heights = c(1, 0.8))

ggexport(ccRCC.HIF.fig.no.titles, width=7, height=5,
		filename=file.path(run.dir, "result/DGE/ccRCC_non_ccRCC/ccRCC.nonccRCC_DGE.log2ct.no.titles.figure.pdf"))
