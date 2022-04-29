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

sample.info <- sample.info.with.clinical.info.dt[Histology%in%c('ccRCC', 'ACD-RCC'),]
sample.info[,exprClass:=Histology]
samples.to.keep <- sample.info$tumor.sample.id

## > nrow(cts)
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
groups.f <- factor(x = groups, levels=c('ccRCC', 'ACD-RCC'), ordered = TRUE)

y <- DGEList(cts, group=groups.f)

y <- scaleOffset(y, normMat)
# filtering

keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]

pdf(file=file.path(run.dir, 'result/DGE/ccRCC_ACD.RCC/MDS.pdf'))
plotMDS(y)
dev.off()

#
# TMM normalization and log transformation
#
y <- estimateDisp(y)

pdf(file=file.path(run.dir, 'result/DGE/ccRCC_ACD.RCC/BCV.pdf'))
plotBCV(y)
dev.off()

fit <- glmFit(y)
lrt <- glmLRT(fit)
topTags(lrt)
## > topTags(lrt)
#Coefficient:  y$samples$group.L 
#logFC    logCPM       LR       PValue          FDR
#CA9      -6.683282  8.196256 92.45509 6.887018e-22 1.243176e-17
#ADM      -3.822228  9.552496 74.18028 7.129893e-18 6.435085e-14
#DNAH11   -7.193443  4.981269 72.27175 1.875129e-17 1.128265e-13
#PFKFB2    2.615395  4.766955 69.14521 9.147414e-17 4.127999e-13
#FAM102A   1.407674  5.806817 61.83308 3.738415e-15 1.349643e-11
#MIR210HG -2.996257  5.118196 61.42750 4.593508e-15 1.381957e-11
#PPFIA4   -3.859971  6.811417 59.97122 9.625468e-15 2.313954e-11
#HILPDA   -4.008712  8.421163 59.84650 1.025518e-14 2.313954e-11
#FUT11    -1.671065  5.331781 58.68736 1.848249e-14 3.706971e-11
#ANGPTL4  -4.143344 10.124125 57.10928 4.122558e-14 7.441630e-11

summary(decideTests(lrt))
#y$samples$group.L
#Down                 941
#NotSig             15777
#Up                  1333


pdf(file=file.path(run.dir, 'result/DGE/ccRCC_ACD.RCC/MD.pdf'))
plotMD(lrt)
abline(h=c(-1, 1), col="blue")
dev.off()

library(EnhancedVolcano)

volcanodata <- as.data.table(topTags(lrt, n=nrow(lrt$table)), keep.rownames=TRUE)
setnames(volcanodata,
	old=paste('table.', c('rn', 'logFC', 'logCPM', 'LR', 'PValue', 'FDR'), sep=""),
	new = c("Gene", "logFC",  "logCPM", "LR", "pvalue", "FDR"))

pdf(file=file.path(run.dir, 'result/DGE/ccRCC_ACD.RCC/ccRCC_ACD.RCC_Volcano.pdf'))
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

save( list=c('volcanodata', 'sample.info', 'lrt', 'y'),
		file=file.path(run.dir, 'result/DGE/ccRCC_ACD.RCC/ccRCC_ACD.RCC_DE.Rdata'))

fwrite(
		volcanodata[,list(Gene, logFC, logCPM, LR, pvalue, FDR)],
		file = file.path(run.dir, 'result/DGE/ccRCC_ACD.RCC/ccRCC_ACD.RCC_DGE_list.tsv'),
		sep = '\t')


extracted.cpm.dt <- as.data.table(extracted.cpm, keep.rownames=TRUE)
setnames(extracted.cpm.dt, old='rn', new='gene')

extracted.cpm.dt.long <- melt(
		data = extracted.cpm.dt,
		id.vars = c('gene'),
		measure.vars = colnames(extracted.cpm.dt[,-c('gene')]),
		variable.name = 'tumor.sample.id',
		value.name = 'ct')

extracted.cpm.dt.long[,log2.ct:=log2(ct)]

extracted.cpm.dt.long <- merge(
		x = data.table(
				tumor.sample.id = rownames(y$samples),
				exprClass = y$samples$group),
		y = extracted.cpm.dt.long,
		by = c('tumor.sample.id'))

extracted.cpm.dt.by.gene <- dcast(
		extracted.cpm.dt.long,
		formula = exprClass + tumor.sample.id  ~ gene,
		value.var = c('ct', 'log2.ct'))

extracted.cpm.dt.by.gene[,Histology:=exprClass]

makeTransparent<-function(someColor, alpha=100)
{
	newColor<-col2rgb(someColor)
	apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
						blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

t.test(x=extracted.cpm.dt.by.gene[exprClass=='non-ccRCC',]$log2.ct_EPAS1, y=extracted.cpm.dt.by.gene[exprClass=='ccRCC',]$log2.ct_EPAS1)

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
	filename=file.path(run.dir, "result/isofox/result/DE/ccRCC_non_ccRCC/ccRCC.nonccRCC_DGE.figure.pdf"))

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
filename=file.path(run.dir, "result/isofox/result/DE/ccRCC_non_ccRCC/ccRCC.nonccRCC_DGE.no.titles.figure.pdf"))


#ggexport(HIF1A_AS3.boxplot.p + ggtitle('A'), width=2.5, height=3, filename=file.path(run.dir, "result/isofox/result/DE/ccRCC_non_ccRCC/ccRCC_DE.boxplot.figure.pdf"))

#ggexport(ccRCC.HIF.fig, width=4.5, height=5, filename=file.path(run.dir, "result/isofox/result/DE/ccRCC_non_ccRCC/ccRCC_DE.figure.pdf"))

HIF1A.dt <- data.table(
	names)
pdf(file=file.path(run.dir, 'result/isofox/result/DE/ccRCC_non_ccRCC/ccRCC_nonccRCC_HIF1A_vsAS3.pdf'))
ggscatterplot()
mean(cpm(y)['HIF1A',][which(cpm(y)['HIF1A-AS3',]==0)])
#[1] 124.9022
mean(cpm(y)['HIF1A',][which(cpm(y)['HIF1A-AS3',]>0)])
#[1] 62.91495


# For ordering sashimi plot
HIF1A_AS3 <- cpm(y)[c('HIF1A-AS3'),]
HIF1A_AS3 <- HIF1A_AS3[order(HIF1A_AS3)]
HIF1A_AS3
IWK001_T     IWK008_T     IWK013_T     IWK016_T     IWK018_T     IWK019_T     IWK026_T     IWK028_T     IWK020_T     IWK017_T     IWK024_T     IWK029_T     IWK002_T     IWK015_T     IWK021_T     IWK012_T     IWK009_T     IWK022_T     IWK010_T     IWK025_T     IWK005_T     IWK048_T 
0.00000000   0.00000000   0.00000000   0.00000000   0.00000000   0.00000000   0.00000000   0.05551224  11.46560233  15.83670133  38.74375676  49.51869832  50.20472527  55.64860962  64.63038692  67.02070790  70.34304684  86.39830920  89.70442026  93.82434979  98.80286559 176.98471861 
> 
