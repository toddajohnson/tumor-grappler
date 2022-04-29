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
library(openxlsx)

source( file.path(run.dir, "config_files/common_config.R") )

#wgs.run.dir <- '/Users/toddjohnson/OneDrive/Documents/RIKEN/Projects/IWK_ESRD_RCC/IWK_WGS_HMF_20210726'
#run.dir <- '/Users/toddjohnson/OneDrive/Documents/RIKEN/Projects/IWK_ESRD_RCC/IWK_RNAseq_20210805'

#ref.dir <- "/home/tjohnson/reference/HMF/38/dbs/ensembl_data_cache"

#HOME.dir <- "~/HGC_mounts/HGC/"
HOME.dir <- "/home/tjohnson"

isofox.dir <- file.path(run.dir, 'result', 'isofox')

load(file = paste(isofox.dir, "/isofox.gene_expression.edgeR.normalized.Rdata", sep=""))

#load(file=file.path(isofox.dir, 'imported.tx.gene.expression.Rdata'))
load(file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))

load(file = paste(wgs.run.dir, "/result/sample_summaries/clinical_data_with_colors.Rdata", sep=""))

library(org.Hs.eg.db)

e2s = toTable(org.Hs.egSYMBOL)
e2s <- as.data.table(e2s)
setkey(e2s, 'gene_id')

kidney.gene.info.df <- read.xlsx(file.path(run.dir, '/result/TCGA/Lindgren_2017/Lindgren et al 2017_TableS4.xlsx'))
kidney.gene.info.dt <- as.data.table(kidney.gene.info.df)
kidney.gene.info.dt[,ENTREZ_ID:=as.character(ENTREZ_ID)]
kidney.gene.info.dt.filt <- kidney.gene.info.dt[which(ENTREZ_ID%in%e2s$gene_id),]

nrow(kidney.gene.info.dt)
##1171
nrow(kidney.gene.info.dt.filt)
## 1167

e2s.ls <- e2s$symbol
names(e2s.ls) <- e2s$gene_id

kidney.gene.info.dt.filt[,gene:=e2s.ls[ENTREZ_ID]]

kidney.genes.high.FC <- kidney.gene.info.dt.filt[Included.Figure.3=='yes',]

high.FC.kidney.genes.ct <- nrow(kidney.genes.high.FC)

high.FC.kidney.genes.in.ESRD.filt <- kidney.genes.high.FC[which(gene%in%rownames(mat.all.samples.filt)),]$gene
high.FC.kidney.genes.in.ESRD.filt.ct <- length(high.FC.kidney.genes.in.ESRD.filt)

print(paste('There were ', high.FC.kidney.genes.in.ESRD.filt.ct, '/', high.FC.kidney.genes.ct, ' in the filtered ESRD RNA-seq', sep=''))
##print(paste('There were ', high.FC.kidney.genes.in.ESRD.filt, '/', high.FC.kidney.genes.ct, sep=''))
##[1] "There were 427/708"

high.FC.kidney.genes.in.ESRD.unfilt <- kidney.genes.high.FC[which(gene%in%rownames(mat.all.samples)),]$gene
high.FC.kidney.genes.in.ESRD.unfilt.ct <- length(high.FC.kidney.genes.in.ESRD.unfilt)

kidney.genes.high.FC[which(!gene%in%rownames(mat.all.samples)),]

print(paste('There were ', high.FC.kidney.genes.in.ESRD.unfilt, '/', high.FC.kidney.genes.ct, ' in the un-filtered ESRD RNA-seq', sep=''))
##"There were 688/705 in the un-filtered ESRD RNA-seq"


#mat.all.samples.filt <- mat.all.samples[sds.all.samples > 1.5, ]
mat.all.samples.filt <- mat.all.samples[which(rownames(mat.all.samples)%in%high.FC.kidney.genes.in.ESRD.unfilt), ]

install.packages("factoextra")
library(factoextra)

kidney.genes.pca <- prcomp(t(mat.all.samples.filt))

pdf(onefile=TRUE, file=file.path(run.dir, 'result/PCA/kidney_genes_ESRD_RCC_PCA.pdf'))

fviz_eig(kidney.genes.pca)

fviz_pca_ind(kidney.genes.pca,
		col.ind = "cos2", # Color by the quality of representation
		gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
		repel = TRUE     # Avoid text overlapping
)

fviz_pca_var(kidney.genes.pca,
		col.var = "contrib", # Color by contributions to the PC
		gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
		repel = TRUE     # Avoid text overlapping
)

dev.off()

mat.filt <- mat.all.samples[,colnames(mat.all.samples)[which(!colnames(mat.all.samples)%in%c('IWK006_T'))]]

sds.filt <- apply(mat.filt , 1, sd)
mat.filt <- mat.filt[sds.filt > 1.5, ]

# scale per gene
mat.all.samples.filt <- mat.all.samples.filt %>% t %>% scale %>% t
mat.all.samples <- mat.all.samples %>% t %>% scale %>% t
mat.filt <- mat.filt %>% t %>% scale %>% t


mat.all.samples.filt

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

sample.info <- sample.info.with.clinical.info.dt
sample.info[,exprClass:=cls3]

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

comparisons.dt <- data.table(
	cls.col = c('cls3', 'cls8'),
	cls.1 = c('B', 'B',),
	cls.2 = c('C', 'H', ))
		
		
# Creating a DGEList object for use in edgeR.
lrt.fits.ls <- lapply(
	X = comparisons.ls,
	FUN = function( curr.comparision ){
		groups <- sample.info$exprClass
		names(groups) <- sample.info$tumor.sample.id
		groups <- groups[colnames(cts)]
		groups.f <- factor(x = groups, levels=c('B', 'C', 'A'), ordered = TRUE)
		
		y <- DGEList(cts, group=groups.f)
		
		y <- scaleOffset(y, normMat)
# filtering
		
		keep <- filterByExpr(y)
		y <- y[keep,,keep.lib.sizes=FALSE]
		
#
# TMM normalization and log transformation
#
		y <- estimateDisp(y)
		
		fit <- glmFit(y)
		lrt <- glmLRT(fit)
		
	})
topTags(lrt)
## > topTags(lrt)

## Coefficient:  y$samples$group.L 
##         logFC     logCPM        LR       PValue          FDR
## HIF1A-AS3  8.956940  5.4730263 260.07831 1.649854e-58 2.793863e-54
## VEGFA      3.528596 11.6556107 149.71792 1.998092e-34 1.691785e-30
## SYT14     -6.238924  0.9644833 148.52748 3.637707e-34 2.053364e-30
## CCND1      2.226365  9.0729975 139.62632 3.213139e-32 1.360282e-28
## PCSK6      3.018219  5.9993382 124.68659 5.960247e-29 2.018616e-25
## NDUFA4L2   4.770442 10.1637754 113.55987 1.626939e-26 4.591765e-23
## EHD2       1.821529  7.7436289 106.56780 5.536558e-25 1.339373e-21
## CP         5.644333  7.8879619 104.40242 1.651291e-24 3.495371e-21
## SLC47A2   -5.287568  7.4110902  97.94118 4.309968e-23 8.109445e-20
## STC2       2.605020  6.4266925  92.67768 6.154377e-22 9.460382e-19


summary(decideTests(lrt))
## y$samples$group.L
## Down                2183
## NotSig             12782
## Up                  1969

check <- data.table(sample=rownames(y$samples),
	HIF1A=cpm(y)['HIF1A',], EPAS1=cpm(y)['EPAS1',], HIF1A_AS3=cpm(y)['HIF1A-AS3',], SLC6A3=cpm(y)['SLC6A3',], SLC47A2=cpm(y)['SLC47A2',], SLC22A7=cpm(y)['SLC22A7',])

check <- merge(
	x=sample.info[,list(sample=tumor.sample.id, Histology, years.dialysis, cls3, cls8)],
	y=check)

check[,list(mean.years.dialysis=mean(years.dialysis)), by=list(cls3)]

#> check[,list(mean.years.dialysis=mean(years.dialysis)), by=list(cls3, cls8)][order(mean.years.dialysis),]
#	cls3 cls8 mean.years.dialysis
#1:    A    A           0.4916667
#2:    B    F           4.2930556
#3:    B    G           4.3548611
#4:    B    C           7.8064815
#5:    C    D          11.4791667
#6:    B    E          12.8277778
#7:    C    B          14.5259259
#8:    C    H          22.7935185

#> check[,list(mean.years.dialysis=mean(years.dialysis)), by=list(cls3)][order(mean.years.dialysis),]
#	cls3 mean.years.dialysis
#1:    A           0.4916667
#2:    B           7.0357143
#3:    C          16.8645833

#> check[order(cls3, cls8),]
#	sample Histology years.dialysis cls3 cls8     HIF1A      EPAS1   HIF1A_AS3       SLC6A3      SLC47A2      SLC22A7
#1: IWK040_T     chRCC    0.491666667    A    A  81.23395  239.67461   0.3265322 0.000000e+00    0.0000000   0.00000000
#2: IWK017_T     ccRCC    3.150000000    B    C  39.24830  642.76631  15.8287080 6.124459e+00    6.1244592   1.85458977
#3: IWK020_T     ccRCC    7.672222222    B    C  37.79696  360.20592  11.2605268 1.457758e+01    1.2220727   1.65852720
#4: IWK022_T     ccRCC    0.294444444    B    C  72.27684  827.28455  86.6816589 5.354334e+02    4.3322778   0.50543241
#5: IWK024_T     ccRCC    9.025000000    B    C  38.57374  663.97055  39.0682800 2.909008e+02    4.5268990   0.49453519
#6: IWK025_T     ccRCC    4.666666667    B    C  82.85696  820.17782  89.2199276 3.563264e+02    0.5533019   0.09221698
#7: IWK029_T     ccRCC   22.030555556    B    C  46.68443  996.30162  50.2077883 1.408240e+02    0.6606288   0.84413679
#8: IWK015_T     ccRCC   25.647222222    B    E  72.69429  996.40290  56.3216991 1.070767e+01    2.0956911   0.09823552
#9: IWK048_T      pRCC    0.008333333    B    E  90.96399  363.93202 174.8546944 2.194620e+02    0.9887390   3.95495611
#10: IWK009_T     ccRCC    3.083333333    B    F 110.40058  833.68743  70.3717189 3.792364e+02    0.3953467   0.09883668
#11: IWK021_T     ccRCC    5.502777778    B    F 185.75407 1029.21776  65.1358625 1.558282e+02    0.2032320   0.05080801
#12: IWK002_T     ccRCC    2.216666667    B    G  37.63852  343.05644  50.6252472 7.521383e+02    9.9986300   2.70077936
#13: IWK005_T     ccRCC    6.619444444    B    G  21.33223  421.54958  98.0425387 5.234967e+02    1.8094300   1.90466321
#14: IWK010_T     ccRCC    0.416666667    B    G  26.38959  475.36310  88.8332603 1.345318e+03    1.0015024  11.36705190
#15: IWK012_T     ccRCC    8.166666667    B    G  21.70034  459.47699  67.2758985 4.200877e+02    3.7697702  11.74428400
#16: IWK019_T      pRCC   12.813888889    C    B  63.19611  295.27490   0.0000000 2.023355e-01  212.4522499 299.01811104
#17: IWK026_T   ACD-RCC   18.447222222    C    B  86.20720  198.50608   0.0000000 6.413221e-01  518.3233119 410.81746740
#18: IWK028_T   ACD-RCC   12.316666667    C    B  56.88166  205.91828   0.0555485 1.110970e-01  835.8382351 248.52397581
#19: IWK016_T    ccpRCC   16.525000000    C    D 286.43150  204.40287   0.0000000 2.229039e-01   27.0669013   1.43295360
#20: IWK018_T   ACD-RCC    6.433333333    C    D  77.85379  230.94227   0.0000000 1.027095e-01   72.1020572   0.10270948
#21: IWK001_T      pRCC   15.833333333    C    H  31.88849   56.73100   0.0000000 7.263892e-02   27.9659851   0.07263892
#22: IWK008_T     ccRCC   28.102777778    C    H 117.93229   77.09250   0.0000000 7.398512e-01  665.9647482   0.09864683
#23: IWK013_T   ACD-RCC   24.444444444    C    H 219.50258   82.74032   0.0000000 7.854013e-01 1515.8927602   0.00000000
#	sample Histology years.dialysis cls3 cls8     HIF1A      EPAS1   HIF1A_AS3       SLC6A3      SLC47A2      SLC22A7

library(EnhancedVolcano)

# Import edgeR data
## Coefficient:  y$samples$group.L 
## logFC     logCPM        LR       PValue          FDR

volcanodata <- as.data.table(topTags(lrt, n=nrow(lrt$table)), keep.rownames=TRUE)
setnames(volcanodata,
	old=paste('table.', c('rn', 'logFC', 'logCPM', 'LR', 'PValue', 'FDR'), sep=""),
	new = c("Gene", "logFC",  "logCPM", "LR", "pvalue", "FDR"))

pdf(file=file.path(run.dir, 'result/isofox/result/DE/ccRCC_non_ccRCC/ccRCC_nonccRCC_Volcano.pdf'))
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

pdf(file=file.path(run.dir, 'result/isofox/result/DE/ccRCC_non_ccRCC/ccRCC_nonccRCC_Volcano_byFDR.pdf'))
EnhancedVolcano(
		volcanodata, 
		x = 'logFC',
		y = 'FDR',
		lab = volcanodata$Gene,
#		AdjustedCutoff = 0.05, 
#		LabellingCutoff = 0.05, 
		FCCutoff = 2.0,
		title = "ccRCC vs. non-ccRCC edgeR")
#		col = c("grey30", "forestgreen", "royalblue", "red2"))
dev.off()


save( list=c('volcanodata', 'sample.info', 'lrt', 'y'),
		file=file.path(run.dir, 'result/isofox/result/DE/ccRCC_non_ccRCC/ccRCC_nonccRCC_DE.Rdata'))

fwrite(
		volcanodata[,list(Gene, logFC, logCPM, LR, pvalue, FDR)],
		file = file.path(run.dir, 'result/isofox/result/DE/ccRCC_non_ccRCC/ccRCC_nonccRCC_DGE_list.tsv'),
		sep = '\t')

extracted.cpm <- cpm(y, log=FALSE)[c('HIF1A', 'HIF1A-AS3', 'VHL', 'EPAS1'),]
extracted.log.cpm <- cpm(y, log=TRUE)[c('HIF1A', 'HIF1A-AS3', 'VHL', 'EPAS1'),]

volcanodata[Gene%in%c('HIF1A', 'HIF1A-AS3', 'VHL', 'EPAS1'),]
## > volcanodata[Gene%in%c('HIF1A', 'HIF1A-AS2', 'VHL', 'EPAS1'),]
#	Gene      logFC   logCPM         LR       pvalue          FDR adjust.method        comparison test
#1: HIF1A-AS3  8.9569400 5.473026 260.078314 1.649854e-58 2.793863e-54            BH y$samples$group.L  glm
#2:     EPAS1  1.3993977 8.913715  35.144971 3.060508e-09 2.115373e-07            BH y$samples$group.L  glm
#3:       VHL -0.5659945 5.034772  11.526025 6.862859e-04 5.468972e-03            BH y$samples$group.L  glm
#4:     HIF1A -0.6207052 6.343994   5.408615 2.003760e-02 7.213367e-02            BH y$samples$group.L  glm

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
		filename=file.path(run.dir, "result/isofox/result/DE/ccRCC_non_ccRCC/ccRCC.nonccRCC_DGE.log2ct.figure.pdf"))

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
		filename=file.path(run.dir, "result/isofox/result/DE/ccRCC_non_ccRCC/ccRCC.nonccRCC_DGE.log2ct.no.titles.figure.pdf"))


###

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
