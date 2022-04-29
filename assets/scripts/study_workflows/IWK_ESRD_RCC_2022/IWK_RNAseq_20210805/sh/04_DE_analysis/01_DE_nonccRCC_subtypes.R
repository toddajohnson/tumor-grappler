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

sample.info <- sample.info.with.clinical.info.dt[cls8%in%c('B', 'H'),]
clusters.B.H.sample.ct <- nrow(sample.info)
sample.info[,exprClass:=ifelse(cls8=='H', 'non-ccRCC_H', ifelse(cls8=='B', 'non-ccRCC_B', 'UNK'))]
samples.to.keep <- sample.info$tumor.sample.id

## > nrow(cts)
## [1] 39357
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
groups.f <- factor(x = groups, levels=c('non-ccRCC_H', 'non-ccRCC_B'), ordered = TRUE)

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
topTags(lrt)
## > topTags(lrt)
## Coefficient:  y$samples$group.L 
## logFC   logCPM        LR       PValue          FDR
#SLC22A7    8.839865 7.326197 275.53181 7.068385e-62 1.200353e-57
#FABP1      6.692916 5.723530 124.45575 6.695568e-29 5.685207e-25
#XPNPEP2    5.934763 6.681573 121.55172 2.893592e-28 1.637966e-24
#CT69      -5.708411 4.008859 104.75319 1.383373e-24 5.873111e-21
#LINC00671  3.889205 4.757531  81.33160 1.908503e-19 6.482040e-16
#GLTPD2     3.842572 3.605725  78.71018 7.192403e-19 2.035690e-15
#ASPDH      4.171311 3.916485  74.84633 5.088157e-18 1.151846e-14
#GPR143    -2.938970 7.071047  74.68168 5.530705e-18 1.151846e-14
#SLC44A4    4.266439 5.790508  74.48681 6.104471e-18 1.151846e-14
#SLC5A12    6.065283 7.303365  71.80670 2.373454e-17 3.835108e-14

o <- order(lrt$table$PValue)
cpm(y)[o[1:10],]
#			IWK001_T     IWK008_T    IWK013_T    IWK019_T    IWK026_T  IWK028_T
#SLC22A7     0.07205343   0.09894163   0.0000000 287.3531367 396.6067769 244.77352
#FABP1       0.07205343   0.00000000   0.2986491  51.1649304 148.9474441 110.70874
#XPNPEP2     0.57642747   0.49470816   0.6968480 123.0680213 103.5986159 377.19736
#CT69       20.24701501  33.34333017  46.7883631   0.1296135   0.1954691   0.00000
#LINC00671   0.36026717   1.92936183   1.0286803  47.1793152  47.7270353  58.58155
#GLTPD2      0.50437404   0.54417898   0.5309318  16.0396710  31.8940393  21.16812
#ASPDH       0.64848091   0.49470816   0.3318324  17.8218567  44.9253118  23.73893
#GPR143    308.60485862 187.09862714 281.2942938  10.4662904  12.5425997  20.62114
#SLC44A4     3.09829767   0.89047469   0.9623139  74.3333441 131.5506953 112.62317
#SLC5A12     1.58517555   0.59364980   0.3318324  82.6286083 448.5689758 407.55475

summary(decideTests(lrt))
y$samples$group.L
##Down                 310
##NotSig             16224
##Up                   448


library(EnhancedVolcano)

# Import edgeR data
## Coefficient:  y$samples$group.L 
## logFC     logCPM        LR       PValue          FDR

volcanodata <- as.data.table(topTags(lrt, n=nrow(lrt$table)), keep.rownames=TRUE)
setnames(volcanodata,
	old=paste('table.', c('rn', 'logFC', 'logCPM', 'LR', 'PValue', 'FDR'), sep=""),
	new = c("Gene", "logFC",  "logCPM", "LR", "pvalue", "FDR"))

pdf(file=file.path(run.dir, 'result/DGE/non-ccRCC_B_H/non-ccRCC_B_H_Volcano.pdf'))
EnhancedVolcano(
		volcanodata, 
		x = 'logFC',
		y = 'pvalue',
		lab = volcanodata$Gene,
#		AdjustedCutoff = 0.05, 
#		LabellingCutoff = 0.05, 
#		FCCutoff = 2.0,
		title = "non-ccRCC Clusters B vs. H edgeR")
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

ggexport( curr.p, width=2.2, height=2.2, filename=file.path(run.dir, 'result/DGE/non-ccRCC_B_H/non-ccRCC_B_H_Volcano.for.figure.2.2in.pdf'))
ggexport( curr.p, width=2.0, height=2.0, filename=file.path(run.dir, 'result/DGE/non-ccRCC_B_H/non-ccRCC_B_H_Volcano.for.figure.2.0in.pdf'))


pdf(file=file.path(run.dir, 'result/DGE/non-ccRCC_B_H/non-ccRCC_B_H_Volcano.no_title.pdf'))
EnhancedVolcano(
		volcanodata, 
		x = 'logFC',
		y = 'pvalue',
		lab = volcanodata$Gene,
#		AdjustedCutoff = 0.05, 
#		LabellingCutoff = 0.05, 
#		FCCutoff = 2.0,
		title = NULL)
#		col = c("grey30", "forestgreen", "royalblue", "red2"))
dev.off()

## pdf(file=file.path(run.dir, 'result/DGE/non-ccRCC_B_H/non-ccRCC_B_H_Volcano_byFDR.pdf'))
## EnhancedVolcano(
##         volcanodata, 
##         x = 'logFC',
##         y = 'FDR',
##         lab = volcanodata$Gene,
## #		AdjustedCutoff = 0.05, 
## #		LabellingCutoff = 0.05, 
##         FCCutoff = 2.0,
##         title = "ccRCC vs. non-ccRCC edgeR")
## #		col = c("grey30", "forestgreen", "royalblue", "red2"))
## dev.off()


save( list=c('volcanodata', 'sample.info', 'lrt', 'y'),
		file=file.path(run.dir, 'result/DGE/non-ccRCC_B_H/non-ccRCC_B_H_DE.Rdata'))

fwrite(
		volcanodata[,list(Gene, logFC, logCPM, LR, pvalue, FDR)],
		file = file.path(run.dir, 'result/DGE/non-ccRCC_B_H/non-ccRCC_B_H_DGE_list.tsv'),
		sep = '\t')

library(AnnotationHub)
library(ensembldb)
ah <- AnnotationHub()

## query AnnotationHub for available Ensembl gtf files for Ensembl release 102
query(ah, c("Homo sapiens", "release-104"))

## get the version 104 gtf:
gtf <- ah[["AH92109"]]

## generate the annotation database
DbFile <- ensDbFromGRanges(gtf, organism="Homo_sapiens", version=104, genomeVersion="GRCh38")

## we can either generate a database package using the makeEnsembldbPackage
## , or directly load the data
Edb <- EnsDb(DbFile)

## you can then use e.g. genes to get all annotations from all genes
genes.obj <- genes(Edb)

genes.dt <- as.data.table(genes.obj)

volcanodata <- merge(
	x = genes.dt[,list(gene_id, Gene=gene_name, symbol, chr=seqnames, start, end, width, strand, gene_biotype)],
	y = volcanodata,
	by = c('Gene'))


volcanodata[,target:=ifelse(logFC > 0.5 & FDR<0.25, 'Cluster B', ifelse(logFC < -0.5 & FDR<0.25, 'Cluster H', 'Not DE'))]
volcanodata[,Gene.3:=substring(Gene, 1, 3)]
volcanodata[,L.ribosomal:=0L]
volcanodata[,S.ribosomal:=0L]
volcanodata[Gene.3=='RPL',L.ribosomal:=1L]
volcanodata[Gene.3=='RPS',S.ribosomal:=1L]

cls.B.H.Ribosomal.summary <- volcanodata[,list(gene.ct=.N, L.ct=sum(L.ribosomal), S.ct=sum(S.ribosomal)),
	by=list(target)][target!='Not DE',]
cls.B.H.Ribosomal.summary[,ribosomal.genes:=L.ct+S.ct]
cls.B.H.Ribosomal.summary[,ribosomal.gene.propn:=ribosomal.genes/gene.ct]
#	target gene.ct L.ct S.ct ribosomal.genes ribosomal.gene.propn
#1: Cluster B     952    0    1               1           0.00105042
#2: Cluster H     844   64   18              82           0.09715640

cls.B.H.biotype.summary[,total.ct:=sum(biotype.ct), by=list(target)]
cls.B.H.biotype.summary[,biotype.propn:=biotype.ct/total.ct]

save( list=c('volcanodata'),
		file=file.path(run.dir, 'result/DGE/non-ccRCC_B_H/non-ccRCC_B_H_DE.with.Ensembl.gene.info.Rdata'))

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
#length(gct.dt$NAME)
#[1] 20229
#> nrow(ensembl.genes.in.file)
#[1] 20229

ensembl.genes.in.file.no.loc <- ensembl.genes.in.file[is.na(Chromosome)]
ensembl.genes.in.file.loc <- ensembl.genes.in.file[!is.na(Chromosome)]

mat.all.samples.natural.cpm.dt.loc <- mat.all.samples.natural.cpm.dt[ensembl.genes.in.file.loc$GeneName,]

gene.ct <- nrow(ensembl.genes.in.file.loc)

sample.info.with.clinical.info.dt.clusters.B.H <- sample.info[cls8%in%c('B', 'H'),]
clusters.B.H.cls <- file.path(isofox.dir, "result/GSEA/nonccRCC_clusters_B_H/clusters.B.H.cls")

clusters.B.H.sample.ct <- nrow(sample.info.with.clinical.info.dt.clusters.B.H)

writeLines( c(
				paste(c(clusters.B.H.sample.ct, '2', '1'), collapse=' '),
				paste( c("#", unique(sample.info.with.clinical.info.dt.clusters.B.H$exprClass)), collapse=' '),
				paste( c(sample.info.with.clinical.info.dt.clusters.B.H$exprClass), collapse=' ')), clusters.B.H.cls)


clusters.B.H.gct.dt <- mat.all.samples.natural.cpm.dt.loc[,c('NAME', 'Description', sample.info.with.clinical.info.dt.clusters.B.H$tumor.sample.id), with=FALSE]
clusters.B.H.gct.file <- file.path(isofox.dir, "result/GSEA/nonccRCC_clusters_B_H/isofox_gene_expression.clusters.B.H.gct")

writeLines(c('#1.2', paste(c(gene.ct, clusters.B.H.sample.ct), collapse='\t')), clusters.B.H.gct.file)
fwrite(clusters.B.H.gct.dt, clusters.B.H.gct.file, sep='\t', col.names=TRUE, append=TRUE)
