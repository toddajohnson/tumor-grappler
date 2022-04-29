#!/usr/bin/env Rscript

# module use /usr/local/package/modulefiles/
# module load R/4.0.2
# . ~/.R-4.1.0_setup.sh

options(width=350)
options(width=210)
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
library(magick)
library(ggpubr)

source( file.path(run.dir, "config_files/common_config.R") )

#wgs.run.dir <- '/Users/toddjohnson/OneDrive/Documents/RIKEN/Projects/IWK_ESRD_RCC/IWK_WGS_HMF_20210726'
#run.dir <- '/Users/toddjohnson/OneDrive/Documents/RIKEN/Projects/IWK_ESRD_RCC/IWK_RNAseq_20210805'
#wgs.run.dir <- '/Users/tajohnson/OneDrive/Documents/RIKEN/Projects/IWK_ESRD_RCC/IWK_WGS_HMF_20210726'
#run.dir <- '/Users/tajohnson/OneDrive/Documents/RIKEN/Projects/IWK_ESRD_RCC/IWK_RNAseq_20210805'

#ref.dir <- "/home/tjohnson/reference/HMF/38/dbs/ensembl_data_cache"

#HOME.dir <- "~/HGC_mounts/HGC/"
HOME.dir <- "/home/tjohnson"

isofox.dir <- file.path(run.dir, 'result', 'isofox')

load(file=file.path(isofox.dir, 'imported.tx.gene.expression.Rdata'))
load(file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))

gene.mat <- fread(file.path(isofox.dir, "isofox.gene_expression_matrix.csv"))
tx.mat <- fread(file.path(isofox.dir, "isofox.transcript_expression_matrix.csv"))


check <- tx.mat[GeneName%in%c('HIF1A-AS3', 'HIF1A', 'EPAS1', 'VHL'),]
check.dt <- melt(check, id.vars=c('GeneId', 'GeneName', 'TransName'), value.name='tx.cts', variable.name='tumor.sample.id', value.factor=FALSE, variable.factor=FALSE)
check.dt <- dcast(check.dt, formula= tumor.sample.id + GeneId + GeneName ~ TransName, value.var='tx.cts')

txi <- copy(imported.gene.expr.data)

adjTPM <- txi$abundance
# adjTPM should sum to close to 1 million
#colSums(adjTPM)
##	IWK001_T IWK002_T IWK005_T IWK006_T IWK008_T IWK009_T IWK010_T IWK012_T IWK013_T IWK015_T IWK016_T IWK017_T IWK018_T IWK019_T IWK020_T IWK021_T IWK022_T IWK024_T IWK025_T IWK026_T IWK028_T IWK029_T IWK039_T IWK040_T IWK048_T IWK049_T 
##	1006125  1004020  1008519  1059575  1049016  1003132  1000226  1000241  1026875  1002011  1012465  1000152  1069968  1003368  1001236  1001430  1000358  1009295  1003912  1017387  1029253  1001414  1079054  1009931  1012630  1022590 

cts <- txi$counts
normMat <- txi$length

## > nrow(cts)
# 39357
## [1] 37640
## > ncol(cts)
## [1] 26

#non.zero.cts <- apply(X=cts, MAR=2, FUN=function(x){length(which(x>0))})
#non.zero.cts.mean <- mean(non.zero.cts)
#non.zero.cts.sd <- sd(non.zero.cts)
#non.zero.cts.outlier.cutoff <- non.zero.cts.mean-1.5*non.zero.cts.sd

#samples.to.remove <- names(non.zero.cts[which(non.zero.cts<non.zero.cts.outlier.cutoff )])
## sample with high mitochondrial RNA content, 2 no tumor samples
samples.to.remove <- c('IWK006_T', 'IWK039_T', 'IWK049_T')

cols.to.keep <- colnames(cts)
cols.to.keep <- cols.to.keep[which(!cols.to.keep%in%samples.to.remove)]

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

groups <- sample.info.with.clinical.info.dt$CancerType
names(groups) <- sample.info.with.clinical.info.dt$sample.id
groups <- groups[colnames(cts)]
groups.f <- factor(x = groups, ordered = FALSE)

y <- DGEList(cts, group=groups.f)

y <- scaleOffset(y, normMat)
# filtering

keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]

# TMM normalization and log transformation
mat.all.samples.natural.cpm <- cpm(y, log=FALSE)

mat.all.samples <- cpm(y, log=TRUE)

# remove low-variance genes
sds.all.samples <- apply(mat.all.samples , 1, sd)
mat.all.samples.filt <- mat.all.samples[sds.all.samples > 1.5, ]

mat.filt <- mat.all.samples[,colnames(mat.all.samples)[which(!colnames(mat.all.samples)%in%c('IWK006_T'))]]

sds.filt <- apply(mat.filt , 1, sd)
mat.filt <- mat.filt[sds.filt > 1.5, ]

# scale per gene
mat.all.samples.filt <- mat.all.samples.filt %>% t %>% scale %>% t
mat.all.samples <- mat.all.samples %>% t %>% scale %>% t
mat.filt <- mat.filt %>% t %>% scale %>% t

# clustering of tumors
cluster_columns <- mat.filt %>%
	t %>%
	dist %>%
	hclust(method="ward.D2")

k.ls <- c(3, 6, 8)
cluster.names.ls <- LETTERS[1:max(k.ls)]
k.colors.ls <- get_palette(palette='jco', k=5+max(k.ls))[6:(5+max(k.ls))]
names(k.colors.ls) <- cluster.names.ls
#
# prepare data for cluster annotation A, B, C
#
cls <- cutree(cluster_columns, k=3)  # cut dendrogram at k = 5
n <- names(cls)
cls_map <- c("C", "B", "A")
cls <- cls_map[cls]  # change cluster names from 1-3 to A-C
names(cls) <- n
colors.3.ls <- c("A" = "#4DAF4A", "B" = "#984EA3", "C" = "#FF7F00")

cls8 <- cutree(cluster_columns, k=8)
n8 <- names(cls8)
cls8_map <- LETTERS[8:1]
cls8 <- cls8_map[cls8]
names(cls8) <- n8

sample.info.with.clinical.info.dt[,DialysisPeriod:=cut(years.dialysis, breaks=2, labels=c('0-14 years', '15-29 years'))]

# prepare data for histology annotation
his <- sample.info.with.clinical.info.dt[,list(id=sample.id, CancerType, DialysisPeriod)]
his[,Clear.cell:=ifelse(CancerType=='Clear cell', TRUE, FALSE)]
his[,ACD.RCC:=ifelse(CancerType=='ACD-RCC', TRUE, FALSE)]
his[,Papillary:=ifelse(CancerType=='Papillary', TRUE, FALSE)]
his[,Chromophobe:=ifelse(CancerType=='Chromophobe', TRUE, FALSE)]
his[,Papillary.Clear.papillary:=ifelse(CancerType=='Clear cell papillary', TRUE, FALSE)]

setnames(his,
	old = c('Clear.cell', 'ACD.RCC', 'Papillary.Clear.papillary'),
	new = c('Clear cell', 'ACD-RCC', 'Clear cell papillary'))

his <- as.data.frame(his)
rownames(his) <- his$id
his <- his[names(cls),]

# prepare annotation tracks
bw <- c("TRUE" = "black", "FALSE" = "gray80")
bluered = c('0-14 years' = 'blue', '15-29 years' = 'red')

top_annotation <- HeatmapAnnotation(
	'Cluster (k=3)' = cls,
	'Cluster (k=8)' = cls8,
	'Dialysis period' = his$DialysisPeriod,
	`Clear cell` = his$"Clear cell",
	`ACD-RCC` = his$"ACD-RCC",
	Papillary = his$Papillary,
	`Clear cell papillary` = his$"Clear cell papillary",
	Chromophobe = his$Chromophobe,
	
	col = list(
			'Cluster (k=3)' = colors.3.ls,
			'Cluster (k=8)' = k.colors.ls,
			'Dialysis period' = bluered,
			'Clear cell' = bw,
			 `ACD-RCC` = bw,
			 Papillary = bw,
			 `Clear cell papillary` = bw,
			 Chromophobe = bw),
	 annotation_name_gp = gpar(fontsize = 8),
	 show_legend = FALSE,
	 simple_anno_size = unit(0.25, "cm"))


# draw heatmap
col <- colorRamp2(c(-.5, 0, .5), c("#377EB8", "white", "#E41A1C"))
h <- Heatmap(
	mat.filt,
	name = "Expression",
	heatmap_legend_param = list(title_gp = gpar(fontsize=9), labels_gp = gpar(fontsize=8)),
	col = col,
#	show_row_names=F,
	show_row_names=T,
	km = 4,
	cluster_columns = cluster_columns,
	clustering_method_rows="ward.D2",
	top_annotation = top_annotation,
	column_names_gp = gpar(fontsize = 9),
	height = unit(7.3, "cm")
)


lgd_list = list(
		Legend(labels = names(colors.3.ls), title = "Cluster (k=3)", title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 8), legend_gp = gpar(fill=colors.3.ls)),
		Legend(labels = names(k.colors.ls), title = "Cluster (k=8)", title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 8), , legend_gp = gpar(fill=k.colors.ls)),
		Legend(labels = c("0-14 years", "15-29 years"), title = "Dialysis period", title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 8), , legend_gp = gpar(fill = bluered)),
		Legend(labels = c('TRUE', 'FALSE'), title = "Histology", title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 8), legend_gp = gpar(fill = bw)))


HM <- Heatmap(
		mat.filt,
		km = 2)
HM <- draw(HM)
#pdf( file.path(isofox.dir, "figures", "heatmap_twolevel_cls.pdf"), width=7, height=5)
HM <- draw(
		h,
		annotation_legend_side = "left",
		annotation_legend_list = lgd_list,
		padding = unit(c(2, 2, 2, 8), "mm")
)
r.dend <- row_dend(HM)  #Extract row dendrogram
rcl.list <- row_order(HM)  #Extract clusters (output is a list)
lapply(rcl.list, function(x) length(x))  #check/confirm size clusters
## $`4`
## [1] 541
## 
## $`1`
## [1] 1078
## 
## $`3`
## [1] 626
## 
## $`2`
## [1] 984

for (i in 1:length(row_order(HM))){
	if (i == 1) {
				clu <- t(t(row.names(mat.filt[row_order(HM)[[i]],])))
				out <- cbind(clu, paste("cluster", i, sep=""))
				colnames(out) <- c("GeneID", "Cluster")
				} else {
				clu <- t(t(row.names(mat.filt[row_order(HM)[[i]],])))
				clu <- cbind(clu, paste("cluster", i, sep=""))
				out <- rbind(out, clu)
				}
	}
out  <- as.data.table(out)
fwrite( out[Cluster=='cluster1',], file=file.path(run.dir, 'result/adhoc/selected_genes/chRCC.cluster.genes.tsv'), sep='\t')

ggexport(HM, width=7, height=5, filename=file.path(isofox.dir, "figures", "heatmap_twolevel_cls.20220413.pdf"))

#dev.off()

## pdf( file.path(isofox.dir, "figures", "heatmap_twolevel_cls.20220323.pdf"), width=7, height=5)
## draw(
##         h,
##         annotation_legend_side = "left",
##         annotation_legend_list = lgd_list,
##         padding = unit(c(2, 2, 2, 8), "mm")
## )
## dev.off()

save( list=c('y', 'mat.all.samples.natural.cpm', 'mat.all.samples', 'mat.all.samples.filt', 'mat.filt', 'his', 'cls', 'cls8'),
	file = paste(isofox.dir, "/isofox.gene_expression.edgeR.normalized.Rdata", sep=""))

# export CLS file for GSEA
#
outfile <- file.path(isofox.dir, "result", "heatmap.cls")
cat(length(cls),  length(unique(cls)), "1\n", file=outfile)
cat("#", unique(cls), "\n", file=outfile, append=T)
cat(cls, file=outfile, append=T)
