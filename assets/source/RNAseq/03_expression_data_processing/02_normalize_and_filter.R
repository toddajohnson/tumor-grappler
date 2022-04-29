#!/home/tjohnson/.local/bin/Rscript

options(width=350)

working_dir <- getwd()
gpl_prefix <- basename(working_dir)

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "~/workspace/runs", study.dir )

library(data.table)
library(parallel)

source( file.path(run.dir, "config_files/common_config.R") )
source( file.path(RNAseq_cg_pipeline_dir, "03_expression_data_processing/expression_processing_functions.R") )

print( paste("Attempting to load imported gene expression data for ", study.name, sep="") )
# process.tximport.data.ls has slots for
# sample.info, sample.ids.ls, imported.tx.expr.data, imported.gene.expr.data
load(file = file.path(run.dir, 'result/processed_expression_data/imported.tx.gene.expression.Rdata'))

normalized.expression.data.ls <- normalize_and_filter_counts(
	input.processed.tximport.data.ls = processed.tximport.data.ls,
	input.grouping.column = 'CancerType.short',
	input.samples.to.remove = c('IWK006_T', 'IWK039_T', 'IWK049_T') )

#normalized.expression.data.ls
#sample.ids.filtered, expr_filtered_DGEList, mat.all.samples.natural.cpm, mat.all.samples
save( list=c('normalized.expression.data.ls'),
		file = file.path(run.dir, 'result/processed_expression_data/gene_expression.edgeR.normalized.Rdata'))

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
