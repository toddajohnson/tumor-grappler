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
sample.ct <- nrow(sample.info)
sample.info[,exprClass:=ifelse(years.dialysis>14.1, 'Long', ifelse(years.dialysis<14.1, 'Short', 'UNK'))]
samples.to.keep <- sample.info$tumor.sample.id

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
groups.f <- factor(x = groups, levels=c('Short', 'Long'), ordered = TRUE)

y <- DGEList(cts, group=groups.f)

y <- scaleOffset(y, normMat)
# filtering

keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]

pdf(file=file.path(run.dir, 'result/DGE/DialysisPeriod_adj_for_cls3/MDS.pdf'))
plotMDS(y)
dev.off()



tumor_class <- sample.info$cls3
dialysis_period_class <- factor(sample.info$exprClass, levels=c('Short', 'Long'), ordered = TRUE)

design <- model.matrix(~tumor_class+dialysis_period_class)
rownames(design) <- colnames(y)

y <- estimateDisp(y)

pdf(file=file.path(run.dir, 'result/DGE/DialysisPeriod_adj_for_cls3/BCV.pdf'))
plotBCV(y)
dev.off()

fit <- glmFit(y, design)
lrt <- glmLRT(fit)
topTags(lrt)
## Coefficient:  y$samples$group.L 
##         logFC     logCPM       LR       PValue        FDR
## CNTN6    3.7786132 3.10330152 21.39212 3.743073e-06 0.03543724
## AEBP1   -2.3168043 8.42174703 21.24597 4.039583e-06 0.03543724
## GAS2     1.6558515 1.15143183 19.73182 8.910536e-06 0.05211178
## ZNF365   2.4580046 1.34669809 17.99009 2.220576e-05 0.09740002
## CLGN     2.1798369 0.03694713 16.61941 4.568103e-05 0.16029474
## PLCH2   -2.1749391 5.87848619 15.81275 6.992968e-05 0.17761990
## CLDN9   -1.8921895 0.27936959 15.78759 7.086573e-05 0.17761990
## POGLUT2  0.8921113 3.98467369 15.37037 8.836299e-05 0.19379109
## ITGBL1  -2.6467560 2.95979699 14.35983 1.509894e-04 0.29434552
## TRABD2A  1.5402903 2.25125398 13.76882 2.067389e-04 0.36272348


o <- order(lrt$table$PValue)
cpm(y)[o[1:10],]
##     IWK001_T     IWK002_T     IWK005_T   IWK008_T    IWK009_T   IWK010_T     IWK012_T    IWK013_T   IWK015_T   IWK016_T     IWK017_T     IWK018_T   IWK019_T    IWK020_T     IWK021_T     IWK022_T
## CNTN6   20.8591529   0.00000000   0.00000000  1.0038744   0.0000000   0.049847   0.14306967 127.4317112 12.2531226 19.0851285   1.61933663    1.3030577  0.1618742   0.0000000   0.04981946   0.03547538
## AEBP1   28.1019144 835.96595049 223.44108684  2.9614293 373.3220805  35.540912 296.48803938  28.9898206 86.4419510  8.4510945 918.71785004 1286.9699962 74.7534977 581.5676828 493.16280502 402.57461261
## GAS2     0.3621381   1.57485344   1.66215474  1.3552304   0.9758779   1.395716   1.66914611   4.5174720 13.8804905  4.9583912   1.87502136    0.4510585  0.7769961   0.5274592   0.64765294   1.63186748
## ZNF365   5.0699330   0.16873430   0.09498027  0.9536806   0.0000000   0.249235   0.04768989   7.0524563 31.0476258  0.4365879   0.21307061    1.8042338  0.9064954   0.1318648   1.59422262   0.28380304
## CLGN     2.1728284   0.05624477   0.18996054  3.6139477   0.5855267   0.049847   0.04768989   3.3474793  4.5949210  1.9646456   0.08522824    0.2004704  0.3884980   0.9230535   0.09963891   0.07095076
## PLCH2    2.9695322 119.74510587  66.67615004  3.3127854   4.0498932 113.003153  26.18174901   8.2874487 13.3380345  3.3991487  62.38707419   26.0611550 23.6336307  26.7685526   4.88230679  59.81149073
## CLDN9    0.5794209   1.23738484   0.28494081  0.2007749   4.6354199   0.747705   0.85841800   0.0649996  0.4148193  0.5301425   1.36365189    0.7016465  0.7769961   0.7911887   1.29530588   0.70950760
## POGLUT2 18.4690417   8.94291772   9.97292842 33.6297909   8.5877253   5.832099   9.15645867  23.1073570 28.3672552 50.3011639  16.27859449    8.6202282 16.3492925  15.2084059  11.40865566  15.68011797
## ITGBL1   1.7382627   0.67493719   0.37992108  0.1505812   3.5131603   5.582864   0.19075956   1.4949907  1.0849119  0.6860667  45.72495260   65.2030050  2.5252373   6.1976452  12.80360045   0.67403222
## TRABD2A 21.8007119   0.61869242   0.42741122  8.0309948   1.3174351   0.199388   0.52458878  19.9548765  7.8815659  5.8003822   0.98012480    4.6108197  1.9424902   3.0328902   0.34873620   0.63855684
## IWK024_T    IWK025_T   IWK026_T    IWK028_T     IWK029_T     IWK048_T
## CNTN6     0.03730799   0.4314842  5.7698812   0.4393526   0.10721819   0.15149826
## AEBP1   430.98194392 489.8783858 20.1945843  95.5042683 151.32060878 618.98399847
## GAS2      1.23116379   1.7738795  2.5830282   0.5491907   1.14366072   1.21198605
## ZNF365    0.14923198   0.2876561  0.5031873   1.4828150   0.07147879   1.47710799
## CLGN      0.11192398   1.1026818  1.1405579   0.2196763   0.00000000   0.03787456
## PLCH2    49.99271164  13.7116088 20.9661382 455.2241978  17.90543812 148.35466695
## CLDN9     0.44769593   0.8150257  0.1677291   6.3156934   0.10721819   1.89372820
## POGLUT2   8.91661051  11.4103597 15.4646235   5.6017455  16.08272885   8.93839709
## ITGBL1    3.61887539   4.7942688  1.0063746   5.3820692   0.78626674   2.23459927
## TRABD2A   2.53694358   1.6300514 14.7937071   1.9221676   2.00140626   2.23459927
summary(decideTests(lrt))
## y$samples$group.L
## Down                 1
## NotSig             17543
## Up                   1


pdf(file=file.path(run.dir, 'result/DGE/DialysisPeriod_adj_for_cls3/MD.pdf'))
plotMD(lrt)
abline(h=c(-1, 1), col="blue")
dev.off()

library(EnhancedVolcano)

volcanodata <- as.data.table(topTags(lrt, n=nrow(lrt$table)), keep.rownames=TRUE)
setnames(volcanodata,
	old=paste('table.', c('rn', 'logFC', 'logCPM', 'LR', 'PValue', 'FDR'), sep=""),
	new = c("Gene", "logFC",  "logCPM", "LR", "pvalue", "FDR"))

pdf(file=file.path(run.dir, 'result/DGE/DialysisPeriod_adj_for_cls3/DialysisPeriod_Volcano.pdf'))
EnhancedVolcano(
		volcanodata, 
		x = 'logFC',
		y = 'pvalue',
		lab = volcanodata$Gene,
		pCutoff = 0.001,
		FCcutoff = 1,
		title = "Long vs. Short dialysis period edgeR")
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
				pCutoff = 0.001,
				FCcutoff = 1,
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

ggexport( curr.p, width=2.2, height=2.2, filename=file.path(run.dir, 'result/DGE/DialysisPeriod_adj_for_cls3/DialysisPeriod_Volcano.for.figure.2.2in.pdf'))
ggexport( curr.p, width=2.0, height=2.0, filename=file.path(run.dir, 'result/DGE/DialysisPeriod_adj_for_cls3/DialysisPeriod_Volcano.for.figure.2.0in.pdf'))


save( list=c('volcanodata', 'sample.info', 'lrt', 'y'),
	file=file.path(run.dir, 'result/DGE/DialysisPeriod_adj_for_cls3/DialysisPeriod_DE.Rdata'))

fwrite(
		volcanodata[,list(Gene, logFC, logCPM, LR, pvalue, FDR)],
		file = file.path(run.dir, 'result/DGE/DialysisPeriod_adj_for_cls3/DialysisPeriod_DGE_list.tsv'),
		sep = '\t')

volcanodata[FDR<0.25,]
nrow(volcanodata[FDR<0.25,])

GSEA.ranked.gene.list.dt <- fread(file.path(run.dir, '/result/GSEA/April02/Long_Short_dialysis_Period.MSigDB.gene_permutation.Gsea.1648900368712/ranked_gene_list_Long_versus_Short_1648900368712.tsv'))
GSEA.ranked.gene.list.dt.filt <- GSEA.ranked.gene.list.dt[SCORE > 0.78 | SCORE < -0.78,]

#genes.ls <- c('MCM4', 'MCM2', 'UBE2E1', 'GOT2', 'C1QBP')
genes.ls <- GSEA.ranked.gene.list.dt.filt$NAME

#volcanodata.sig <- volcanodata[FDR<0.2 | Gene%in%genes.ls,]
volcanodata.sig <- volcanodata[Gene%in%genes.ls,]

extracted.cpm <- cpm(y, log=FALSE)[volcanodata.sig$Gene,]
extracted.log.cpm <- cpm(y, log=TRUE)[volcanodata.sig$Gene,]

sig.genes.ls <- volcanodata.sig$Gene
names(sig.genes.ls) <- sig.genes.ls
volcanodata.sig
## > volcanodata[Gene%in%c('HIF1A', 'HIF1A-AS2', 'VHL', 'EPAS1'),]
#		Gene      logFC      logCPM         LR       pvalue        FDR adjust.method              comparison test
#1:      AEBP1 -2.3168043  8.42174703 21.2459724 4.039583e-06 0.03543724            BH dialysis_period_class.L  glm
#2:       CLGN  2.1798369  0.03694713 16.6194138 4.568103e-05 0.16029474            BH dialysis_period_class.L  glm
#3:    POGLUT2  0.8921113  3.98467369 15.3703725 8.836299e-05 0.19379109            BH dialysis_period_class.L  glm

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
		x = extracted.cpm.dt.long,
		y = extracted.log.cpm.dt.long,
		by = c('tumor.sample.id', 'gene'))

extracted.cpm.dt.long <- merge(
	x = sample.info[,list(tumor.sample.id, Histology, exprClass, years.dialysis)],
	y = extracted.cpm.dt.long,
	by = c('tumor.sample.id'))

extracted.cpm.dt.by.gene <- dcast(
		extracted.cpm.dt.long,
		formula = exprClass + Histology + tumor.sample.id  ~ gene,
		value.var = c('ct', 'log2.ct'))

#extracted.cpm.dt.by.gene[,Histology:=exprClass]

exprClass.colors <- c('brown2', 'aquamarine4')
names(exprClass.colors) <- c("ccRCC", "non-ccRCC")

extracted.cpm.dt.long[,exprClass.color:=exprClass.colors[exprClass]]

my_comparisons1 <- list(
		c("Long", "Short"))

curr.alpha <- 1.0
jitter.size <- 1.0

sig.genes.ls.split <- split(sig.genes.ls, cut(1:length(sig.genes.ls), breaks=23, labels=FALSE))

arranged.plots.ls <- lapply(
	X = sig.genes.ls.split,
	FUN = function( curr.split.genes ){
		boxplots.ls <- lapply(
			X = curr.split.genes,
			FUN = function( curr.gene ){
				curr.cpm.dt <- extracted.cpm.dt.long[gene==curr.gene,list(tumor.sample.id, ct, exprClass, exprClass.color)]
		
				ggboxplot(
					curr.cpm.dt,
					x = "exprClass", y = "ct", 
					add = "jitter",
					add.params = list(shape = 16, size = jitter.size, color = "exprClass"),
					title = curr.gene,
					ylab = "Counts per million",
					xlab = "Dialysis period class") +
						scale_color_manual(values = alpha(exprClass.colors, curr.alpha)) +
						stat_compare_means(comparisons = my_comparisons1, method="wilcox.test", vjust=0) +
						coord_cartesian(clip = 'off') +
						rremove('legend') +
						rremove('x.title') +
						rotate_x_text(45) +
						theme(plot.margin = margin(0.3, 0.1, 0.1, 0.1, "cm")) +
						theme(
							title = element_text(size=10, face='bold', color='black'),
							axis.title = element_text(size=8), axis.text = element_text(size=8))
			})
		if(length(boxplots.ls)==4){
			ggarrange(
					boxplots.ls[[1]],
					boxplots.ls[[2]],
					boxplots.ls[[3]],
					boxplots.ls[[4]],
					nrow = 2, ncol=2)
		}else if (length(boxplots.ls)==3){
		ggarrange(
				boxplots.ls[[1]],
				boxplots.ls[[2]],
				boxplots.ls[[3]],
				nrow = 2, ncol=2)
		}

})

arranged.scatterplots.ls <- lapply(
	X = sig.genes.ls.split,
	FUN = function( curr.split.genes ){
		plots.ls <- lapply(
			X = curr.split.genes,
			FUN = function( curr.gene ){
				curr.cpm.dt <- extracted.cpm.dt.long[gene==curr.gene,
					list(tumor.sample.id, ct, log2.ct, years.dialysis, Histology)]
				
				curr.plot <- ggscatter(
					data = curr.cpm.dt,
					x = 'years.dialysis',
					y  = "log2.ct",
					shape = 16,
					color = "Histology",
					add = "reg.line", 
					title = curr.gene,
					add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
					xlab = "Years of dialysis",
					ylab = "log2(cpm)")
				curr.x.range <- ggplot_build(curr.plot)$panel$ranges[[1]]$x.range
				curr.y.range <- ggplot_build(curr.plot)$panel$ranges[[1]]$y.range
				
				curr.plot +
						stat_cor(method = "pearson", label.x = curr.x.range[[2]]*0.15, label.y = curr.y.range*0.9, size=2.5) +
						rremove("legend") +
						scale_color_manual(values = alpha(exprClass.colors, curr.alpha)) +
						theme(title = element_text(size=10, face='bold', color='black'),
					legend.key.size = unit(0.6,"line"),
					legend.text = element_text(size=8),
				axis.title = element_text(size=8), axis.text = element_text(size=8))
				})
		if(length(plots.ls)==4){
			ggarrange(
					plots.ls[[1]],
					plots.ls[[2]],
					plots.ls[[3]],
					plots.ls[[4]],
					nrow = 2, ncol=2)
		}else if (length(plots.ls)==3){
			ggarrange(
					plots.ls[[1]],
					plots.ls[[2]],
					plots.ls[[3]],
					nrow = 2, ncol=2)
		}
		
	})
## VHL.HIF1A_AS3.scatter.p <- ggscatter(
##                 data = extracted.cpm.dt.by.gene[order(-exprClass),],
##                 x = 'ct_HIF1A-AS3',
##                 y  = "ct_VHL",
##                 shape = 16,
## #				color = 'black',
##                 color = "Histology",
## #				title = 'C',
##                 add = "reg.line", 
##                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
##                 xlab = "HIF1A-AS3 CPM",
##                 ylab = "VHL CPM") +
##         stat_cor(method = "pearson", label.x = 15, label.y = 45, size=2.5) +
##         rremove("legend") +
##         scale_color_manual(values = alpha(exprClass.colors, curr.alpha)) +
##         theme(title = element_text(size=10, face='bold', color='black'),
##                 legend.key.size = unit(0.6,"line"),
##                 legend.text = element_text(size=8),
##                 axis.title = element_text(size=8), axis.text = element_text(size=8))

DialysisPeriod.boxplots <- ggarrange(
		boxplots.ls[['CNTN6']],
		boxplots.ls[['AEBP1']],
		boxplots.ls[['ZNF365']],
		boxplots.ls[['CLGN']],
		boxplots.ls[['MCM4']],
		boxplots.ls[['MCM2']],
		boxplots.ls[['UBE2E1']],
		boxplots.ls[['GOT2']],
		boxplots.ls[['C1QBP']],
		ncol=3, nrow=3)

#ggexport(DialysisPeriod.boxplots, width=6, height=6, filename=file.path(run.dir, "result/DGE/DialysisPeriod_adj_for_cls3/DialysisPeriod_DGE.boxplots.pdf"))
ggexport(DialysisPeriod.boxplots, width=8, height=8, filename=file.path(run.dir, "result/DGE/DialysisPeriod_adj_for_cls3/DialysisPeriod_DGE.boxplots.added.pdf"))

ggexport(arranged.plots.ls, width=5, height=5, filename=file.path(run.dir, "result/DGE/DialysisPeriod_adj_for_cls3/DialysisPeriod_DGE.boxplots.all.sig.pdf"))


ggexport(arranged.scatterplots.ls, width=5, height=5, filename=file.path(run.dir, "result/DGE/DialysisPeriod_adj_for_cls3/DialysisPeriod_DGE.scatterplots.all.sig.pdf"))


curr.gene <- 'CNTN6'
#curr.gene <- 'AEBP1'
#curr.gene <- 'IGDCC4'

curr.cpm.dt <- data.table(
		log.cpm=cpm(y, log=TRUE)[curr.gene,],
		cpm=cpm(y, log=FALSE)[curr.gene,],
		sample=colnames(cpm(y)))

setkey(sample.info, tumor.sample.id)
curr.cpm.dt[,years.dialysis:=sample.info[sample,]$years.dialysis]

plot( x=curr.cpm.dt$years.dialysis, y=curr.cpm.dt$log.cpm, pch=19)

cor.test(x=curr.cpm.dt$log.cpm, y=curr.cpm.dt$years.dialysis)
## Pearson's product-moment correlation
## 
##         data:  curr.cpm and sample.info[curr.cpm.dt$sample, ]$years.dialysis
##         t = 3.8758, df = 20, p-value = 0.0009404
##         alternative hypothesis: true correlation is not equal to 0
##         95 percent confidence interval:
##         0.3223210 0.8435997
##         sample estimates:
##         cor 
##         0.654929 
