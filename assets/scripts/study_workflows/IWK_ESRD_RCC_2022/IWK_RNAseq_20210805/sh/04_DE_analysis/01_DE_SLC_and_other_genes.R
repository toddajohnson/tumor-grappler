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


load(file=file.path(isofox.dir, 'imported.tx.gene.expression.Rdata'))
#load(file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))
load(file = paste(wgs.run.dir, "/result/sample_summaries/clinical_data_with_colors.Rdata", sep=""))

#gene.mat <- fread(file.path(isofox.dir, "isofox.gene_expression_matrix.csv"))
#tx.mat <- fread(file.path(isofox.dir, "isofox.transcript_expression_matrix.csv"))

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


Young.dt <- read.xlsx(file.path(run.dir, 'result/adhoc/kidney_and_RCC_genes.xlsx'), sheet='Young_marker_list')
Liao.dt <- read.xlsx(file.path(run.dir, 'result/adhoc/kidney_and_RCC_genes.xlsx'), sheet='Liao_cluster_genes')
Lindren.dt <- read.xlsx(file.path(run.dir, 'result/adhoc/kidney_and_RCC_genes.xlsx'), sheet='Lindren_TableS4')

Young.dt <- as.data.table(Young.dt)
setnames(Young.dt, old='Positive_marker_mRNA', new='genes')

extract.genes <- function(curr.genes){
	strsplit(curr.genes, split='; ', fixed=TRUE)[[1]]
}
Young.genes.dt <- Young.dt[,list(gene=extract.genes(genes)), by=list(Cluster_ID, Alias)]
Young.gene.Clusters <- Young.genes.dt[,list(Clusters=paste(Cluster_ID, collapse=';'), Aliases=paste(Alias, collapse=';')), by=list(gene)]

Liao.dt <- as.data.table(Liao.dt)

kidney.cell.genes <- unique(Liao.dt$gene, Young.gene.Clusters$gene)

load(file=file.path(run.dir, 'result/TCGA/DGE/cls.3.ccRCC.Clus.B.chRCC/cls.3.ccRCC.Clus.B.chRCC.Clus.C_DE.Rdata'))
chRCC.top.20.genes <- volcanodata[FDR<1e-100 & Gene%in%rownames(cts),][order(-logFC),][1:20,]$Gene
chRCC.top.10.genes <- volcanodata[FDR<1e-100 & Gene%in%rownames(cts),][order(-logFC),][1:10,]$Gene
chRCC.top.5.genes <- volcanodata[FDR<1e-100 & Gene%in%rownames(cts),][order(-logFC),][1:5,]$Gene

cls.3.B.C.DGE.dt <- fread(file.path(run.dir, 'result/DGE/ccRCC_non_ccRCC/ccRCC_nonccRCC_DGE_list.tsv'))
cls.3.B.GSEA.dt <- fread(file.path(run.dir, 'result/GSEA/April02/cls3.B.C.MSigDB.gene_permutation.Gsea.1648900266798/gsea_report_for_B_1648900266798.tsv'))
cls.3.B.GSEA.dt <- cls.3.B.GSEA.dt[1:20,]
cls.3.B.GSEA.dt <- cls.3.B.GSEA.dt[which(cls.3.B.GSEA.dt[,'FWER p-val', with=FALSE]<0.25),]

cls.3.B.GSEA.sig.pathway.names <- cls.3.B.GSEA.dt$NAME
cls.3.B.MSigDB.GSEA.sig.pathway.ct <- length(cls.3.B.GSEA.sig.pathway.names)

cls.3.B.MSigDB.GSEA.sig.pathway.gene.info.dt.ls <- lapply(
		X = cls.3.B.GSEA.sig.pathway.names,
		FUN = function( curr.name ){
			print(paste('Extracting ', curr.name, sep=""))
			curr.dt <- fread(file.path(run.dir, 'result/GSEA/April02/cls3.B.C.MSigDB.gene_permutation.Gsea.1648900266798', paste(curr.name, '.tsv', sep='')))
			setnames(curr.dt,
					old=c('SYMBOL', 'TITLE', 'RANK IN GENE LIST', 'RANK METRIC SCORE', 'RUNNING ES'),
					new = c('gene', 'description', 'rank', 'rank.score', 'running.ES'))
			curr.dt[,geneset:=curr.name]
			curr.dt[which(curr.dt[,'CORE ENRICHMENT', with=FALSE]=='Yes'),list(gene, description, rank, rank.score, running.ES, geneset)]
		})

cls.3.B.MSigDB.GSEA.sig.pathway.gene.info.dt <- rbindlist(cls.3.B.MSigDB.GSEA.sig.pathway.gene.info.dt.ls)
cls.3.B.MSigDB.GSEA.Hypoxia.genes <- cls.3.B.MSigDB.GSEA.sig.pathway.gene.info.dt[geneset=='HALLMARK_HYPOXIA',]$gene
cls.3.B.MSigDB.GSEA.sig.pathway.gene.summary <- cls.3.B.MSigDB.GSEA.sig.pathway.gene.info.dt[,list(max.rank=max(rank), min.rank.score=min(rank.score), occurrences=.N), by =list(gene)]

cls.3.B.MSigDB.GSEA.sig.genes <- cls.3.B.MSigDB.GSEA.sig.pathway.gene.summary[occurrences>=0.2*cls.3.B.MSigDB.GSEA.sig.pathway.ct,][order(min.rank.score),]$gene

length(cls.3.B.MSigDB.GSEA.sig.genes)
#27
# 94
cls.3.B.MSigDB.top.20.genes <- cls.3.B.MSigDB.GSEA.sig.genes[1:20]
cls.3.B.MSigDB.top.10.genes <- cls.3.B.MSigDB.GSEA.sig.genes[1:10]


cls.3.C.GSEA.dt <- fread(file.path(run.dir, 'result/GSEA/April02/cls3.B.C.MSigDB.gene_permutation.Gsea.1648900266798/gsea_report_for_C_1648900266798.tsv'))
cls.3.C.GSEA.dt <- cls.3.C.GSEA.dt[1:20,]
cls.3.C.GSEA.dt <- cls.3.C.GSEA.dt[which(cls.3.C.GSEA.dt[,'FWER p-val', with=FALSE]<0.25),]

cls.3.C.GSEA.sig.pathway.names <- cls.3.C.GSEA.dt$NAME
cls.3.C.MSigDB.GSEA.sig.pathway.ct <- length(cls.3.C.GSEA.sig.pathway.names)

cls.3.C.MSigDB.GSEA.sig.pathway.gene.info.dt.ls <- lapply(
		X = cls.3.C.GSEA.sig.pathway.names,
		FUN = function( curr.name ){
			print(paste('Extracting ', curr.name, sep=""))
			curr.dt <- fread(file.path(run.dir, 'result/GSEA/April02/cls3.B.C.MSigDB.gene_permutation.Gsea.1648900266798', paste(curr.name, '.tsv', sep='')))
			setnames(curr.dt,
					old=c('SYMBOL', 'TITLE', 'RANK IN GENE LIST', 'RANK METRIC SCORE', 'RUNNING ES'),
					new = c('gene', 'description', 'rank', 'rank.score', 'running.ES'))
			curr.dt[,geneset:=curr.name]
			curr.dt[which(curr.dt[,'CORE ENRICHMENT', with=FALSE]=='Yes'),list(gene, description, rank, rank.score, running.ES, geneset)]
		})

cls.3.C.MSigDB.GSEA.sig.pathway.gene.info.dt <- rbindlist(cls.3.C.MSigDB.GSEA.sig.pathway.gene.info.dt.ls)
cls.3.C.MSigDB.GSEA.Hypoxia.genes <- cls.3.C.MSigDB.GSEA.sig.pathway.gene.info.dt[geneset=='HALLMARK_HYPOXIA',]$gene
cls.3.C.MSigDB.GSEA.sig.pathway.gene.summary <- cls.3.C.MSigDB.GSEA.sig.pathway.gene.info.dt[,list(max.rank=max(rank), min.rank.score=min(rank.score), occurrences=.N), by =list(gene)]

cls.3.C.MSigDB.GSEA.sig.genes <- cls.3.C.MSigDB.GSEA.sig.pathway.gene.summary[occurrences>=0.2*cls.3.C.MSigDB.GSEA.sig.pathway.ct,][order(min.rank.score),]$gene

length(cls.3.C.MSigDB.GSEA.sig.genes)
#120
cls.3.C.MSigDB.top.20.genes <- cls.3.C.MSigDB.GSEA.sig.genes[1:20]
cls.3.C.MSigDB.top.10.genes <- cls.3.C.MSigDB.GSEA.sig.genes[1:10]



cls.8.B.H.DGE.dt <- fread(file.path(run.dir, 'result/DGE/non-ccRCC_B_H/non-ccRCC_B_H_DGE_list.tsv'))
cls.8.H.BaderLab.GSEA.dt <- fread(file.path(run.dir, 'result/GSEA/April02/cld8.B.H.BaderLab_GO_All_no_GO.gene_permutation.Gsea.1648900461725/gsea_report_for_H_1648900461725.tsv'))
cls.8.H.BaderLab.GSEA.dt <- cls.8.H.BaderLab.GSEA.dt[1:20,]
cls.8.H.BaderLab.GSEA.dt <- cls.8.H.BaderLab.GSEA.dt[which(cls.8.H.BaderLab.GSEA.dt[,'FWER p-val', with=FALSE]<0.25),]

cls.8.H.BaderLab.GSEA.sig.pathway.names <- cls.8.H.BaderLab.GSEA.dt$NAME
cls.8.H.BaderLab.GSEA.sig.pathway.ct <- length(cls.8.H.BaderLab.GSEA.sig.pathway.names)

cls.8.H.BaderLab.GSEA.sig.pathway.gene.info.dt.ls <- lapply(
	X = cls.8.H.BaderLab.GSEA.sig.pathway.names,
	FUN = function( curr.name ){
		print(paste('Extracting ', curr.name, sep=""))
		curr.dt <- fread(file.path(run.dir, 'result/GSEA/April02/cld8.B.H.BaderLab_GO_All_no_GO.gene_permutation.Gsea.1648900461725', paste(curr.name, '.tsv', sep='')))
		setnames(curr.dt,
			old=c('SYMBOL', 'TITLE', 'RANK IN GENE LIST', 'RANK METRIC SCORE', 'RUNNING ES'),
			new = c('gene', 'description', 'rank', 'rank.score', 'running.ES'))
		curr.dt[,geneset:=curr.name]
		curr.dt[which(curr.dt[,'CORE ENRICHMENT', with=FALSE]=='Yes'),list(gene, description, rank, rank.score, running.ES, geneset)]
	})

cls.8.H.BaderLab.GSEA.sig.pathway.gene.info.dt <- rbindlist(cls.8.H.BaderLab.GSEA.sig.pathway.gene.info.dt.ls)
cls.8.H.BaderLab.GSEA.sig.pathway.gene.summary <- cls.8.H.BaderLab.GSEA.sig.pathway.gene.info.dt[,list(min.rank=min(rank), max.rank.score=max(rank.score), occurrences=.N), by =list(gene)]

cls.8.H.BaderLab.GSEA.sig.genes <- cls.8.H.BaderLab.GSEA.sig.pathway.gene.summary[occurrences>=0.2*cls.8.H.BaderLab.GSEA.sig.pathway.ct,][order(-max.rank.score),]$gene
length(cls.8.H.BaderLab.GSEA.sig.genes)
# 78
cls.8.H.BaderLab.top.20.genes <- cls.8.H.BaderLab.GSEA.sig.genes[1:20]
cls.8.H.BaderLab.top.10.genes <- cls.8.H.BaderLab.GSEA.sig.genes[1:10]


cls.8.B.BaderLab.GSEA.dt <- fread(file.path(run.dir, 'result/GSEA/April02/cld8.B.H.BaderLab_GO_All_no_GO.gene_permutation.Gsea.1648900461725/gsea_report_for_B_1648900461725.tsv'))
cls.8.B.BaderLab.GSEA.dt <- cls.8.B.BaderLab.GSEA.dt[1:20,]
cls.8.B.BaderLab.GSEA.dt <- cls.8.B.BaderLab.GSEA.dt[which(cls.8.B.BaderLab.GSEA.dt[,'FWER p-val', with=FALSE]<0.25),]

cls.8.B.BaderLab.GSEA.sig.pathway.names <- cls.8.B.BaderLab.GSEA.dt$NAME
cls.8.B.BaderLab.GSEA.sig.pathway.ct <- length(cls.8.B.BaderLab.GSEA.sig.pathway.names)

cls.8.B.BaderLab.GSEA.sig.pathway.gene.info.dt.ls <- lapply(
	X = cls.8.B.BaderLab.GSEA.sig.pathway.names,
	FUN = function( curr.name ){
		print(paste('Extracting ', curr.name, sep=""))
		curr.dt <- fread(file.path(run.dir, 'result/GSEA/April02/cld8.B.H.BaderLab_GO_All_no_GO.gene_permutation.Gsea.1648900461725', paste(curr.name, '.tsv', sep='')))
		setnames(curr.dt,
				old=c('SYMBOL', 'TITLE', 'RANK IN GENE LIST', 'RANK METRIC SCORE', 'RUNNING ES'),
				new = c('gene', 'description', 'rank', 'rank.score', 'running.ES'))
		curr.dt[,geneset:=curr.name]
		curr.dt[which(curr.dt[,'CORE ENRICHMENT', with=FALSE]=='Yes'),list(gene, description, rank, rank.score, running.ES, geneset)]
	})

cls.8.B.BaderLab.GSEA.sig.pathway.gene.info.dt <- rbindlist(cls.8.B.BaderLab.GSEA.sig.pathway.gene.info.dt.ls)
cls.8.B.BaderLab.GSEA.sig.pathway.gene.summary <- cls.8.B.BaderLab.GSEA.sig.pathway.gene.info.dt[,list(max.rank=max(rank), min.rank.score=min(rank.score), occurrences=.N), by =list(gene)]

cls.8.B.BaderLab.GSEA.sig.genes <- cls.8.B.BaderLab.GSEA.sig.pathway.gene.summary[occurrences>=0.2*cls.8.B.BaderLab.GSEA.sig.pathway.ct,][order(min.rank.score),]$gene
length(cls.8.H.BaderLab.GSEA.sig.genes)
# 78
cls.8.B.BaderLab.top.20.genes <- cls.8.B.BaderLab.GSEA.sig.genes[1:20]
cls.8.B.BaderLab.top.10.genes <- cls.8.B.BaderLab.GSEA.sig.genes[1:10]


cls.8.B.MSigDB.GSEA.dt <- fread(file.path(run.dir, 'result/GSEA/April02/cld8.B.H.MSigDB.gene_permutation.Gsea.1648900403910/gsea_report_for_B_1648900403910.tsv'))
cls.8.B.MSigDB.GSEA.dt <- cls.8.B.MSigDB.GSEA.dt[1:20,]
cls.8.B.MSigDB.GSEA.dt <- cls.8.B.MSigDB.GSEA.dt[which(cls.8.B.MSigDB.GSEA.dt[,'FWER p-val', with=FALSE]<0.25),]

cls.8.B.MSigDB.GSEA.sig.pathway.names <- cls.8.B.MSigDB.GSEA.dt$NAME
cls.8.B.MSigDB.GSEA.sig.pathway.ct <- length(cls.8.B.MSigDB.GSEA.sig.pathway.names)

cls.8.B.MSigDB.GSEA.sig.pathway.gene.info.dt.ls <- lapply(
		X = cls.8.B.MSigDB.GSEA.sig.pathway.names,
		FUN = function( curr.name ){
			print(paste('Extracting ', curr.name, sep=""))
			curr.dt <- fread(file.path(run.dir, 'result/GSEA/April02/cld8.B.H.MSigDB.gene_permutation.Gsea.1648900403910', paste(curr.name, '.tsv', sep='')))
			setnames(curr.dt,
					old=c('SYMBOL', 'TITLE', 'RANK IN GENE LIST', 'RANK METRIC SCORE', 'RUNNING ES'),
					new = c('gene', 'description', 'rank', 'rank.score', 'running.ES'))
			curr.dt[,geneset:=curr.name]
			curr.dt[which(curr.dt[,'CORE ENRICHMENT', with=FALSE]=='Yes'),list(gene, description, rank, rank.score, running.ES, geneset)]
		})

cls.8.B.MSigDB.GSEA.sig.pathway.gene.info.dt <- rbindlist(cls.8.B.MSigDB.GSEA.sig.pathway.gene.info.dt.ls)
cls.8.B.MSigDB.GSEA.sig.pathway.gene.summary <- cls.8.B.MSigDB.GSEA.sig.pathway.gene.info.dt[,list(max.rank=max(rank), min.rank.score=min(rank.score), occurrences=.N), by =list(gene)]

cls.8.B.MSigDB.GSEA.sig.genes <- cls.8.B.MSigDB.GSEA.sig.pathway.gene.summary[occurrences>=0.2*cls.8.B.MSigDB.GSEA.sig.pathway.ct,][order(min.rank.score),]$gene
length(cls.8.B.MSigDB.GSEA.sig.genes)
# 55
cls.8.B.MSigDB.top.20.genes <- cls.8.B.MSigDB.GSEA.sig.genes[1:20]
cls.8.B.MSigDB.top.10.genes <- cls.8.B.MSigDB.GSEA.sig.genes[1:10]


cls.8.H.MSigDB.GSEA.dt <- fread(file.path(run.dir, 'result/GSEA/April02/cld8.B.H.MSigDB.gene_permutation.Gsea.1648900403910/gsea_report_for_H_1648900403910.tsv'))
cls.8.H.MSigDB.GSEA.dt <- cls.8.H.MSigDB.GSEA.dt[1:20,]
cls.8.H.MSigDB.GSEA.dt <- cls.8.H.MSigDB.GSEA.dt[which(cls.8.H.MSigDB.GSEA.dt[,'FWER p-val', with=FALSE]<0.25),]

cls.8.H.MSigDB.GSEA.sig.pathway.names <- cls.8.H.MSigDB.GSEA.dt$NAME
cls.8.H.MSigDB.GSEA.sig.pathway.ct <- length(cls.8.H.MSigDB.GSEA.sig.pathway.names)

cls.8.H.MSigDB.GSEA.sig.pathway.gene.info.dt.ls <- lapply(
		X = cls.8.H.MSigDB.GSEA.sig.pathway.names,
		FUN = function( curr.name ){
			print(paste('Extracting ', curr.name, sep=""))
			curr.dt <- fread(file.path(run.dir, 'result/GSEA/April02/cld8.B.H.MSigDB.gene_permutation.Gsea.1648900403910', paste(curr.name, '.tsv', sep='')))
			setnames(curr.dt,
					old=c('SYMBOL', 'TITLE', 'RANK IN GENE LIST', 'RANK METRIC SCORE', 'RUNNING ES'),
					new = c('gene', 'description', 'rank', 'rank.score', 'running.ES'))
			curr.dt[,geneset:=curr.name]
			curr.dt[which(curr.dt[,'CORE ENRICHMENT', with=FALSE]=='Yes'),list(gene, description, rank, rank.score, running.ES, geneset)]
		})

cls.8.H.MSigDB.GSEA.sig.pathway.gene.info.dt <- rbindlist(cls.8.H.MSigDB.GSEA.sig.pathway.gene.info.dt.ls)
cls.8.H.MSigDB.GSEA.sig.pathway.gene.summary <- cls.8.H.MSigDB.GSEA.sig.pathway.gene.info.dt[,list(max.rank=max(rank), min.rank.score=min(rank.score), occurrences=.N), by =list(gene)]

cls.8.H.MSigDB.GSEA.sig.genes <- cls.8.H.MSigDB.GSEA.sig.pathway.gene.summary[occurrences>=0.2*cls.8.H.MSigDB.GSEA.sig.pathway.ct,][order(min.rank.score),]$gene
length(cls.8.H.MSigDB.GSEA.sig.genes)
# 68
cls.8.H.MSigDB.top.20.genes <- cls.8.H.MSigDB.GSEA.sig.genes[1:20]
cls.8.H.MSigDB.top.10.genes <- cls.8.H.MSigDB.GSEA.sig.genes[1:10]


HIF.pathway.genes <- c('HIF1A', 'EPAS1', 'HIF1A-AS3', 'ARNT', 'ARNT2', 'HIF3A', 'VHL', 'PBRM1')

genes.to.examine <- data.table(
	source = c(
			rep('HIF.related', length(HIF.pathway.genes)),
			rep('chRCC', length(chRCC.top.20.genes)),
			rep('Kidney', length(kidney.cell.genes)),
			rep('cls.3.B.MSigDB', length(cls.3.B.MSigDB.top.20.genes)),
			rep('cls.3.C.MSigDB', length(cls.3.C.MSigDB.top.20.genes)),
			rep('cls.8.B.MSigDB', length(cls.8.B.MSigDB.top.20.genes)),
			rep('cls.8.H.MSigDB', length(cls.8.H.MSigDB.top.20.genes)),
			rep('cls.8.B.BaderLab', length(cls.8.B.BaderLab.top.20.genes)),
			rep('cls.8.H.BaderLab', length(cls.8.H.BaderLab.top.20.genes))),
	gene = c(
			HIF.pathway.genes,
			chRCC.top.20.genes,
			kidney.cell.genes,
			cls.3.B.MSigDB.top.20.genes,
			cls.3.C.MSigDB.top.20.genes,
			cls.8.B.MSigDB.top.20.genes,
			cls.8.H.MSigDB.top.20.genes,
			cls.8.B.BaderLab.top.20.genes,
			cls.8.H.BaderLab.top.20.genes))

genes.to.examine.10 <- data.table(
		source = c(
				rep('HIF.related', length(HIF.pathway.genes)),
				rep('chRCC', length(chRCC.top.5.genes)),
				rep('Kidney', length(kidney.cell.genes)),
				rep('cls.3.B.MSigDB', length(cls.3.B.MSigDB.top.10.genes)),
				rep('cls.3.C.MSigDB', length(cls.3.C.MSigDB.top.10.genes)),
				rep('cls.8.B.MSigDB', length(cls.8.B.MSigDB.top.10.genes)),
				rep('cls.8.H.MSigDB', length(cls.8.H.MSigDB.top.10.genes)),
				rep('cls.8.B.BaderLab', length(cls.8.B.BaderLab.top.10.genes)),
				rep('cls.8.H.BaderLab', length(cls.8.H.BaderLab.top.10.genes))),
		gene = c(
				HIF.pathway.genes,
				chRCC.top.5.genes,
				kidney.cell.genes,
				cls.3.B.MSigDB.top.10.genes,
				cls.3.C.MSigDB.top.10.genes,
				cls.8.B.MSigDB.top.10.genes,
				cls.8.H.MSigDB.top.10.genes,
				cls.8.B.BaderLab.top.10.genes,
				cls.8.H.BaderLab.top.10.genes))

save( list=c('genes.to.examine.10', 'genes.to.examine', 'kidney.cell.genes', 'chRCC.top.20.genes', 'HIF.pathway.genes',
	'cls.3.B.MSigDB.GSEA.sig.genes', 'cls.3.C.MSigDB.GSEA.sig.genes', 
	'cls.8.B.BaderLab.GSEA.sig.genes', 'cls.8.H.BaderLab.GSEA.sig.genes',
	'cls.8.B.MSigDB.GSEA.sig.genes', 'cls.8.H.MSigDB.GSEA.sig.genes'),
	file = file.path(run.dir, 'result/adhoc/selected_genes/GSEA.extracted.genes.Rdata'))


unique.genes <- unique(genes.to.examine$gene)
length(unique.genes)
# 166

unique.genes.10 <- unique(genes.to.examine.10$gene)
length(unique.genes.10)
#96


## sample with high mitochondrial RNA content, 2 no tumor samples
samples.to.remove <- c('IWK006_T', 'IWK039_T', 'IWK049_T')

cols.to.keep <- colnames(cts)
cols.to.keep <- cols.to.keep[which(!cols.to.keep%in%samples.to.remove)]

clinical.data <- clinical.data[tumor.sample.id%in%cols.to.keep,]
setkey(clinical.data, tumor.sample.id)

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

groups <- clinical.data$Histology
names(groups) <- clinical.data$tumor.sample.id
groups <- groups[colnames(cts)]
groups.f <- factor(x = groups, ordered = FALSE)

y <- DGEList(cts, group=groups.f)

y <- scaleOffset(y, normMat)
# filtering

y.all <- y
all.natural.cpm <- cpm(y, log=FALSE)
all.log.cpm <- cpm(y, log=TRUE)

keep <- filterByExpr(y)
#y <- y[keep,,keep.lib.sizes=FALSE]
y <- y[keep,,keep.lib.sizes=FALSE]


unique.genes.missing <- unique.genes[which(!unique.genes%in%names(keep[keep==TRUE]))]
#unique.genes.missing
#[1] "PALS1"
unique.genes.non.missing <- unique.genes[which(!unique.genes%in%unique.genes.missing)]

# TMM normalization and log transformation
mat.all.samples.natural.cpm <- cpm(y, log=FALSE)
mat.all.samples <- cpm(y, log=TRUE)

# remove low-variance genes
#sds.all.samples <- apply(mat.all.samples , 1, sd)
#mat.all.samples.filt <- mat.all.samples[sds.all.samples > 1.5, ]

#mat.filt <- mat.all.samples[,colnames(mat.all.samples)[which(!colnames(mat.all.samples)%in%c('IWK006_T'))]]

#sds.filt <- apply(mat.filt , 1, sd)
#mat.filt <- mat.filt[sds.filt > 1.5, ]

mat.filt <- mat.all.samples[unique.genes.non.missing,]

# scale per gene
mat.filt <- mat.filt %>% t %>% scale %>% t

# clustering of tumors
cluster_columns <- mat.filt %>%
		t %>%
		dist %>%
		hclust(method="ward.D2")

cluster_rows <- mat.filt %>%
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

clinical.data[,DialysisPeriod:=cut(years.dialysis, breaks=c(0,14,30), labels=c('0-14 years', '15-29 years'))]

# prepare data for histology annotation
his <- clinical.data[,list(id=tumor.sample.id, Histology, DialysisPeriod)]
his[,ccRCC:=ifelse(Histology=='ccRCC', TRUE, FALSE)]
his[,ACD.RCC:=ifelse(Histology=='ACD-RCC', TRUE, FALSE)]
his[,pRCC:=ifelse(Histology=='pRCC', TRUE, FALSE)]
his[,chRCC:=ifelse(Histology=='chRCC', TRUE, FALSE)]
his[,ccpRCC:=ifelse(Histology=='ccpRCC', TRUE, FALSE)]

setnames(his,
		old = c('ACD.RCC'),
		new = c('ACD-RCC'))

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
		ccRCC = his$ccRCC,
		`ACD-RCC` = his$"ACD-RCC",
		pRCC = his$pRCC,
		ccpRCC = his$ccpRCC,
		chRCC = his$chRCC,
		
		col = list(
				'Cluster (k=3)' = colors.3.ls,
				'Cluster (k=8)' = k.colors.ls,
				'Dialysis period' = bluered,
				ccRCC = bw,
				`ACD-RCC` = bw,
				pRCC = bw,
				ccpRCC = bw,
				chRCC = bw),
		annotation_name_gp = gpar(fontsize = 8),
		show_legend = FALSE,
		simple_anno_size = unit(0.25, "cm"))


# draw heatmap
#col <- colorRamp2(c(-.5, 0, .5), c("#377EB8", "white", "#E41A1C"))
col <- colorRamp2(c(-.5, 0, .5), c("#377EB8", "white", "#E41A1C"))
col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))


h <- Heatmap(
		mat.filt,
		name = "Expression",
		heatmap_legend_param = list(title_gp = gpar(fontsize=9), labels_gp = gpar(fontsize=8)),
		col = col,
		show_row_names=T,
		row_names_gp = gpar(fontsize = 2),
		row_split = 6,
		cluster_rows = cluster_rows,
		cluster_columns = cluster_columns,
#		clustering_method_rows="ward.D2",
		top_annotation = top_annotation,
		column_names_gp = gpar(fontsize = 9),
		height = unit(7.3, "cm")
)


lgd_list = list(
		Legend(labels = names(colors.3.ls), title = "Cluster (k=3)", title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 8), legend_gp = gpar(fill=colors.3.ls)),
		Legend(labels = names(k.colors.ls), title = "Cluster (k=8)", title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 8), , legend_gp = gpar(fill=k.colors.ls)),
		Legend(labels = c("0-14 years", "15-29 years"), title = "Dialysis period", title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 8), , legend_gp = gpar(fill = bluered)),
		Legend(labels = c('TRUE', 'FALSE'), title = "Histology", title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 8), legend_gp = gpar(fill = bw)))


#pdf( file.path(run.dir, "result/adhoc/selected_genes/figures/heatmap_twolevel_cls.pdf"), width=7, height=5)
pdf( file.path(run.dir, "result/adhoc/selected_genes/figures/heatmap_twolevel_top.20.genes.each.source.cls.pdf"), width=7, height=5)
draw(
		h,
		annotation_legend_side = "left",
		annotation_legend_list = lgd_list,
		padding = unit(c(2, 2, 2, 8), "mm")
)
dev.off()

save( list=c('y', 'mat.all.samples.natural.cpm', 'mat.all.samples', 'mat.filt', 'his', 'cls', 'cls8', 'cluster_rows', 'genes.to.examine',
		'unique.genes', 'unique.genes.non.missing'),
		file = paste(run.dir, "/result/adhoc/selected_genes/selected_genes.gene_expression.edgeR.normalized.20.top.genes.Rdata", sep=""))


##
unique.genes.10.missing <- unique.genes[which(!unique.genes.10%in%names(keep[which(keep==TRUE)]))]
#unique.genes.missing
#[1] "PALS1"
unique.genes.non.missing <- unique.genes.10[which(!unique.genes.10%in%unique.genes.10.missing)]

mat.filt <- mat.all.samples[unique.genes.non.missing,]

# scale per gene
mat.filt <- mat.filt %>% t %>% scale %>% t

mat.filt <- mat.filt[gene.to.keep,]
# clustering of tumors
cluster_columns <- mat.filt %>%
		t %>%
		dist %>%
		hclust(method="ward.D2")

cluster_rows <- mat.filt %>%
		dist %>%
		hclust(method="ward.D2")

rows.cls <- cutree(cluster_rows, k=4)  # cut dendrogram at k = 5
rows.n <- names(rows.cls)
rows.colors.4.ls <- c("A1" = "grey", "2" = "black", "3" = "grey", "4" = "black")


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

clinical.data[,DialysisPeriod:=cut(years.dialysis, breaks=c(0,14,30), labels=c('0-14 years', '15-29 years'))]

# prepare data for histology annotation
his <- clinical.data[,list(id=tumor.sample.id, Histology, DialysisPeriod)]
his[,ccRCC:=ifelse(Histology=='ccRCC', TRUE, FALSE)]
his[,ACD.RCC:=ifelse(Histology=='ACD-RCC', TRUE, FALSE)]
his[,pRCC:=ifelse(Histology=='pRCC', TRUE, FALSE)]
his[,chRCC:=ifelse(Histology=='chRCC', TRUE, FALSE)]
his[,ccpRCC:=ifelse(Histology=='ccpRCC', TRUE, FALSE)]

setnames(his,
		old = c('ACD.RCC'),
		new = c('ACD-RCC'))

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
		ccRCC = his$ccRCC,
		`ACD-RCC` = his$"ACD-RCC",
		pRCC = his$pRCC,
		ccpRCC = his$ccpRCC,
		chRCC = his$chRCC,
		
		col = list(
				'Cluster (k=3)' = colors.3.ls,
				'Cluster (k=8)' = k.colors.ls,
				'Dialysis period' = bluered,
				ccRCC = bw,
				`ACD-RCC` = bw,
				pRCC = bw,
				ccpRCC = bw,
				chRCC = bw),
		annotation_name_gp = gpar(fontsize = 8),
		show_legend = FALSE,
		simple_anno_size = unit(0.25, "cm"))

# draw heatmap
col <- colorRamp2(c(-.5, 0, .5), c("#377EB8", "white", "#E41A1C"))
col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
h <- Heatmap(
		mat.filt,
		name = "Expression",
		heatmap_legend_param = list(title_gp = gpar(fontsize=9), labels_gp = gpar(fontsize=8)),
		col = col,
#		col = col_fun,
		show_row_names=T,
		row_names_gp = gpar(fontsize = 3),
		cluster_rows = cluster_rows,
		cluster_columns = cluster_columns,
#		clustering_method_rows="ward.D2",
		top_annotation = top_annotation,
		row_split = 6,
		column_names_gp = gpar(fontsize = 9),
		height = unit(7.3, "cm")
)


lgd_list = list(
		Legend(labels = names(colors.3.ls), title = "Cluster (k=3)", title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 8), legend_gp = gpar(fill=colors.3.ls)),
		Legend(labels = names(k.colors.ls), title = "Cluster (k=8)", title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 8), , legend_gp = gpar(fill=k.colors.ls)),
		Legend(labels = c("0-14 years", "15-29 years"), title = "Dialysis period", title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 8), , legend_gp = gpar(fill = bluered)),
		Legend(labels = c('TRUE', 'FALSE'), title = "Histology", title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 8), legend_gp = gpar(fill = bw)))


#pdf( file.path(run.dir, "result/adhoc/selected_genes/figures/heatmap_twolevel_top.10.genes.each.source.cls.pdf"), width=7, height=5)
pdf( file.path(run.dir, "result/adhoc/selected_genes/figures/heatmap_twolevel_top.10.genes.each.source.clss.pdf"), width=7, height=5)
draw(
		h,
		annotation_legend_side = "left",
		annotation_legend_list = lgd_list,
		padding = unit(c(2, 2, 2, 8), "mm")
)
dev.off()


save( list=c('y', 'mat.all.samples.natural.cpm', 'mat.all.samples', 'mat.filt', 'his', 'cls', 'cls8', 'cluster_rows', 'genes.to.examine',
		'unique.genes.10', 'unique.genes.non.missing'),
		file = paste(run.dir, "/result/adhoc/selected_genes/selected_genes.gene_expression.edgeR.normalized.10.top.genes.Rdata", sep=""))

which(cluster_rows$labels%in%cls.8.B.MSigDB.top.20.genes)

HIF.pathway.genes,
chRCC.top.20.genes,
kidney.cell.genes,
top.20.cls.3.B.GSEA.Hypoxia.genes,
cls.8.B.MSigDB.top.20.genes,
cls.8.B.BaderLab.top.20.genes,
cls.8.H.BaderLab.top.20.genes))

which(cluster_rows$labels%in%cls.8.H.BaderLab.top.20.genes)