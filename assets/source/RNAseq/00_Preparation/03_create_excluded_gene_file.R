#!/usr/bin/env Rscript

# module use /usr/local/package/modulefiles/
# module load R/4.0.2
# . ~/.R-4.1.0_setup.sh

options(width=200)

library(data.table)
library(parallel)

ensembl_cache_dir <- "/home/tjohnson/reference/HMF_38/dbs/linx"
RNAseq_cg_pipeline_dir="/home/tjohnson/workspace/scripts/cancer_genomics_pipeline/v1/RNAseq"

genes.dt <- fread(paste(ensembl_cache_dir, "/ensembl_gene_data.csv", sep=""))

genes.excluded.dt <- rbind(
	genes.dt[substring(GeneName, 1, 5)=='MTATP',list(GeneId, GeneName)],
	genes.dt[substring(GeneName, 1, 4)=='MTND',list(GeneId, GeneName)],
	genes.dt[substring(GeneName, 1, 4)=='MTCO',list(GeneId, GeneName)],
	genes.dt[substring(GeneName, 1, 4)=='RNA5',list(GeneId, GeneName)])
nrow(genes.excluded.dt)
# 888

fwrite(x = genes.excluded.dt,
	file = file.path( RNAseq_cg_pipeline_dir, 'common_config_files/excluded_genes_20211116.csv'),
	col.names = TRUE, sep = ',')
