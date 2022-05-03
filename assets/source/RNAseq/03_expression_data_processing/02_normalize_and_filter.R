args = commandArgs(trailingOnly=TRUE)
run.dir <- args[[1]]

options(width=350)

source( file.path(run.dir, "config_files/common_config.R") )
source( file.path(RNAseq_cg_pipeline_dir, "03_expression_data_processing/expression_processing_functions.R") )

library(data.table)
library(parallel)

print( paste("Attempting to load imported gene expression data for ", study.name, sep="") )
# process.tximport.data.ls has slots for
# sample.info, sample.ids.ls, imported.tx.expr.data, imported.gene.expr.data
load(file = file.path(run.dir, 'result/processed_expression_data/imported.tx.gene.expression.Rdata'))

normalized.expression.data.ls <- normalize_and_filter_counts(
	input.processed.tximport.data.ls = processed.tximport.data.ls,
	input.grouping.column = grouping.column,
	input.samples.to.remove = samples.to.exclude.from.normalization )

#normalized.expression.data.ls contains
#sample.ids.filtered, expr_filtered_DGEList, mat.all.samples.natural.cpm, mat.all.samples
print(paste('Saving normalized expression data for ', study.name, sep=""))

save( list=c('normalized.expression.data.ls'),
	file = file.path(run.dir, 'result/processed_expression_data/gene_expression.edgeR.normalized.Rdata'))
