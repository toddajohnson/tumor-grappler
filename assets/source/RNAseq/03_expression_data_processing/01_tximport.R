args = commandArgs(trailingOnly=TRUE)
run.dir <- args[[1]]

library(data.table)

print( paste( 'Loading common scripts and function', sep='') )
source( file.path(run.dir, "config_files/common_config.R") )
source( file.path(RNAseq_cg_pipeline_dir, "03_expression_data_processing/expression_processing_functions.R") )

print( paste( 'Loading file processing information for ', study.name, sep='') )
load(file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))

print( paste( 'Running tximport for ', study.name, sep='') )
processed.tximport.data.ls <- tximport_processing(
	input.sample.info = sample.info.with.clinical.info.dt,
	input.read.length = 126)

print( paste( 'Saving imported transcript and gene expression data for ', study.name, sep='') )
save(list=c('processed.tximport.data.ls'),
	file = file.path(run.dir, 'result/processed_expression_data/imported.tx.gene.expression.Rdata'))