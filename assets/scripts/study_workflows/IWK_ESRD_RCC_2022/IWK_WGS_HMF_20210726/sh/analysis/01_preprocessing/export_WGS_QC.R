#!/usr/bin/env Rscript

options(width=220)
working_dir <- getwd()
gpl_prefix <- basename(working_dir)

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "~/workspace/runs", study.dir )

source( file.path( run.dir, "config_files/common_config.R" ) )


library(data.table)
library(parallel)
library(openxlsx)

load(file = paste(run.dir, "/result_summaries/merged_file_consistency_summary.Rdata", sep=""))

sample.wgs.qc.info.dt <- sample_consistency_summary_merged[,
	list(sorted.bam.paired.read.total, sorted.bam.mapped.paired.read.total, UNMAPPED_READS, PERCENT_DUPLICATION),
	by=list(subject.id, sample.id)]

sample.wgs.qc.info.dt[,propn.unmapped:=UNMAPPED_READS/sorted.bam.mapped.paired.read.total]

setnames(sample.wgs.qc.info.dt,
	old = c('subject.id', 'sample.id',
		'sorted.bam.paired.read.total', 'sorted.bam.mapped.paired.read.total', 'UNMAPPED_READS', 'propn.unmapped', 'PERCENT_DUPLICATION'),
	new = c("Subject ID", "Sample ID",
		'Total Reads', 'Total mapped reads', 'Unmapped reads', 'Propn. unmapped', 'Propn. duplicates'))

wb <- createWorkbook(title = paste(study.name, " study: WGS quality", sep=""))

modifyBaseFont(wb, fontSize = 12, fontName = "Calibri")
addWorksheet(wb, sheetName = "WGS quality", gridLines = TRUE)

writeDataTable(wb, sheet = "WGS quality",
		startCol = 1, startRow = 1,
		x = sample.wgs.qc.info.dt,
		colNames = TRUE, rowNames = FALSE,
		withFilter = FALSE,
		bandedRows=FALSE,
		tableStyle = 'none')

saveWorkbook(wb,
		file = paste(run.dir, "/result/sample_summaries/WGS_QC_summary.xlsx", sep=""),
		overwrite = TRUE)