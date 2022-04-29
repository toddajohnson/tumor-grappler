#!/usr/bin/env Rscript

# . ~/.R-4.1.0_setup.sh

options(width=350)
working_dir <- getwd()

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "~/workspace/runs", study.dir )

suppressPackageStartupMessages({
			library(data.table)
			library(parallel)
			library(maftools)
			library(sigminer)})

library(stringr)
library(openxlsx)

## arguments that could be added to common_config.R
remove.copy.number.noise.samples <- TRUE
tumor.depth.cutoff <- 10

source( file.path(run.dir, "config_files/common_config.R") )

HOME.dir <- "/home/tjohnson"

fig.dir <- file.path( run.dir, "result", "SigProfiler", "figures")

sig.type.ls <- c("SBS", "DBS", "ID", "CN_S_40", "CN_S_48", "CN_W")
names(sig.type.ls) <- sig.type.ls

sig.class.ls <- c("SBS", "DBS", "ID", "CN_S", "CN_S", "CN_W")
names(sig.class.ls) <- sig.type.ls

sig.mode.ls <- c("SBS", "DBS", "ID", "copynumber", "copynumber", "copynumber")
names(sig.mode.ls) <- sig.type.ls

curr.categorical.annotation.col <- 'Histology'
curr.numerical.annotation.col <- 'years.dialysis'

load(file = file.path( run.dir, "result/sample_summaries", "clinical_data_with_colors.Rdata" ))
clinical.data$cat.annotation <- clinical.data[[curr.categorical.annotation.col]]
clinical.data$num.annotation <- clinical.data[[curr.numerical.annotation.col]]

SBS.denovo.dt <- fread(file.path(run.dir, "result/SigProfiler", paste(study.name, "_mutation_signature_analysis/results/SBS/SBS96/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt", sep="")))
setnames(SBS.denovo.dt, old='Samples', new='tumor.sample.id')
SBS.denovo.dt <- merge(
		x = clinical.data[,list(tumor.sample.id, gender, purity, ploidy, cat.annotation, num.annotation)],
		y = SBS.denovo.dt,
		by = c('tumor.sample.id'))

SBS.dt <- fread(file.path(run.dir, "result/SigProfiler/", paste(study.name, "_mutation_signature_analysis/results/SBS/SBS96/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities_refit.txt", sep="")))
setnames(SBS.dt, old='Samples', new='tumor.sample.id')
SBS.dt <- merge(
		x = clinical.data[,list(tumor.sample.id, gender, purity, ploidy, cat.annotation, num.annotation)],
		y = SBS.dt,
		by = c('tumor.sample.id'))

sig.dt <- SBS.dt
sig.denovo.dt <- SBS.denovo.dt
save( list=c('sig.denovo.dt', 'sig.dt'),
		file = file.path(run.dir, "result/SigProfiler/", paste(study.name, "_mutation_signature_analysis/results/SBS.refit.activities.Rdata", sep="")))


DBS.denovo.dt <- fread(file.path(run.dir, "result/SigProfiler/", paste(study.name, "_mutation_signature_analysis/results/DBS/DBS78/DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Activities/DBS78_De-Novo_Activities_refit.txt", sep="")))
setnames(DBS.denovo.dt, old='Samples', new='tumor.sample.id')
DBS.denovo.dt <- merge(
		x = clinical.data[,list(tumor.sample.id, gender, purity, ploidy, cat.annotation, num.annotation)],
		y = DBS.denovo.dt,
		by = c('tumor.sample.id'))

DBS.dt <- fread(file.path(run.dir, "result/SigProfiler/", paste(study.name, "_mutation_signature_analysis/results/DBS/DBS78/DBS78/Suggested_Solution/COSMIC_DBS78_Decomposed_Solution/Activities/COSMIC_DBS78_Activities_refit.txt", sep="")))
setnames(DBS.dt, old='Samples', new='tumor.sample.id')
DBS.dt <- merge(
		x = clinical.data[,list(tumor.sample.id, gender, purity, ploidy, cat.annotation, num.annotation)],
		y = DBS.dt,
		by = c('tumor.sample.id'))

sig.dt <- DBS.dt
sig.denovo.dt <- DBS.denovo.dt
save( list=c('sig.denovo.dt', 'sig.dt'),
		file = file.path(run.dir, "result/SigProfiler/", paste(study.name, "_mutation_signature_analysis/results/DBS.refit.activities.Rdata", sep="")))


ID.denovo.dt <- fread(file.path(run.dir, "result/SigProfiler/", paste(study.name, "_mutation_signature_analysis/results/ID/ID83/ID83/Suggested_Solution/ID83_De-Novo_Solution/Activities/ID83_De-Novo_Activities_refit.txt", sep="")))
setnames(ID.denovo.dt, old='Samples', new='tumor.sample.id')
ID.denovo.dt <- merge(
		x = clinical.data[,list(tumor.sample.id, gender, purity, ploidy, cat.annotation, num.annotation)],
		y = ID.denovo.dt,
		by = c('tumor.sample.id'))

ID.dt <- fread(file.path(run.dir, "result/SigProfiler/", paste(study.name, "_mutation_signature_analysis/results/ID/ID83/ID83/Suggested_Solution/COSMIC_ID83_Decomposed_Solution/Activities/COSMIC_ID83_Activities_refit.txt", sep="")))
setnames(ID.dt, old='Samples', new='tumor.sample.id')
ID.dt <- merge(
		x = clinical.data[,list(tumor.sample.id, gender, purity, ploidy, cat.annotation, num.annotation)],
		y = ID.dt,
		by = c('tumor.sample.id'))

sig.dt <- ID.dt
sig.denovo.dt <- ID.denovo.dt
save( list=c('sig.denovo.dt', 'sig.dt'),
		file = file.path(run.dir, "result/SigProfiler/", paste(study.name, "_mutation_signature_analysis/results/ID.refit.activities.Rdata", sep="")))


load( file = file.path(run.dir, "result/SigProfiler/", paste(study.name, "_mutation_signature_analysis/results/SBS.refit.activities.Rdata", sep="")))
load( file = file.path(run.dir, "result/SigProfiler/", paste(study.name, "_mutation_signature_analysis/results/DBS.refit.activities.Rdata", sep="")))
load( file = file.path(run.dir, "result/SigProfiler/", paste(study.name, "_mutation_signature_analysis/results/ID.refit.activities.Rdata", sep="")))

ID.cols <- copy(colnames(ID.dt))
ID.cols <- ID.cols[grep('ID', ID.cols)]
merged.COSMIC.dt <- merge(
		x = SBS.dt,
		y = ID.dt[,c('tumor.sample.id', ID.cols), with=FALSE],
		by = c('tumor.sample.id'))

DBS.cols <- copy(colnames(DBS.dt))
DBS.cols <- DBS.cols[grep('DBS', DBS.cols)]
merged.COSMIC.dt <- merge(
		x = merged.COSMIC.dt,
		y = DBS.dt[,c('tumor.sample.id', DBS.cols), with=FALSE],
		by = c('tumor.sample.id'))

save(list=c('merged.COSMIC.dt'),
		file = file.path(run.dir, "result/SigProfiler/", paste(study.name, "_mutation_signature_analysis/results/merged.COSMIC.refit.activities.Rdata", sep="")))

wb <- createWorkbook(title = paste(study.name, " study: SBS, DBS, and ID COSMIC signature activities", sep=""))
addWorksheet(wb, sheetName = "BTC COSMIC signatures", gridLines = TRUE)
writeDataTable(wb,
		sheet = "BTC COSMIC signatures",
		startCol = 1, startRow = 1,
		x = merged.COSMIC.dt)

saveWorkbook(wb,
		file = paste(run.dir, "/result/SigProfiler/", study.name, "_mutation_signature_analysis/results/merged.COSMIC.refit.activities.xlsx", sep=""),
		overwrite = TRUE)


merged.COSMIC.dt.long <- melt(
	data = merged.COSMIC.dt,
	id.vars = c('tumor.sample.id', 'gender', 'purity', 'ploidy', 'cat.annotation', 'num.annotation'),
	variable.name = 'sig.class',
	value.name = 'mutation.ct',
	variable.factor = FALSE,
	value.factor = FALSE)

genome.length <- 2900000000
genome.length.Mb <- genome.length/1e6
merged.COSMIC.dt.long[,list(median.mutation.ct=median(ifelse(mutation.ct==0, as.numeric(NA), as.numeric(mutation.ct)), na.rm=TRUE)), by=list(sig.class)][,
	list(mutations.per.Mb=median.mutation.ct/genome.length.Mb), by=list(sig.class)]

sig.class mutations.per.Mb
1:      SBS1      0.048275862
2:      SBS4      0.302758621
3:      SBS5      0.250000000
4:     SBS12      0.347758621
5:     SBS23      0.059310345
6:     SBS40      1.110344828
7:       ID1      0.010172414
8:       ID3      0.048275862
9:       ID5      0.142758621
10:       ID8      0.045172414
11:      ID12      0.077931034
12:      DBS2      0.004137931
13:      DBS4      0.006551724
14:      DBS9      0.012068966




wb2 <- createWorkbook(title = paste(study.name, " study: SBS, DBS, and ID COSMIC signature TMB", sep=""))
addWorksheet(wb2, sheetName = "BTC COSMIC TMB", gridLines = TRUE)
writeDataTable(wb2,
		sheet = "BTC COSMIC TMB",
		startCol = 1, startRow = 1,
		x = merged.COSMIC.dt.long[,list(median.mutation.ct=median(ifelse(mutation.ct==0, as.numeric(NA), as.numeric(mutation.ct)), na.rm=TRUE)), by=list(sig.class)][,
				list(mutations.per.Mb=median.mutation.ct/genome.length.Mb), by=list(sig.class)])


addWorksheet(wb2, sheetName = "BTC COSMIC subtype TMB", gridLines = TRUE)
writeDataTable(wb2,
		sheet = "BTC COSMIC subtype TMB",
		startCol = 1, startRow = 1,
		x = merged.COSMIC.dt.long[,list(median.mutation.ct=median(ifelse(mutation.ct==0, as.numeric(NA), as.numeric(mutation.ct)), na.rm=TRUE)), by=list(Histology=cat.annotation, sig.class)][,
				list(mutations.per.Mb=median.mutation.ct/genome.length.Mb), by=list(Histology, sig.class)])

saveWorkbook(wb2,
		file = paste(run.dir, "/result/SigProfiler/", study.name, "_mutation_signature_analysis/results/merged.COSMIC.refit.activities.TMB.xlsx", sep=""),
		overwrite = TRUE)
