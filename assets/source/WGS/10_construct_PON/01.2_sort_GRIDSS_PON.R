#!/usr/bin/env Rscript

# module use /usr/local/package/modulefiles/
#module load R/4.0.2

library(data.table)

working_dir <- getwd()

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "~/workspace/runs", study.dir )

source( file.path( run.dir, "config_files/common_config.R" ) )

ref_data_path <- "/home/tjohnson/reference/HMF/38"
breakend_filename <- file.path(ref_data_path, 'dbs/gridss_pon/gridss_pon_single_breakend.JP.149x.bed')
breakpoint_filename <- file.path(ref_data_path, 'dbs/gridss_pon/gridss_pon_breakpoint.JP.149x.bedpe')

breakend_filename.unsorted <- file.path(ref_data_path, 'dbs/gridss_pon/gridss_pon_single_breakend.JP.149x.unsorted.bed')
breakpoint_filename.unsorted <- file.path(ref_data_path, 'dbs/gridss_pon/gridss_pon_breakpoint.JP.149x.unsorted.bedpe')

system(paste('mv ', breakend_filename, ' ', breakend_filename.unsorted, sep=''))
system(paste('mv ', breakpoint_filename, ' ', breakpoint_filename.unsorted, sep=''))

breakend_pon.dt <- fread(breakend_filename.unsorted)
breakpoint_pon.dt <- fread(breakpoint_filename.unsorted)

breakend_pon.dt <- breakend_pon.dt[order(V1, V2, V3),]
breakpoint_pon.dt <- breakpoint_pon.dt[order(V1, V2, V3, V4, V5, V6),]

fwrite(breakend_pon.dt,
	file = breakend_filename,
	sep = '\t',
	col.names = FALSE)
	
fwrite(breakpoint_pon.dt,
	file = breakpoint_filename,
	sep = '\t',
	col.names = FALSE)