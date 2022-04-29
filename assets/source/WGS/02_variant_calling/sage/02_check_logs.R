#!/usr/bin/env Rscript

# module use /usr/local/package/modulefiles/
# module load R/4.0.2

options(width=200)

working_dir <- getwd()
gpl_prefix <- basename(working_dir)

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "/home/tjohnson/workspace/runs", study.dir )

source( file.path( run.dir, "config_files/common_config.R" ) )

suppressPackageStartupMessages({
	library(data.table)
	library(parallel)})

load(file = paste(run.dir, "/config_files/sage_run_info.Rdata", sep=""))

base.log.dir <- paste(run.dir, "/log/02_variant_calling/sage", sep="")

germline.log.dir <- paste(base.log.dir, "/germline/01_sage_combined_germline_calls", sep="")
somatic.log.dir <- paste(base.log.dir, "/somatic/01_sage_combined_somatic_calls", sep="")

germline.log.files <- list.files(germline.log.dir)
germline.previous.subdir.idxs <- grep('previous', germline.log.files, fixed=TRUE)
if (length(germline.previous.subdir.idxs)>0){
	germline.log.files <- germline.log.files[-germline.previous.subdir.idxs]
}


somatic.log.files <- list.files(somatic.log.dir)
somatic.previous.subdir.idxs <- grep('previous', somatic.log.files, fixed=TRUE)
if (length(somatic.previous.subdir.idxs)>0){
	somatic.log.files <- somatic.log.files[-somatic.previous.subdir.idxs]
}

setwd(germline.log.dir)

germline.run.indexes <- as.integer(unlist(lapply(X=germline.log.files, FUN=function(curr.file){strsplit(curr.file, split=".", fixed=TRUE)[[1]][[3]]})))
names(germline.log.files) <- germline.run.indexes

germline.completion.status.ls <- lapply(
	X = germline.run.indexes,
	FUN = function(curr.run.index){
		curr.file <- germline.log.files[[curr.run.index]]
		curr.log <- scan(curr.file, what="character", sep="\n")
		data.table(
			sage.run.index = curr.run.index,
			last.log.line = curr.log[[length(curr.log)]],
			completion.status = grep("Completed in", curr.log[[length(curr.log)]], fixed=TRUE))
})
germline.completion.status.dt <- rbindlist(germline.completion.status.ls)
germline.completion.status.dt <- germline.completion.status.dt[order(sage.run.index),]

sage.germline.run.completion.info.dt <- merge(
	x = sage.germline.run.info.dt[,list(sage.run.index=as.integer(1:nrow(sage.germline.run.info.dt)),
		subject.id, subject.index, tumor.name, reference.name)],
	y = germline.completion.status.dt,
	by = c("sage.run.index"),
	all.x = TRUE)


setwd(somatic.log.dir)

somatic.run.indexes <- as.integer(unlist(lapply(X=somatic.log.files, FUN=function(curr.file){strsplit(curr.file, split=".", fixed=TRUE)[[1]][[3]]})))
names(somatic.log.files) <- somatic.run.indexes

somatic.completion.status.ls <- lapply(
	X = somatic.run.indexes,
	FUN = function(curr.run.index){
		curr.file <- somatic.log.files[[curr.run.index]]
		curr.log <- scan(curr.file, what="character", sep="\n")
		data.table(
			sage.run.index = curr.run.index,
			last.log.line = curr.log[[length(curr.log)]],
			completion.status = grep("Completed in", curr.log[[length(curr.log)]], fixed=TRUE))
	})
somatic.completion.status.dt <- rbindlist(somatic.completion.status.ls)
somatic.completion.status.dt <- somatic.completion.status.dt[order(sage.run.index),]

sage.somatic.run.completion.info.dt <- merge(
		x = sage.somatic.run.info.dt[,list(sage.run.index=as.integer(1:nrow(sage.somatic.run.info.dt)),
						subject.id, subject.index, tumor.name, reference.name)],
		y = somatic.completion.status.dt,
		by = c("sage.run.index"),
		all.x = TRUE)

print( "Table of SAGE germline calling status" )
print( table(sage.germline.run.completion.info.dt$completion.status) )

print( "Table of SAGE somatic calling status" )
print( table(sage.somatic.run.completion.info.dt$completion.status) )

save( list=c("germline.completion.status.dt", "sage.somatic.run.completion.info.dt"),
	file = paste(run.dir, "/result_summaries/sage_germline_run_completeion_info.Rdata", sep=""))
