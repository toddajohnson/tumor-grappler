#!/usr/bin/env Rscript

# module use /usr/local/package/modulefiles/
#module load R/4.0.2

working_dir <- getwd()

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "~/workspace/runs", study.dir )

## arguments that could be added to common_config.R
remove.copy.number.noise.samples <- TRUE
tumor.depth.cutoff <- 10

source( file.path( run.dir, "config_files/common_config.R" ) )

suppressPackageStartupMessages({
	library(data.table)
	library(parallel)})

purple.dir <- paste(run.dir, "/result/", gpl_prefix, sep="")
destination.dir <- paste(run.dir, "/result/SigProfiler/", study.name, "_mutation_signature_analysis", sep="")

if ( length(list.dirs(destination.dir, recursive=FALSE))==0 ){
	system(paste("mkdir -p ", destination.dir, sep=""))
}

load( file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))

tumor.normal.pair.info <- tumor_normal_pair_info.GPL

if(exists('version.date.for.driver.catalog.file')==TRUE){
	load(file = paste(run.dir, "/result_summaries/", gpl_prefix, "/candidate.driver.and.variant.info_", version.date.for.driver.catalog.file, ".Rdata", sep=""))
	load( file = paste(run.dir, "/result_summaries/", gpl_prefix, "/driver.catalog.germline.somatic.variant.SV.info.", version.date.for.driver.catalog.file, ".Rdata", sep=""))
	
}else{
	load(file = paste(run.dir, "/result_summaries/", gpl_prefix, "/candidate.driver.and.variant.info.Rdata", sep=""))
	load( file = paste(run.dir, "/result_summaries/", gpl_prefix, "/driver.catalog.germline.somatic.variant.SV.info.Rdata", sep=""))
}

starting.sample.ct <- nrow(purple.qc.dt.combined)

purple.qc.dt.combined[,tumor.depth:=Purity*AmberMeanDepth]
purple.qc.dt.combined.QCd <- purple.qc.dt.combined[grep('FAIL_NO_TUMOR', QCStatus, invert=TRUE),]
purple.qc.dt.combined.QCd <- purple.qc.dt.combined.QCd[grep('FAIL_CONTAMINATION', QCStatus, invert=TRUE),]

nonfail.sample.ct <- nrow(purple.qc.dt.combined.QCd)

#modified QC categories 3/11/2022
if (remove.copy.number.noise.samples==TRUE){
	purple.qc.dt.combined.QCd <- purple.qc.dt.combined.QCd[grep('WARN_HIGH_COPY_NUMBER_NOISE', QCStatus, invert=TRUE),]
}

if (remove.deleted.genes.samples==TRUE){
	purple.qc.dt.combined.QCd <- purple.qc.dt.combined.QCd[grep('WARN_DELETED_GENES', QCStatus, invert=TRUE),]
}

if (remove.gender.mismatch.samples==TRUE){
	purple.qc.dt.combined.QCd <- purple.qc.dt.combined.QCd[grep('WARN_GENDER_MISMATCH', QCStatus, invert=TRUE),]
}

purple.qc.dt.combined.QCd[,low.purity:=0L]
purple.qc.dt.combined.QCd[grep('WARN_LOW_PURITY', QCStatus, fixed=TRUE),low.purity:=1L]
purple.qc.dt.combined.QCd <- purple.qc.dt.combined.QCd[which(!(low.purity==1L & tumor.depth<tumor.depth.cutoff)),]

purple.qc.dt.combined.QCd <- purple.qc.dt.combined.QCd[!tumor.sample.id%in%extra.samples.to.exclude,]

QC.pass.samples <- purple.qc.dt.combined.QCd$tumor.sample.id

tumor.sample.ids.ls <- QC.pass.samples
names(tumor.sample.ids.ls) <- tumor.sample.ids.ls
sample.ct <- length(tumor.sample.ids.ls)

print( paste( study.name, ' study sample ct.: ', starting.sample.ct, '/', nonfail.sample.ct, '/', sample.ct, ' (starting/non-fail/QC+)', sep=''))

threads <- min( c(18, sample.ct) )

message( paste("Copying somatic variant vcf files for ", study.name, sep=""))

catch.output.ls <- mclapply(
	X = tumor.sample.ids.ls,
	FUN = function( curr.sample.id ){
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		curr.normal.sample.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$sample.id.N

		message(paste("Copying somatic vcf data for ", curr.sample.id, sep=""))
		file.gz <- paste(purple.dir, "/", curr.subject.id, "/purple/", curr.sample.id, "/", curr.sample.id, ".purple.somatic.vcf.gz", sep="")
		system( paste('gunzip -c ', file.gz, ' > ', destination.dir, "/", curr.sample.id, ".vcf", sep="") )
	},
	mc.cores = threads)

message( paste("Copied PURPLE processed somatic variant data for ", study.name, sep=""))
