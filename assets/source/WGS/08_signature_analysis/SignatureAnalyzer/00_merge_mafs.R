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
destination.dir <- paste(run.dir, "/result/SignatureAnalyzer", sep="")

if ( length(list.dirs(destination.dir, recursive=FALSE))==0 ){
	system(paste("mkdir -p ", destination.dir, sep=""))
}

load( file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))

tumor.normal.pair.info <- tumor_normal_pair_info.GPL

load(file = paste(run.dir, "/result_summaries/", gpl_prefix, "/candidate.driver.and.variant.info.Rdata", sep=""))
load( file = paste(run.dir, "/result_summaries/", gpl_prefix, "/driver.catalog.germline.somatic.variant.SV.info.Rdata", sep=""))

purple.qc.dt.combined[,tumor.depth:=Purity*AmberMeanDepth]
purple.qc.dt.combined <- purple.qc.dt.combined[grep('FAIL_NO_TUMOR', QCStatus, invert=TRUE),]
purple.qc.dt.combined <- purple.qc.dt.combined[grep('FAIL_CONTAMINATION', QCStatus, invert=TRUE),]
nrow(purple.qc.dt.combined)
##	35

if (remove.copy.number.noise.samples==TRUE){
	purple.qc.dt.combined <- purple.qc.dt.combined[grep('WARN_HIGH_COPY_NUMBER_NOISE', QCStatus, invert=TRUE),]
}

purple.qc.dt.combined <- purple.qc.dt.combined[QCStatus=='PASS' | (QCStatus=='WARN_LOW_PURITY' & tumor.depth>=tumor.depth.cutoff),]
purple.qc.dt.combined <- purple.qc.dt.combined[!tumor.sample.id%in%extra.samples.to.exclude,]

nrow(purple.qc.dt.combined)
##	33

tumor.sample.ids.ls <- purple.qc.dt.combined$tumor.sample.id
names(tumor.sample.ids.ls) <- tumor.sample.ids.ls
sample.ct <- length(tumor.sample.ids.ls)

threads <- min( c(18, sample.ct) )

message( paste("Copying somatic variant maf files for ", study.name, sep=""))

maf.dt.ls <- mclapply(
	X = tumor.sample.ids.ls,
	FUN = function( curr.sample.id ){
		curr.subject.id <- tumor.normal.pair.info[sample.id.T==curr.sample.id,]$subject.id
		
		message(paste("Copying somatic maf for ", curr.sample.id, sep=""))
		curr.maf <- file.path( run.dir, "result", gpl_prefix, curr.subject.id, "purple", curr.sample.id, paste(curr.sample.id, ".purple.somatic.vep.maf", sep=""))
		
		fread(file = curr.maf)
	},
	mc.cores = threads)

maf.dt <- rbindlist(maf.dt.ls)


merged.maf.file <- file.path( destination.dir, "merged.maf")

fwrite( data.table(version="#version 2.4"), merged.maf.file, sep="\t", append=FALSE, col.names = FALSE)	
fwrite( maf.dt, merged.maf.file, sep="\t", append=TRUE, col.names = TRUE)

message( paste("Copied PURPLE processed somatic variant data for ", study.name, sep=""))
