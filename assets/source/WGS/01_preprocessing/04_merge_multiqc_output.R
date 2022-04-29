#!/usr/bin/env Rscript

# module use /usr/local/package/modulefiles/
# module load R/4.0.2

#####
## export PATH="/home/tjohnson/.local/bin:${PATH}"
## Run multiqc . in ../../result/fastqc and ../../result/bam/stats
####

## Copy from pipeline directory 
## Depending on input fastq, I always seem to have to change the Sample Id to link tables below

options(width=350)
options(datatable.print.nrows=2500)

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

dir.create( file.path(run.dir, '/result_summaries'), showWarnings = FALSE, recursive = FALSE, mode = "0777")

load(file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))

fastqc.dir <- paste(run.dir, "/result/fastqc", sep="")
bamstats.dir <- paste(run.dir, "/result/bam/stats", sep="")

message(paste("Reading multiqc output for ", study.name, sep=""))
# ~/.multiqc_config.yaml needs to have max_table_rows set high enough to handle the number of fastq files
# or multiqc_general_stats.txt will not be created
multiqc.fastq.cts.dt <- fread(file=paste(fastqc.dir, "/multiqc_data/multiqc_general_stats.txt", sep=""))
multiqc.sorted.bam.general.stats.dt <- fread(file=paste(bamstats.dir, "/multiqc_data/multiqc_general_stats.txt", sep=""))
multiqc.sorted.bam.samtools.stats.dt <- fread(file=paste(bamstats.dir, "/multiqc_data/multiqc_samtools_stats.txt", sep=""))

#sample.read.fastq.info.dt.with.lanes[,read.index:=as.integer(1:nrow(sample.read.fastq.info.dt.with.lanes))]

sample.read.fastq.info.dt.with.lanes <- merge(
	x = sample.read.fastq.info.dt.with.lanes,
	y = readpair.by.lane.file.info[,list(subject.id, sample.id, FCID, Lane, fastq.files.key, sample.FCID.Lane.idx)],
	by = c("subject.id", "sample.id", "FCID", "Lane", "fastq.files.key"),
	all = TRUE)

extract_sample_from_filename <- function( curr.file ){
	curr.file.ls <- strsplit(curr.file, split='.', fixed=TRUE)[[1]][[1]]
}
#sample.read.fastq.info.dt.with.lanes[,Sample:=substring(fastq.file, 1, nchar(fastq.file)-9)]
#sample.read.fastq.info.dt.with.lanes[,Sample:=strsplit(fastq.file, split='.', fixed=TRUE)[[1]][[1]], by=1:nrow(sample.read.fastq.info.dt.with.lanes)]
sample.read.fastq.info.dt.with.lanes[,Sample:=paste(sample.id, "_", FCID, "_L00", ifelse(sample.FCID.Lane.idx==2, paste(Lane, "_", sample.FCID.Lane.idx, sep=""), Lane), "_R", read.number, sep=""), by=1:nrow(sample.read.fastq.info.dt.with.lanes)]

message(paste("Merging multiqc output for ", study.name, " fastq files", sep=""))


fastq.stats.dt <- merge(
	x = sample.read.fastq.info.dt.with.lanes[,list(
		subject.index, sample.index, subject.id, sample.id, smpl.class=sample.class.short, FCID, Lane, fastq.files.key, sample.FCID.Lane.idx, read.number, Sample)],
	y = multiqc.fastq.cts.dt,
	by = c("Sample"),
	all = TRUE)

setnames(fastq.stats.dt, old=c(
	'FastQC_mqc-generalstats-fastqc-percent_duplicates',
	'FastQC_mqc-generalstats-fastqc-percent_gc',
	'FastQC_mqc-generalstats-fastqc-avg_sequence_length',
	'FastQC_mqc-generalstats-fastqc-percent_fails',
	'FastQC_mqc-generalstats-fastqc-total_sequences'),
	new=c('pct.duplicates', 'pct.gc', 'avg_read_length', 'pct.failed', 'total.read.ct'))

fastq.stats.dt[,propn.duplicates:=pct.duplicates/100]

fastq.stats.dt[,total.read.ct:=as.integer(total.read.ct)]
fastq.stats.dt[,unq.read.ct:=round(total.read.ct*(1-propn.duplicates))]
fastq.stats.dt[,dupl.read.ct:=round(total.read.ct*propn.duplicates)]
fastq.stats.dt[,Sample:=paste(sample.id, '_', FCID, '_L00', ifelse(sample.FCID.Lane.idx==2, paste(Lane, "_", sample.FCID.Lane.idx, sep=""), Lane), sep='')]

message(paste("Merging multiqc output for ", study.name, " sorted bam files", sep=""))

processing.stats.dt <- merge(
	x = fastq.stats.dt[,
		list( unq.read.ct=sum(as.numeric(unq.read.ct)), dupl.read.ct=sum(as.numeric(dupl.read.ct)), total.read.ct=sum(as.numeric(total.read.ct))),
		by=list(Sample, subject.index, sample.index, subject.id, sample.id, smpl.class, FCID, Lane, sample.FCID.Lane.idx)],
	y = multiqc.sorted.bam.samtools.stats.dt[,list(Sample, reads_mapped, reads_mapped_and_paired, reads_paired, supplementary_alignments)],
	by = c("Sample"))


processing.stats.dt[,reads_mapped.propn.fastq:=reads_mapped/total.read.ct]
processing.stats.dt[,reads_mapped_and_paired.propn.fastq:=reads_mapped_and_paired/total.read.ct]

save(processing.stats.dt,
	file = paste( run.dir, '/result_summaries/fastq_sorted_bam_statistics.Rdata', sep=""))

message(paste("Finished processing multiqc output for ", study.name, sep=""))

print( processing.stats.dt[,list(subject.index, sample.index, subject.id, sample.id, smpl.class, FCID, Lane, sample.FCID.Lane.idx, total.read.ct, reads_mapped.propn.fastq, reads_mapped_and_paired.propn.fastq)] )
