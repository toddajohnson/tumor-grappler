#!/usr/bin/env Rscript

# module use /usr/local/package/modulefiles/
# module load R/4.0.2

options(width=350)
working_dir <- getwd()
gpl_prefix <- basename(working_dir)

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "~/workspace/runs", study.dir )

source( file.path( run.dir, "config_files/common_config.R" ) )

suppressPackageStartupMessages({
	library(data.table)
	library(parallel)})


load( file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))

load( file = paste( run.dir, '/result_summaries/fastq_sorted_bam_statistics.Rdata', sep=""))

## try multiqc on bam.metrics files
processing.stats.dt[,consistency.status:=ifelse(reads_mapped_and_paired.propn.fastq>0.95, 1, 0)]

sample_consistency_summary <- processing.stats.dt[, list(orig.fastq.total=sum(total.read.ct),
		sorted.bam.paired.read.total=sum(reads_paired),
		sorted.bam.mapped.paired.read.total=sum(reads_mapped_and_paired),
		sorted.bam.paired.read.total.frac.fastq=mean(reads_mapped_and_paired.propn.fastq),
		sorted.bam.consistency.status=min(consistency.status)),
	by=list(subject.index, sample.index, subject.id, sample.id)]

sample.info.dt <- sample.info.dt[,
	list(Sample_Name=paste(Sample_Name, collapse=";"), file.ct=sum(file.ct),
		total.file.size=sum(total.file.size), total.file.size.GB=sum(total.file.size.GB)),
	by=list(subject.id, sample.id, subject.index, sample.type, sample.class, sample.class.short)]
#sample.info.dt <- sample.info.dt[order(subject.index, sample.class.short),]
sample.info.dt <- sample.info.dt[order(subject.index, sample.id),]
sample.info.dt[,sample.index:=as.integer(1:nrow(sample.info.dt))]

bam_dir <- paste( run.dir, "/result/bam_markdup", sep="")
bam_final_dir <- paste( run.dir, "/result/bam_final", sep="")

sample.ids.ls <- unique(sample.info.dt$sample.id)
names(sample.ids.ls) <- sample.ids.ls

final_bam_summary.dt.ls <- lapply(
	X = sample.ids.ls,
	FUN = function( curr.sample.id ){
		print(curr.sample.id)
		curr.subject.id <- unique(sample.info.dt[sample.id==curr.sample.id,]$subject.id);
		curr.subject.index <- unique(sample.info.dt[sample.id==curr.sample.id,]$subject.index);
		curr.sample.index <- sample.info.dt[sample.id==curr.sample.id,]$sample.index
		
		get_metrics <- function(curr.file){
			curr.data <- scan(curr.file, what="character", sep="\n")
			hist.row <- which(curr.data=="## HISTOGRAM")
			curr.dt <- fread(metrics_file, skip=3, nrow=hist.row-(3+1))
			orig.cols <- copy(colnames(curr.dt))
			curr.dt[,subject.index:=curr.subject.index]
			curr.dt[,sample.index:=curr.sample.index]
			curr.dt[,subject.id:=curr.subject.id]
			curr.dt[,sample.id:=curr.sample.id]
			curr.dt[,c('subject.index', 'sample.index', 'subject.id', 'sample.id', orig.cols), with=FALSE]
		}
		
		md5_file <- paste( bam_dir, "/", curr.sample.id, ".markdup.bam.md5", sep="")
		metrics_file <- paste( bam_dir, "/", curr.sample.id, ".markdup.metrics", sep="")
		
		if (!file.exists(metrics_file) ){
			md5_file <- paste( bam_final_dir, "/", curr.sample.id, ".bam.md5", sep="")
			metrics_file <- paste( bam_final_dir, "/", curr.sample.id, ".metrics", sep="")
			
			if (!file.exists(metrics_file)){
				curr.dt <- data.table(
					subject.index = curr.subject.index,
					sample.index = curr.sample.index,
					subject.id=curr.subject.id,
					sample.id=curr.sample.id,
					LIBRARY = "MISSING",
					UNPAIRED_READS_EXAMINED = -99,
					READ_PAIRS_EXAMINED = -99,
					UNMAPPED_READS = -99,
					UNPAIRED_READ_DUPLICATES = -99,
					READ_PAIR_DUPLICATES = -99,
					READ_PAIR_OPTICAL_DUPLICATES = -99,
					PERCENT_DUPLICATION = -99,
					ESTIMATED_LIBRARY_SIZE = -99,
					finished = FALSE,
					finalized = -99)
			}else{
				curr.dt <- get_metrics(metrics_file)
				curr.dt[,finished:=file.exists(md5_file)]
				curr.dt[,finalized:=1L]
			}
		}else{
			curr.dt <- get_metrics(metrics_file)
			curr.dt[,finished:=file.exists(md5_file)]
			curr.dt[,finalized:=0L]
		}
		curr.dt
	})


final_bam_summary.dt <- rbindlist(final_bam_summary.dt.ls)
final_bam_summary.dt <- final_bam_summary.dt[,list(
	UNPAIRED_READS_EXAMINED = sum(UNPAIRED_READS_EXAMINED),
	READ_PAIRS_EXAMINED = sum(READ_PAIRS_EXAMINED),
	UNMAPPED_READS = sum(UNMAPPED_READS),
	UNPAIRED_READ_DUPLICATES = sum(UNPAIRED_READ_DUPLICATES),
	READ_PAIR_DUPLICATES = sum(READ_PAIR_DUPLICATES),
	READ_PAIR_OPTICAL_DUPLICATES = sum(READ_PAIR_OPTICAL_DUPLICATES),
	PERCENT_DUPLICATION = mean(PERCENT_DUPLICATION),
	ESTIMATED_LIBRARY_SIZE = max(ESTIMATED_LIBRARY_SIZE)),
	by=list(subject.id, sample.id, subject.index, sample.index, finished, finalized)]

sample_consistency_summary_merged <- merge(
	x = sample_consistency_summary,
	y = final_bam_summary.dt,
	by = c("subject.index", "sample.index", "subject.id", "sample.id"))

sample_consistency_summary_merged[,PAIRED_READS_EXAMINED:=2*READ_PAIRS_EXAMINED]
sample_consistency_summary_merged[,bam.markdup.sorted.read.propn:=PAIRED_READS_EXAMINED/sorted.bam.paired.read.total]
sample_consistency_summary_merged[,bam.consistency.status:=ifelse(bam.markdup.sorted.read.propn > 0.98, 1, 0)]
sample_consistency_summary_merged[,list(subject.id, sample.id, sorted.bam.paired.read.total, READ_PAIRS_EXAMINED, PAIRED_READS_EXAMINED,
	bam.markdup.sorted.read.propn, bam.consistency.status, finished, finalized)]

message("Printing samples table with inconsistent status")
sample_consistency_summary_merged[bam.consistency.status==0, list(subject.id, sample.id, sorted.bam.paired.read.total,  READ_PAIRS_EXAMINED, PAIRED_READS_EXAMINED, bam.markdup.sorted.read.propn, bam.consistency.status)]

message("Printing conistency status count tabulation")
table(sample_consistency_summary_merged$bam.consistency.status)

message("Printing samples table with finished but not finalized")
sample_consistency_summary_merged[finished==TRUE & finalized==0L,
	list(subject.id, sample.id, sample.index, sorted.bam.paired.read.total,  READ_PAIRS_EXAMINED, PAIRED_READS_EXAMINED, bam.markdup.sorted.read.propn, bam.consistency.status)]

save( list=c("sample_consistency_summary_merged"),
	file = paste(run.dir, "/result_summaries/merged_file_consistency_summary.Rdata", sep=""))

fwrite(sample_consistency_summary_merged[,list(subject.id, sample.id,
	read_pair_ct=as.integer(READ_PAIRS_EXAMINED*(1-PERCENT_DUPLICATION)),
	merged.sample.index=as.integer(1:nrow(sample_consistency_summary_merged)))],
	file = paste(run.dir, "/config_files/sample_final_bam_read_cts.csv", sep=""), sep=",")
