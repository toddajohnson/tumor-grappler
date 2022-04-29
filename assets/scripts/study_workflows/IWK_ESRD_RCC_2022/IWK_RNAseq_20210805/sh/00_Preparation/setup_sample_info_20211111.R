# module use /usr/local/package/modulefiles/
# module load R/4.0.2


#!/usr/bin/env Rscript

#module use /usr/local/package/modulefiles/
#module load R/4.0.2

working_dir <- getwd()
gpl_prefix <- basename(working_dir)

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "~/workspace/runs", study.dir )

source( file.path( run.dir, "config_files/common_config.R" ) )


options(width=200)
options(datatable.print.nrows=250)
library(data.table)
library(parallel)
#library(readxl)

load(file=file.path(wgs.run.dir, "config_files/fastq_file_info.Rdata"))
rm(list=c("sample.read.fastq.info.dt", "readpair.file.info", "readpair.by.lane.file.info",
	"sample.info.dt", "tumor_normal_pair_info", "sample.read.fastq.info.dt.with.lanes",
	"tumor_normal_pair_info.GPL"))

sample.info.dt <- fread( file.path(run.dir, 'config_files/20210802_IWK_RNA_fastq_list.txt'))

setnames(sample.info.dt, old=c("sample.id.T", "Sample_Name.T", "R1", "R2"), new=c("sample.id", "Sample_Name", "Read1", "Read2"))

sample.info.dt[,sample.type:=strsplit(Sample_Name, split='_', fixed=TRUE)[[1]][[length(strsplit(Sample_Name, split='_', fixed=TRUE)[[1]])]], by=1:nrow(sample.info.dt)]
sample.info.dt[,sample.class:=ifelse(sample.type%in%c("N", "Normal", "Blood"), "Reference", ifelse(sample.type%in%c("T", "Tumor", "Cancer"), "Tumor", "Unknown"))]
sample.info.dt[,sample.class.short:=ifelse(sample.class=='Reference', "R", ifelse(sample.class=='Tumor', "T", "UNK"))]

subject.ids.ls <- unique(sample.info.dt$subject.id)
subject.ids.ls <- subject.ids.ls[order(subject.ids.ls)]
subject.index.ls <- as.integer(1:length(subject.ids.ls))
names(subject.index.ls) <- subject.ids.ls

sample.info.dt[,subject.index:=subject.index.ls[subject.id]]

sample.index.info <- sample.info.dt[,list(file.ct=.N), by=list(subject.id, subject.index, sample.id)][order(subject.index, sample.id),]
sample.index.info[,sample.index:=as.integer(1:nrow(sample.index.info))]

sample.info.dt <- merge(
	x = sample.info.dt,
	y = sample.index.info[,list(subject.id, sample.id, sample.index)],
	by = c('subject.id', 'sample.id'))

sample.info.dt[,Platform:="Hiseq"]

subject.info <- sample.index.info[,list(sample.ct=.N), by=list(subject.id, subject.index)]

save( list=c("sample.info.dt", "subject.info"),
	file = paste(run.dir, "/config_files/fastq_file_intermediate_processing_info.Rdata", sep=""))


split_runs <- function( curr.passed.read, curr.passed.read.number ){
	curr.read.ls <- strsplit(curr.passed.read, split="/", fixed=TRUE)[[1]]
	curr.read.ls.l <- length(curr.read.ls)
	curr.dt <- data.table(
		read.number = curr.passed.read.number,
		run.path = paste(curr.read.ls[1:(curr.read.ls.l-2)], collapse="/"),
		fastq.dir = paste(curr.read.ls[(curr.read.ls.l-1):(curr.read.ls.l-1)], collapse="/"),
		fastq.file = curr.read.ls[[curr.read.ls.l]])
	curr.dt
}

sample.read.fastq.info.dt <- rbind(
	sample.info.dt[,split_runs(Read1, 1L), by=list(subject.index, sample.index, subject.id, sample.id, Sample_Name, sample.type, sample.class, sample.class.short, Platform, Read=Read1)],
	sample.info.dt[,split_runs(Read2, 2L), by=list(subject.index, sample.index, subject.id, sample.id, Sample_Name, sample.type, sample.class, sample.class.short, Platform, Read=Read2)])

sample.read.fastq.info.dt <- sample.read.fastq.info.dt[order(subject.index, sample.index, subject.id, sample.type, sample.class, sample.class.short, Platform, fastq.dir, read.number),]

sample.read.fastq.info.dt[,complete.path:=paste(run.path, "/", fastq.dir, "/", fastq.file, sep="")]

fastq.dir.info <- sample.read.fastq.info.dt[,list(file.ct=.N), by=list(run.path, fastq.dir, full.fastq.dir=file.path(run.path, fastq.dir))]

save( list=c("sample.info.dt", "subject.info", "sample.read.fastq.info.dt", "fastq.dir.info"),
		file = paste(run.dir, "/config_files/fastq_file_intermediate_processing_info.Rdata", sep=""))


get_dir_dsmls_status <- function( curr.dir ){
	curr.status <- system(paste("dsmls ", curr.dir, sep=""), intern=TRUE)
	curr.dir.status.line <- grep(curr.dir, curr.status)
	curr.status.tail.line.ct <- length(curr.status) - curr.dir.status.line

	fread( paste(curr.status[(curr.dir.status.line+1):length(curr.status)], collapse='\n'),
		header=FALSE,
		colClasses = c("double", "double", "double", "character", "character"),
		col.names=c('ActS', 'ResS', 'ResB', 'Fst', 'Fname'))
}

get_file_dsmls_status <- function( curr.file ){
	fread( cmd=paste("dsmls ", curr.file, " | tail -n 1", sep=""),
			header=FALSE,
			colClasses = c("double", "double", "double", "character", "character"),
			col.names=c('ActS', 'ResS', 'ResB', 'Fst', 'Fname'))
}

read.dsmls.status <- fastq.dir.info[,get_dir_dsmls_status(full.fastq.dir), by=list(run.path, fastq.dir)]
read.dsmls.status.on.tape <- read.dsmls.status[Fst=='m']
#read.dsmls.status.on.tape
#Empty data.table (0 rows and 7 cols): run.path,fastq.dir,ActS,ResS,ResB,Fst...

recall_file <- function( curr.file ){
	system( paste("qrecall ", curr.file, sep="") )
}

read.dsmls.status.on.tape[,recall_file(complete.path), by=1:nrow(read.dsmls.status.on.tape)]

#save( list=c("sample.read.fastq.info.dt", "read.dsmls.status", "read.dsmls.status.on.tape"),
#	file = paste(run.dir, "/config_files/fastq_file_dsmls_info.Rdata", sep=""))


extract_header_info <- function( curr.file, curr.compress.type ){
	print( paste('Processing ', curr.file ,sep=''))
	
	curr.file.size <- file.info(curr.file)$size
	curr.file.size.GB <- curr.file.size/(1e9)
	if (curr.compress.type=='gz'){
		curr.header.line <- system( paste( "gunzip -c ", curr.file, " | head -n 1", sep=""), intern=TRUE)
	}else{
		curr.header.line <- system( paste( "bunzip2 -c ", curr.file, " | head -n 1", sep=""), intern=TRUE)
	}

	curr.header.line.split <- strsplit(curr.header.line, split=":", fixed=TRUE)[[1]]
	data.table(
		FCID = curr.header.line.split[[3]],
		Lane = curr.header.line.split[[4]],
		file.size = curr.file.size,
		file.size.GB = curr.file.size.GB)
}

sample.read.fastq.info.dt[,compress_type:=ifelse(substring(fastq.file, nchar(fastq.file)-2, nchar(fastq.file))==".gz", "gz", "bz2")]
sample.read.fastq.info.dt[,sample.origin:="IWK_RCC"]

sample.read.fastq.info.dt.with.lanes <- sample.read.fastq.info.dt[,extract_header_info(curr.file=complete.path, curr.compress.type=compress_type),
	by=list(subject.index, sample.index, subject.id, sample.id, Sample_Name, sample.type, sample.class, sample.class.short, sample.origin, Platform, read.number, run.path, fastq.dir, fastq.file, complete.path, compress_type)]

sample.read.fastq.info.dt.with.lanes <- sample.read.fastq.info.dt.with.lanes[,list(
	subject.id, sample.id, Sample_Name, sample.type, sample.class, sample.class.short, sample.origin, Platform, read.number,  FCID, Lane, compress_type, run.path, fastq.dir, fastq.file, complete.path, file.size, subject.index, sample.index)]

sample.read.fastq.info.dt.with.lanes[,lane.ct.per.sample:=.N, by=list(subject.id, sample.id, FCID, read.number)]



## Some samples have multiple barcodes in the same lane
sample.read.fastq.info.dt.with.lanes[,fastq.files.key:=gsub(pattern='_1.fastq', replacement='_X.fastq', x=fastq.file)]
sample.read.fastq.info.dt.with.lanes[,fastq.files.key:=gsub(pattern='_2.fastq', replacement='_X.fastq', x=fastq.files.key)]

sample.read.fastq.info.dt.with.lanes <- sample.read.fastq.info.dt.with.lanes[order(subject.index, sample.index, subject.id, FCID, Lane, fastq.files.key, read.number),]

sample.info.dt <- sample.read.fastq.info.dt.with.lanes[,
		list(file.ct=.N, total.file.size=sum(file.size)),
		by = list(subject.id, subject.index, sample.id, Sample_Name, sample.type, sample.class, sample.class.short, sample.origin)]
sample.info.dt[,total.file.size.GB:=total.file.size/(1e9)]


readpair.file.info <- merge(
	x = unique(sample.read.fastq.info.dt.with.lanes[read.number==1, list(lanes=paste(Lane,collapse=";")),
		by = list(subject.id, sample.id, Sample_Name, FCID, sample.class.short, compress_type, Platform, lane.ct.per.sample, run.path, fastq.dir, fastq.files.key, R1_FASTQ_FILE=fastq.file)]),
	y = unique(sample.read.fastq.info.dt.with.lanes[read.number==2, list(lanes=paste(Lane,collapse=";")),
		by = list(subject.id, sample.id, Sample_Name, FCID, sample.class.short, compress_type, Platform, lane.ct.per.sample, run.path, fastq.dir, fastq.files.key, R2_FASTQ_FILE=fastq.file)]),
	by = c('subject.id', 'sample.id', 'Sample_Name', 'FCID', 'sample.class.short', 'compress_type', 'Platform', 'lane.ct.per.sample', 'lanes', 'run.path', 'fastq.dir', 'fastq.files.key'))

readpair.file.info[,sample.readpair.ct:=.N, by=list(sample.id)]
readpair.file.info[,sample.readpair.index:=as.integer(1:unique(sample.readpair.ct)), by=list(sample.id)]

fwrite(readpair.file.info[,list(subject.id, sample.id, Sample_Name, FCID, sample.class.short, compress_type, Platform,
	lane.ct.per.sample, lanes, run.path, fastq.dir, R1_FASTQ_FILE, R2_FASTQ_FILE, sample.readpair.index)],
	file = paste(run.dir, "/config_files/fastq_readpair_info.csv", sep=""))


readpair.by.lane.file.info <- merge(
	x = sample.read.fastq.info.dt.with.lanes[read.number==1,
		list(subject.id, sample.id, Sample_Name, FCID, sample.class.short, compress_type, Platform, Lane, run.path, fastq.dir, fastq.files.key, R1_FASTQ_FILE=fastq.file)],
	y = sample.read.fastq.info.dt.with.lanes[read.number==2,
		list(subject.id, sample.id, Sample_Name, FCID, sample.class.short, compress_type, Platform, Lane, run.path, fastq.dir, fastq.files.key, R2_FASTQ_FILE=fastq.file)],
	by = c('subject.id', 'sample.id', 'Sample_Name', 'FCID', 'sample.class.short', 'compress_type', 'Platform', 'Lane', 'run.path', 'fastq.dir', 'fastq.files.key'))

readpair.by.lane.file.info[,sample.FCID.Lane.ct:=.N, by=list(subject.id, sample.id, FCID, Lane)]
readpair.by.lane.file.info[,sample.FCID.Lane.idx:=as.integer(1:unique(sample.FCID.Lane.ct)), by=list(subject.id, sample.id, FCID, Lane)]

readpair.by.lane.file.info[,sample.readpair.ct:=.N, by=list(sample.id)]
readpair.by.lane.file.info[,sample.readpair.index:=as.integer(1:unique(sample.readpair.ct)), by=list(sample.id)]


fwrite(readpair.by.lane.file.info[,list(subject.id, sample.id, Sample_Name, FCID, sample.class.short, compress_type, Platform,
	Lane, run.path, fastq.dir, R1_FASTQ_FILE, R2_FASTQ_FILE, sample.FCID.Lane.idx, sample.readpair.index)],
	file = paste(run.dir, "/config_files/fastq_readpair_by_lane_info.csv", sep=""))


fwrite(subject.info[,list(subject.id, subject.index)],
		file = paste(run.dir, "/config_files/subject_info.csv", sep=""))

fwrite(sample.info.dt[,list(subject.id, sample.id, Sample_Name, sample.class.short)],
	file = paste(run.dir, "/config_files/sample_info.csv", sep=""))


##sample.info.dt[,WGS.sample.id:=paste(subject.id, '_C', sep="")]
sample.info.dt[,WGS.sample.id:=paste(subject.id, '_T', sep="")]
## 20210818 added WGS.sample.id column; need to unify sample IDs
fwrite(sample.info.dt[,list(subject.id, sample.id, Sample_Name, sample.class.short, WGS.sample.id)],
		file = paste(run.dir, "/config_files/sample_info.csv", sep=""))

fwrite(sample.info.dt[,list(WGS.sample.id)],
	file = paste(run.dir, "/config_files/Isofox_sample_list.csv", sep=""))

sample.info.with.clinical.info.dt <- merge(
	x = sample.info.dt,
	y = patient.info.dt[,list(WGS.sample.id=tumor.sample.id, CancerType=histology, CohortName='ESRD RCC', age, gender, years.dialysis)],
	by = 'WGS.sample.id')

fwrite(sample.info.with.clinical.info.dt[,list(WGS.sample.id, CancerType, CohortName)],
	file = paste(run.dir, "/config_files/Isofox_sample_cohort_info.csv", sep=""))


readpair.by.lane.file.info[,ReadGroup.ID:=paste("ID:", FCID, ".", Sample_Name, ".", Lane, sep="")]
readpair.by.lane.file.info[,ReadGroup.string:=paste(ReadGroup.ID, "\tSM:", sample.id, "\tLB:", Sample_Name, "\tPU:ILLUMINA", "\tPL:", ReadGroup.ID, sep="")]

#readpair.by.lane.file.info[,fastq.r1:=file.path(run.dir, 'fastq', subject.id, paste(sample.id, '_', FCID, '_L00', Lane, '_R1.trimmed.fastq', sep=''))]
#readpair.by.lane.file.info[,fastq.r2:=file.path(run.dir, 'fastq', subject.id, paste(sample.id, '_', FCID, '_L00', Lane, '_R2.trimmed.fastq', sep=''))]

readpair.by.lane.file.info[,fastq.r1:=file.path(run.dir, 'fastq', subject.id, paste(sample.id, '_', FCID, '_L00', Lane, '_R1.trimmed.cleaned.fastq', sep=''))]
readpair.by.lane.file.info[,fastq.r2:=file.path(run.dir, 'fastq', subject.id, paste(sample.id, '_', FCID, '_L00', Lane, '_R2.trimmed.cleaned.fastq', sep=''))]

fastq.file.info.for.STAR <- readpair.by.lane.file.info[,list(fastq.files.r1=paste(fastq.r1, collapse=','), fastq.files.r2=paste(fastq.r2, collapse=','), ReadGroups=paste(ReadGroup.string, collapse=' , ')), by=list(subject.id, sample.id)][,list(readFilesIn=paste('--readFilesIn ', fastq.files.r1, ' ', fastq.files.r2, sep=''), outSAMattrRGline=paste('--outSAMattrRGline ', ReadGroups, sep='')), by=list(subject.id, sample.id)]

# switched to "|" separator as last sample separates files within a field by commas
fwrite(fastq.file.info.for.STAR, file=paste(run.dir, "/config_files/STAR_fastq_file_info.txt", sep=""), sep="|", col.names=TRUE)

save( list=c("sample.read.fastq.info.dt", "readpair.file.info", "readpair.by.lane.file.info",
	"sample.info.dt",  "sample.read.fastq.info.dt.with.lanes", "sample.info.with.clinical.info.dt", "fastq.file.info.for.STAR"),
	file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))
