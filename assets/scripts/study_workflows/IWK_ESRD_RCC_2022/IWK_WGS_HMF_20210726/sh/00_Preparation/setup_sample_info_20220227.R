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
library(readxl)


config.info.dt <- data.table(
	target.dir = c(
			"/home/tjohnson/icgc_cgm/premium_nearline/share/03_result/run_MGISEQ_IWK_WGS_Genomon2_20190213",
			"/home/tjohnson/icgc_cgm/premium_nearline/share/03_result/run_MGISEQ_IWK_WGS_Genomon2_20190408",
			"/home/tjohnson/icgc_cgm/premium_nearline/share/03_result/run_MGISEQ_IWK_WGS_Genomon2_20200811",
			"/home/tjohnson/icgc_cgm/premium_nearline/share/05_result/run_MGISEQ_IWK_WGS_Genomon2_20190716/result_bam-1",
			"/home/tjohnson/icgc_cgm/premium_nearline/share/05_result/run_MGISEQ_IWK_WGS_Genomon2_20190716/result_bam-2",
			"/home/tjohnson/icgc_cgm/premium_nearline/share/05_result/run_MGISEQ_IWK_WGS_Genomon2_20190716/result_bam-3",
			"/home/tjohnson/icgc_cgm/premium_nearline/share/05_result/run_MGISEQ_IWK_WGS_Genomon2_20190716/result_bam-4"),
	config.file =c(
			"sample_config//MGISEQ_IWK_WGS_sample_config.csv",
			"sample_config/MGISEQ_IWK_WGS_sample_config.csv",
			"sample_config/MGISEQ_IWK_WGS_sample_config.csv",
			"config/190717_MGI_bam-lonly-1_20190717_092036_172913.csv",
			"config/190717_MGI_bam-lonly-2_20190718_023626_494414.csv",
			"config/190717_MGI_bam-lonly-3_20190719_123230_688154.csv",
			"config/190717_MGI_bam-lonly-4_20190720_170603_634869.csv"),
	fastq.line.ct.dir = c(
			"result/fastq",
			"result/fastq",
			"result/fastq",
			"fastq",
			"fastq",
			"fastq",
			"fastq"))

## need to reparse config.file; (row.past.fastq.info-2) does not work in one file
config.file.info.dt <- rbindlist(apply(
	X = config.info.dt,
	MAR = 1,
	FUN = function( curr.file.info ){
		curr.dir <- curr.file.info[['target.dir']]
		curr.config.file <- curr.file.info[['config.file']]
		curr.fastq.line.ct.dir <- curr.file.info[['fastq.line.ct.dir']]
		
		curr.config.file <- paste(curr.dir, '/', curr.config.file, sep='')
		curr.config.info <- scan(file=curr.config.file, what="character")
		row.past.fastq.info <- which(curr.config.info =='[bam_import]')
		
		if( length(row.past.fastq.info)>0 ){
			curr.config.info <- curr.config.info[2:(row.past.fastq.info-1)];
		}else{
			curr.config.info <- curr.config.info[2:length(curr.config.info)];
		}
		rbindlist(
			lapply(
				X = curr.config.info,
				FUN = function( curr.string ){
					if (curr.string!=''){
						sample.files.R1.string <- strsplit(x=curr.string, split=',')[[1]][[2]]
						sample.files.R2.string <- strsplit(x=curr.string, split=',')[[1]][[3]]
					
						curr.sample.id <- strsplit(x=curr.string, split=',')[[1]][[1]]
					
						data.table(
							Sample_Name = curr.sample.id,
							R1.file = strsplit(x=sample.files.R1.string, split=';')[[1]],
							R2.file = strsplit(x=sample.files.R2.string, split=';')[[1]],
							fastq.line.ct.file = paste(curr.dir, '/', curr.fastq.line.ct.dir, '/', curr.sample.id, '/fastq_line_num.txt', sep=''),
							config.file = curr.config.file)
					}else{
						NULL
					}
				}))
	}))


sample.info.dt <- config.file.info.dt[,list(Sample_Name, Read1=R1.file, Read2=R2.file, config.file, fastq.line.ct.file)]


sample.info.dt[,subject.id:=strsplit(Sample_Name, split='_', fixed=TRUE)[[1]][[1]], by=1:nrow(sample.info.dt)]
sample.info.dt[,subject.id:=ifelse(subject.id=='IWK031', 'IWK017', subject.id)]

sample.info.dt[,sample.id:=strsplit(Sample_Name, split='_', fixed=TRUE)[[1]][[1]], by=1:nrow(sample.info.dt)]

sample.info.dt[,sample.type:=strsplit(Sample_Name, split='_', fixed=TRUE)[[1]][[length(strsplit(Sample_Name, split='_', fixed=TRUE)[[1]])]], by=1:nrow(sample.info.dt)]
sample.info.dt[,sample.class:=ifelse(sample.type%in%c("N", "Normal", "Blood"), "Reference", ifelse(sample.type%in%c("T", "Tumor", "Cancer"), "Tumor", "Unknown"))]
sample.info.dt[,sample.class.short:=ifelse(sample.class=='Reference', "R", ifelse(sample.class=='Tumor', "T", "UNK"))]
sample.info.dt[,sample.id:=paste(sample.id, "_", substring(sample.class.short, 1, 1), sep="")]

subject.ids.ls <- unique(sample.info.dt$subject.id)
subject.ids.ls <- subject.ids.ls[order(subject.ids.ls)]
subject.index.ls <- as.integer(1:length(subject.ids.ls))
names(subject.index.ls) <- subject.ids.ls

sample.info.dt[,subject.index:=subject.index.ls[subject.id]]

sample.index.info <- sample.info.dt[,list(file.ct=.N), by=list(subject.id, subject.index, sample.id)][order(subject.index, sample.id),]
sample.index.info[,sample.index:=as.integer(1:nrow(sample.index.info))]

# Needed to changed subject id for IWK031, as newer tumor from IWK017 (same germline sample)
#sample.index.info[,list(sample.ct=.N), by=list(subject.id)][which(sample.ct<2),]
#subject.id sample.ct
#1:     IWK031         1
#sample.info.dt[subject.id%in%sample.index.info[,list(sample.ct=.N), by=list(subject.id)][which(sample.ct<2),]$subject.id,]

sample.info.dt <- merge(
	x = sample.info.dt,
	y = sample.index.info[,list(subject.id, sample.id, sample.index)],
	by = c('subject.id', 'sample.id'))

sample.info.dt[,Platform:="MGISEQ"]

subject.info <- sample.index.info[,list(sample.ct=.N), by=list(subject.id, subject.index)]

save( list=c("sample.info.dt", "subject.info"),
	file = paste(run.dir, "/config_files/fastq_file_intermediate_processing_info.Rdata", sep=""))


sample.merged.fastq.info.dt <- sample.info.dt[,list(file.ct=.N), by=list(subject.id, subject.index, sample.id, sample.index, Sample_Name, fastq.line.ct.file)]
sample.merged.fastq.info.dt[,fastq.line.ct:=scan(fastq.line.ct.file), by=1:nrow(sample.merged.fastq.info.dt)]

test <- copy(sample.info.dt)
#test[,Read1.split.ct:=length(strsplit(Read1, split='/', fixed=TRUE)[[1]]), by=1:nrow(test)]
#table(test$Read1.split.ct)
##	11  12  13  14 
##	672 263 210  84 

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


extract_header_info <- function( curr.file ){
	print( paste('Processing ', curr.file ,sep=''))
	
	curr.file.size <- file.info(curr.file)$size
	curr.file.size.GB <- curr.file.size/(1e9)
	curr.header.line <- system( paste( "gunzip -c ", curr.file, " | head -n 1", sep=""), intern=TRUE)
	curr.header.line.split <- strsplit(curr.header.line, split=":", fixed=TRUE)[[1]]
	data.table(
		FCID = curr.header.line.split[[3]],
		Lane = curr.header.line.split[[4]],
		file.size = curr.file.size,
		file.size.GB = curr.file.size.GB)
}

sample.read.fastq.info.dt.with.lanes <- sample.read.fastq.info.dt[,extract_header_info(curr.file=complete.path),
	by=list(subject.index, sample.index, subject.id, sample.id, Sample_Name, sample.type, sample.class, sample.class.short, Platform, read.number, run.path, fastq.dir, fastq.file, complete.path)]

sample.read.fastq.info.dt.with.lanes <- sample.read.fastq.info.dt.with.lanes[,list(
	subject.id, sample.id, Sample_Name, sample.type, sample.class, sample.class.short, Platform, read.number,  FCID, Lane, run.path, fastq.dir, fastq.file, complete.path, file.size, subject.index, sample.index)]

sample.read.fastq.info.dt.with.lanes[,lane.ct.per.sample:=.N, by=list(subject.id, sample.id, FCID, read.number)]

sample.read.fastq.info.dt.with.lanes[,compress_type:=ifelse(substring(fastq.file, nchar(fastq.file)-2, nchar(fastq.file))==".gz", "gz", "bz2")]
sample.read.fastq.info.dt.with.lanes[,sample.origin:="IWK_RCC"]

## Some samples have multiple barcodes in the same lane
sample.read.fastq.info.dt.with.lanes[,fastq.files.key:=gsub(pattern='_R1_', replacement='_RX_', x=fastq.file)]
sample.read.fastq.info.dt.with.lanes[,fastq.files.key:=gsub(pattern='_R2_', replacement='_RX_', x=fastq.files.key)]

sample.read.fastq.info.dt.with.lanes <- sample.read.fastq.info.dt.with.lanes[order(subject.index, sample.index, subject.id, FCID, Lane, fastq.files.key, read.number),]

sample.info.dt <- sample.read.fastq.info.dt.with.lanes[,
		list(file.ct=.N, total.file.size=sum(file.size)),
		by = list(subject.id, subject.index, sample.id, Sample_Name, sample.type, sample.class, sample.class.short, sample.origin)]
sample.info.dt[,total.file.size.GB:=total.file.size/(1e9)]

tumor.info <- sample.read.fastq.info.dt.with.lanes[sample.class.short=='T',
	list(sample.ct=.N),
	by=list(subject.id, subject.index, sample.id.T=sample.id, Sample_Name.T=Sample_Name, sample.origin)]
tumor.info[,tumor.index:=as.integer(1:nrow(tumor.info))]

tumor_normal_pair_info <- merge(
	x = tumor.info[,list(subject.id, subject.index, tumor.index, sample.id.T, Sample_Name.T, sample.origin)],
	y = unique(sample.read.fastq.info.dt.with.lanes[sample.class.short%in%c('N', 'R'), list(subject.id, subject.index, sample.id.N=sample.id, Sample_Name.N=Sample_Name)]),
	by = c('subject.id', 'subject.index'),
	all.x = TRUE)


tumor_normal_pair_info[,tumor.status:=ifelse(is.na(sample.id.T), 0, 1)]
tumor_normal_pair_info[,normal.status:=ifelse(is.na(sample.id.N), 0, 1)]
tumor_normal_pair_info[,sample.id.N:=ifelse(is.na(sample.id.N), "none", sample.id.N)]
tumor_normal_pair_info[,Sample_Name.N:=ifelse(is.na(Sample_Name.N), "none", Sample_Name.N)]

tumor_normal_pair_info <- tumor_normal_pair_info[,list(Sample_Name.T=paste(Sample_Name.T, collapse=";"), Sample_Name.N=paste(Sample_Name.N, collapse=";"), Sample_Name.T.ct=sum(tumor.status), Sample_Name.N.ct=sum(normal.status)),
		by=list(subject.id, subject.index, tumor.index, sample.id.T, sample.origin, sample.id.N, normal.status)]

tumor_normal_pair_info <- tumor_normal_pair_info[order(subject.index, tumor.index),]

tumor_normal_pair_info[,sample.id.N.ct:=sum(normal.status), by=list(subject.id, sample.id.N)]

tumor_normal_pair_info[,sample.id.N.index:=ifelse(sample.id.N.ct==0, 0L, as.integer(1:unique(sample.id.N.ct))), by=list(subject.id, sample.id.N)]

## sample.id.N.index needed for sample set that had multiple tumor to single normal sample mapping
fwrite(tumor_normal_pair_info[,list(subject.id, subject.index, tumor.index, sample.id.T, Sample_Name.T, sample.id.N, Sample_Name.N, Sample_Name.T.ct, Sample_Name.N.ct, sample.id.N.index)],
	file = paste(run.dir, "/config_files/tumor_normal_pair_info.csv", sep=""))


readpair.file.info <- merge(
	x = unique(sample.read.fastq.info.dt.with.lanes[read.number==1, list(lanes=paste(Lane,collapse=";")),
		by = list(subject.id, sample.id, Sample_Name, FCID, sample.class.short, compress_type, Platform, lane.ct.per.sample, run.path, fastq.dir, fastq.files.key, R1_FASTQ_FILE=fastq.file)]),
	y = unique(sample.read.fastq.info.dt.with.lanes[read.number==2, list(lanes=paste(Lane,collapse=";")),
		by = list(subject.id, sample.id, Sample_Name, FCID, sample.class.short, compress_type, Platform, lane.ct.per.sample, run.path, fastq.dir, fastq.files.key, R2_FASTQ_FILE=fastq.file)]),
	by = c('subject.id', 'sample.id', 'Sample_Name', 'FCID', 'sample.class.short', 'compress_type', 'Platform', 'lane.ct.per.sample', 'lanes', 'run.path', 'fastq.dir', 'fastq.files.key'))

fwrite(readpair.file.info[,list(subject.id, sample.id, Sample_Name, FCID, sample.class.short, compress_type, Platform,
	lane.ct.per.sample, lanes, run.path, fastq.dir, R1_FASTQ_FILE, R2_FASTQ_FILE)],
	file = paste(run.dir, "/config_files/fastq_readpair_info.csv", sep=""))


readpair.by.lane.file.info <- merge(
	x = sample.read.fastq.info.dt.with.lanes[read.number==1,
		list(subject.id, sample.id, Sample_Name, FCID, sample.class.short, compress_type, Platform, Lane, run.path, fastq.dir, fastq.files.key, R1_FASTQ_FILE=fastq.file)],
	y = sample.read.fastq.info.dt.with.lanes[read.number==2,
		list(subject.id, sample.id, Sample_Name, FCID, sample.class.short, compress_type, Platform, Lane, run.path, fastq.dir, fastq.files.key, R2_FASTQ_FILE=fastq.file)],
	by = c('subject.id', 'sample.id', 'Sample_Name', 'FCID', 'sample.class.short', 'compress_type', 'Platform', 'Lane', 'run.path', 'fastq.dir', 'fastq.files.key'))

readpair.by.lane.file.info[,sample.FCID.Lane.ct:=.N, by=list(subject.id, sample.id, FCID, Lane)]
readpair.by.lane.file.info[,sample.FCID.Lane.idx:=as.integer(1:unique(sample.FCID.Lane.ct)), by=list(subject.id, sample.id, FCID, Lane)]

fwrite(readpair.by.lane.file.info[,list(subject.id, sample.id, Sample_Name, FCID, sample.class.short, compress_type, Platform,
	Lane, run.path, fastq.dir, R1_FASTQ_FILE, R2_FASTQ_FILE, sample.FCID.Lane.idx)],
	file = paste(run.dir, "/config_files/fastq_readpair_by_lane_info.csv", sep=""))


fwrite(subject.info[,list(subject.id, subject.index)],
		file = paste(run.dir, "/config_files/subject_info.csv", sep=""))

save( list=c("sample.read.fastq.info.dt", "readpair.file.info", "readpair.by.lane.file.info",
	"sample.info.dt", "tumor_normal_pair_info", "sample.read.fastq.info.dt.with.lanes"),
	file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))



load(file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))


## Modify and add below depending on sample QC: removed samples that normal has tumor-like CNA
## from QDNAseq CN profiling
load(file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))

fwrite(sample.info.dt[,list(subject.id, sample.id, sample.class.short, bam.dir=file.path(run.dir, 'result/bam_final'), bam.file=paste(sample.id, '.bam', sep=''))],
		file = paste(run.dir, "/config_files/bam_info.csv", sep=""))

tumor.sample.indexes.to.include <- c(1:38)

tumor_normal_pair_info.GPL <- tumor_normal_pair_info[tumor.sample.indexes.to.include,]

fwrite(tumor_normal_pair_info.GPL,
		file = paste(run.dir, "/config_files/tumor_normal_pair_info_GPL.csv", sep=""))

#tumor.diseases <- unique(tumor_normal_pair_info.GPL$sample.origin)
## tumor.diseases
## "lung adenocarcinoma"                 "lung squamous carcinoma"             "small cell lung cancer"              "large cell neuroendocrine carcinoma"
t#umor.disease.abbrev <- c('LUAD', 'LUSC', 'SCLC', 'LCNEC')
#names(tumor.disease.abbrev) <- tumor.diseases

## added 5/27 to allow grouping of samples by disease/tumor class and not just study name
#tumor_normal_pair_info.GPL[,Disease:=tumor.disease.abbrev[sample.origin]]
#tumor_normal_pair_info.GPL[,tumor.class:="Organoid"]

fwrite(tumor_normal_pair_info.GPL[,list(subject.id, tumor.sample.id=sample.id.T, normal.sample.id=sample.id.N)],
	file = file.path(run.dir, "result_summaries", gpl_prefix, "vep_annotation_sample_info.tsv"), sep='\t')


save( list=c("sample.read.fastq.info.dt", "readpair.file.info", "readpair.by.lane.file.info",
				"sample.info.dt", "tumor_normal_pair_info", "sample.read.fastq.info.dt.with.lanes",
				"tumor_normal_pair_info.GPL"),
		file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))

patient.info.file <- paste(run.dir, "/config_files/DialysisAssociatedKidneyCancer.v3.xlsx", sep="")

patient.info <- read_xlsx(path = patient.info.file, sheet='patient_list')
patient.info.dt <- as.data.table(patient.info)

histology.short.map <- c('ACD-RCC', 'ccRCC', 'pRCC', 'ccpRCC', 'chRCC', 'uncRCC')
names(histology.short.map ) <- c('ACD-RCC', 'Clear cell', 'Papillary', 'Clear cell papillary', 'Chromophobe', 'Unclassified')

patient.info.dt[,histology:=ifelse(histology=='Papillary (Clear papillary)', 'Clear cell papillary', histology)]
patient.info.dt[,histology.short:=histology.short.map[histology]]

save( list=c("sample.read.fastq.info.dt", "readpair.file.info", "readpair.by.lane.file.info",
	"sample.info.dt", "tumor_normal_pair_info", "sample.read.fastq.info.dt.with.lanes",
	"tumor_normal_pair_info.GPL", "patient.info.dt"),
	file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))
