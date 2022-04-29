# module use /usr/local/package/modulefiles/
# module load R/4.0.2

working_dir <- getwd()
gpl_prefix <- basename(working_dir)

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "/home/tjohnson/workspace/runs", study.dir )
bam_dir <- paste(run.dir, "/result/bam_final", sep="")

library(data.table)

load( file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep="") )

## for SAGE germline analysis, tumor bam is set to normal sample and reference is set to supporting tumor samples
## tumor= sample(s) in which variants are ascertained, reference = sample(s) that give support for a variant

sample.info.dt[,final.bam:=paste(bam_dir, "/", sample.id, ".bam", sep="")]

sage.germline.run.info.dt <- merge(
	x = sample.info.dt[sample.class.short=='R', 
		 list(tumor.name=paste(sample.id, collapse=","), tumor.bam=paste(final.bam, collapse=",")),
		 by=list(subject.id, subject.index)],
 	y = sample.info.dt[sample.class.short=='T', 
		list(reference.name=paste(sample.id, collapse=","), reference.bam=paste(final.bam, collapse=",")),
		by=list(subject.id, subject.index)],
	by = c('subject.id', 'subject.index'))

sage.germline.run.info.dt[,tumor.only:="FALSE"]

sage.somatic.run.info.dt <- merge(
	x = sample.info.dt[sample.class.short=='T', 
		list(tumor.name=paste(sample.id, collapse=","), tumor.bam=paste(final.bam, collapse=",")),
		by=list(subject.id, subject.index)],
	y = sample.info.dt[sample.class.short=='R', 
		list(reference.name=paste(sample.id, collapse=","), reference.bam=paste(final.bam, collapse=",")),
		by=list(subject.id, subject.index)],
	by = c('subject.id', 'subject.index'),
	all.x = TRUE)
sage.somatic.run.info.dt[,reference.name:=ifelse(is.na(reference.name), 'none', reference.name)]
sage.somatic.run.info.dt[,reference.bam:=ifelse(is.na(reference.bam), 'none', reference.bam)]
sage.somatic.run.info.dt[,tumor.only:=ifelse(reference.name=="none", "TRUE", "FALSE")]


fwrite(sage.germline.run.info.dt[,list(subject.id, tumor.only, tumor.name, tumor.bam, reference.name, reference.bam)],
	sep = "\t",
	file = paste(run.dir, "/config_files/sage_germline_run_info.tsv", sep=""))

fwrite(sage.somatic.run.info.dt[,list(subject.id, tumor.only, tumor.name, tumor.bam, reference.name, reference.bam)],
	sep = "\t",
	file = paste(run.dir, "/config_files/sage_somatic_run_info.tsv", sep=""))

sample.info.dt[,final.bam:=paste(bam_dir, "/", sample.id, ".bam", sep="")]

not.normal.sample.ids <- c('')

sample.info.dt[,sample.class.reassigned:=ifelse(sample.id%in%not.normal.sample.ids, 'T', sample.class.short)];

joint.analysis.file.info <- merge(
		x = sample.info.dt[sample.class.reassigned=='T', 
				list(tumor.name=paste(unique(sample.id), collapse=","), tumor.bam=paste(unique(final.bam), collapse=" ")),
				by=list(subject.id, subject.index)],
		y = sample.info.dt[sample.class.reassigned=='R', 
				list(reference.name=paste(unique(sample.id), collapse=","), reference.bam=paste(unique(final.bam), collapse=" ")),
				by=list(subject.id, subject.index)],
		by = c('subject.id', 'subject.index'),
		all.x = TRUE)
joint.analysis.file.info <- joint.analysis.file.info[order(subject.index),]

joint.analysis.file.info[,joint.labels:=ifelse(is.na(reference.name), tumor.name, paste(reference.name, tumor.name, sep=','))]
joint.analysis.file.info[,joint.bams:=ifelse(is.na(reference.bam), tumor.bam, paste(reference.bam, tumor.bam, sep=' '))]
joint.analysis.file.info <- joint.analysis.file.info[order(subject.index),]

fwrite(joint.analysis.file.info[!is.na(reference.name),list(
	subject.id, subject.index, sample.id.N=reference.name, joint.labels, joint.bams)],
	sep = "\t",
	file=paste(run.dir, '/config_files/joint_tumor_normal_calling_info.tsv', sep=''))

save( list=c("sage.germline.run.info.dt", "sage.somatic.run.info.dt", "joint.analysis.file.info"),
	file = paste(run.dir, "/config_files/sage_run_info.Rdata", sep=""))




