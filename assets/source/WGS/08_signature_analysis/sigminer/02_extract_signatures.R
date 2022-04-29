#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

run.dir <- args[1]
sig.type <- args[2]

options(width=250)

working_dir <- getwd()

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "~/workspace/runs", study.dir )

source( file.path( run.dir, "config_files/common_config.R" ) )

library(sigminer)
library(NMF)

purple.dir <- paste(run.dir, "/result/", gpl_prefix, sep="")

load( file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))

tumor.normal.pair.info <- tumor_normal_pair_info.GPL
tumor.sample.ids.ls <- tumor.normal.pair.info$sample.id.T
names(tumor.sample.ids.ls) <- tumor.sample.ids.ls
sample.ct <- length(tumor.sample.ids.ls)

tumor.normal.pair.info[,somatic.vcf:=file.path(purple.dir, subject.id, "purple", sample.id.T, paste(sample.id.T, ".purple.somatic.vcf.gz", sep="") )]
tumor.normal.pair.info[,cnv.file:=file.path(purple.dir, subject.id, "purple", sample.id.T, paste(sample.id.T, ".purple.segment.tsv", sep="") )]

fig.dir <- file.path(run.dir, "result/sigminer/figures")

if ( file.exists(file.path(run.dir, "result/sigminer"), recursive=FALSE)==FALSE ){
	dir.create(file.path(run.dir, "result/sigminer"))
}
if ( file.exists(file.path(run.dir, "result/sigminer/figures"), recursive=FALSE)==FALSE ){
	dir.create(file.path(run.dir, "result/sigminer/figures"))
}

if (sig.type%in%c('SBS', 'DBS', 'ID') ){
	load( file = paste(run.dir, "/result/sigminer/mt_tally_", sig.type, ".Rdata", sep=""))
	curr.nmf_matrix <- mt_tally$nmf_matrix
}else if (sig.type=='CN_W'){
	load( file = paste(run.dir, "/result/sigminer/CN_W_tally.Rdata", sep=""))
	curr.nmf_matrix <- cn_tally$nmf_matrix
}else if (sig.type=='CN_S_48'){
	load( file = paste(run.dir, "/result/sigminer/CN_S_tally.Rdata", sep=""))
	curr.nmf_matrix <- cn_tally$all_matrices$CN_48
}else if (sig.type=='CN_S_40'){
	load( file = paste(run.dir, "/result/sigminer/CN_S_tally.Rdata", sep=""))
	curr.nmf_matrix <- cn_tally$all_matrices$CN_40
}

load(file = file.path( run.dir, "result", "sigminer", "sample.sex.Rdata" ))
options(sigminer.sex = sex.dt)

pdf( file=file.path(fig.dir, paste(sig.type, ".estimate.survey.pdf", sep="")), onefile=TRUE )

print(paste('Running ', sig.type, ' signature extraction', sep=''))
		
sig_est <- sig_unify_extract(
	curr.nmf_matrix,
	range = 1:10,
	nrun = 50L,
	approach = "bootstrap_nmf",
	n_bootstrap = 20L,
	cores = 24)

p1 <- bp_show_survey(sig_est, add_score=TRUE)
print(p1)
dev.off()

obj_suggested <- bp_get_sig_obj(sig_est, sig_est$suggested)
print(obj_suggested)

save(list=c("sig_est", "obj_suggested"),
	file = paste(run.dir, "/result/sigminer/extracted_signatures_", sig.type, ".Rdata", sep=""))

print(paste('Finished ', sig.type, ' signature extraction', sep=''))