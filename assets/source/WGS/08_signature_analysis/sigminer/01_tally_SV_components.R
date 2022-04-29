#!/usr/bin/env Rscript

## Not yet implemented
args <- commandArgs(TRUE)

run.dir <- args[1]

working_dir <- getwd()
gpl_prefix <- basename(working_dir)

working_dir.split <- strsplit(working_dir, split="/", fixed=TRUE)[[1]]
run.dir.idx <- which(working_dir.split=="sh") - 1
study.dir <- working_dir.split[[run.dir.idx]]
study.name <- strsplit(study.dir, split="_", fixed=TRUE)[[1]][[1]]

run.dir <- file.path( "~/workspace/runs", study.dir )

source( file.path( run.dir, "config_files/common_config.R" ) )

library(sigminer)
library(NMF)

germline.dir <- paste(run.dir, "/result/", gpl_prefix, sep="")
purple.dir <- paste(run.dir, "/result/", gpl_prefix, sep="")
ref.dir <- "/home/tjohnson/reference/HMF_38/dbs/ensembl_data_cache"

load( file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))

tumor.normal.pair.info <- tumor_normal_pair_info.GPL
tumor.sample.ids.ls <- tumor.normal.pair.info$sample.id.T
names(tumor.sample.ids.ls) <- tumor.sample.ids.ls
sample.ct <- length(tumor.sample.ids.ls)

tumor.normal.pair.info[,somatic.vcf:=file.path(purple.dir, subject.id, "purple", sample.id.T, paste(sample.id.T, ".purple.somatic.vcf.gz", sep="") )]
tumor.normal.pair.info[,cnv.file:=file.path(purple.dir, subject.id, "purple", sample.id.T, paste(sample.id.T, ".purple.segment.tsv", sep="") )]

load(file = file.path( run.dir, "result_summaries", "sigminer", "merged_CN.Rdata" ))
load(file = file.path( run.dir, "result_summaries", "sigminer", "sample.sex.Rdata" ))

options(sigminer.sex = sex.dt)

cn_tally_S <- sig_tally(merged.cn, method = "S", add_loh=TRUE, ref_genome="BSgenome.Hsapiens.UCSC.hg38", genome_build="hg38")
cn_tally_W <- sig_tally(merged.cn, method = "W", add_loh=TRUE, ref_genome="BSgenome.Hsapiens.UCSC.hg38", genome_build="hg38")

save(list=c("cn_tally_S", "cn_tally_W"),
	file = paste(run.dir, "/result/sigminer/CN_tally_W.Rdata", sep=""))
