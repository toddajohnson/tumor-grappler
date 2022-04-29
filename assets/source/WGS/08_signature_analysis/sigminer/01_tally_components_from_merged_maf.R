#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

run.dir <- args[1]
sig.type <- args[2]

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

## load( file = paste(run.dir, "/config_files/fastq_file_info.Rdata", sep=""))
## 
## tumor.normal.pair.info <- tumor_normal_pair_info.GPL
## tumor.sample.ids.ls <- tumor.normal.pair.info$sample.id.T
## names(tumor.sample.ids.ls) <- tumor.sample.ids.ls
## sample.ct <- length(tumor.sample.ids.ls)
## 
## tumor.normal.pair.info[,somatic.vcf:=file.path(purple.dir, subject.id, "purple", sample.id.T, paste(sample.id.T, ".purple.somatic.vcf.gz", sep="") )]
## tumor.normal.pair.info[,cnv.file:=file.path(purple.dir, subject.id, "purple", sample.id.T, paste(sample.id.T, ".purple.segment.tsv", sep="") )]

tumor.mafs <- readRDS( file = file.path( run.dir, "result", "maftools", "merged_mafs_nonMT.RDS" ))

# when to use add_trans_bias?
mt_tally <- sig_tally(
	tumor.mafs,
	ref_genome = "BSgenome.Hsapiens.UCSC.hg38",
	use_syn = TRUE,
	mode = sig.type,
#	add_trans_bias = TRUE
	add_trans_bias = FALSE
)

save(list=c("mt_tally"),
	file = paste(run.dir, "/result/sigminer/mt_tally_", sig.type, ".Rdata", sep=""))
