Sys.setenv(LANG = "en_US.UTF-8")
Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

run.dir <- '/home/tjohnson/workspace/runs/IWK_WGS_HMF_20210726'
rnaseq.run.dir="/home/tjohnson/workspace/runs/IWK_RNAseq_20210805"
base.dir <- run.dir
cg_pipeline_dir <- "/home/tjohnson/workspace/scripts/cancer_genomics_pipeline/v1/WGS"

study.name <- "IWK"
gpl_prefix <- "GRIDSS-2.12.0"
gridss.version <- "2.12.0"
sage_version <- "2.8"
snpEff.db <- "GRCh38.manual.104"

run.dirs.ls <- c(run.dir)
names(run.dirs.ls) <- basename(run.dirs.ls)

import_panel_germline <- TRUE
import_genomewide_germline <- TRUE

maf.cutoff <- 0.01
af.cutoff <- 0.99

import.isofox.rna.fusion.match <- FALSE
tsv.sep <- ','

extra.samples.to.exclude <- c('IWK041_T')
exclude.from.CN.analysis <- c('IWK011_T', 'IWK041_T')

remove.copy.number.noise.samples <- FALSE
remove.deleted.genes.samples <- FALSE
remove.gender.mismatch.samples <- FALSE
tumor.depth.cutoff <- 5

version.date.for.driver.catalog.file <- 'ver.20220210'
version.date.for.CN.segments.file <- 'ver_20220131'