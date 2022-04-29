Sys.setenv(LANG = "en_US.UTF-8")
Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

run.dir <- '/home/tjohnson/workspace/runs/IWK_RNAseq_20210805'
wgs.run.dir <- '/home/tjohnson/workspace/runs/IWK_WGS_HMF_20210726'
base.dir <- run.dir
cg_pipeline_dir <- "/home/tjohnson/workspace/scripts/cancer_genomics_pipeline/v1/WGS"
RNAseq_cg_pipeline_dir <- "/home/tjohnson/workspace/scripts/cancer_genomics_pipeline/v1/RNAseq"

snpEff.db <- "GRCh38.manual.104"

run.dirs.ls <- c(run.dir)

study.name <- "IWK"