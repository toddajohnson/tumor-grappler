#!/bin/bash

export LANG=en_US.UTF-8

RUN_DIR="/home/tjohnson/workspace/runs/IWK_RNAseq_20210805"
cg_pipeline_dir="/home/tjohnson/workspace/scripts/cancer_genomics_pipeline/v1/WGS"
RNAseq_cg_pipeline_dir="/home/tjohnson/workspace/scripts/cancer_genomics_pipeline/v1/RNAseq"

adapters_file="${cg_pipeline_dir}/common_config_files/adapters.tsv"

REFERENCE_GENOME="HG38"
#REFERENCE_GENOME="GRCh38"
REFERENCE_GENOME_FA="/home/tjohnson/reference/HMF_GRCh38/refgenomes/Homo_sapiens.GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

STAR_VERSION="2.7.9a"
OVERHANG_BP=125
STAR_GENOME_INDICES_DIR="/home/tjohnson/reference/GRCh38/STAR-2.7.9a_genome_indices_${OVERHANG_BP}bp/GCA_000001405.15_GRCh38_no_alt_analysis_set_gencode_v38"

install_dir="/home/tjohnson/tools"
#ISOFOX_VERSION=1.2
#ISOFOX_VERSION="1.3_beta"

# isofox now uses ensembl cache, so should not need to modify jar
#ISOFOX_VERSION="1.3"
ISOFOX_VERSION="1.4"
export ISOFOX_JAR=$install_dir/hmftools/isofox_v${ISOFOX_VERSION}.jar
#export ISOFOX_JAR=$install_dir/hmftools/isofox_v${ISOFOX_VERSION}.mod.jar

ref_dir="/home/tjohnson/reference/HMF/38"
ensembl_cache_dir="${ref_dir}/dbs/linx"

## excluded rRNA genes within coding genes may cause those to be removed too
## specifically, saw e-cadherin CDH1 was not in latest output
## probably, pre-proessing to remove mtDNA and rRNA genes befor STAR will remove that problem
excluded_genes_list="${RNAseq_cg_pipeline_dir}/common_config_files/excluded_genes_20220402.csv"

#isofox_dir="${ref_dir}/dbs/isofox/without_excluded_genes"
isofox_dir="${ref_dir}/dbs/isofox"
mkdir -p ${isofox_dir}

read_length=126

study_name="IWK"