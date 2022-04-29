#!/bin/bash

export LANG=en_US.UTF-8

RUN_DIR="/home/tjohnson/workspace/runs/IWK_WGS_HMF_20210726"
RNAseq_RUN_DIR="/home/tjohnson/workspace/runs/IWK_RNAseq_20210805"
cg_pipeline_dir="/home/tjohnson/workspace/scripts/cancer_genomics_pipeline/v1/WGS"
RNAseq_cg_pipeline_dir="/home/tjohnson/workspace/scripts/cancer_genomics_pipeline/v1/RNAseq"
pvac_pipeline_dir="/home/tjohnson/workspace/scripts/cancer_genomics_pipeline/v1/pVAC"

study_name="IWK"

adapters_file="none"

REFERENCE_GENOME="HG38"

STAR_VERSION="2.7.9a"

SAGE_VERSION=2.8.mod

AMBER_VERSION=3.5
COBALT_VERSION=1.11
GRIDSS_VERSION=2.12.0
#GRIPSS_VERSION=kt_v1.12
GRIPSS_VERSION=v2.0_beta
PURPLE_VERSION=3.1.mod
LINX_VERSION=1.17

install_dir="/home/tjohnson/tools/gridss-${GRIDSS_VERSION}-purple-linx"

export SAGE_JAR=${install_dir}/hmftools/sage-${SAGE_VERSION}.jar
export AMBER_JAR=$install_dir/hmftools/amber-${AMBER_VERSION}.jar
export COBALT_JAR=$install_dir/hmftools/cobalt-${COBALT_VERSION}.jar
export GRIDSS_JAR=$install_dir/gridss/gridss-${GRIDSS_VERSION}-gridss-jar-with-dependencies.jar
#export GRIPSS_JAR=$install_dir/hmftools/gripss-${GRIPSS_VERSION}.jar
export GRIPSS_JAR=$install_dir/hmftools/gripss_${GRIPSS_VERSION}.jar
export PURPLE_JAR=$install_dir/hmftools/purple_v${PURPLE_VERSION}.jar
export LINX_JAR=$install_dir/hmftools/linx_v${LINX_VERSION}.jar

## used for labelling directories/files
gridss_version="2.12.0"
sage_version="2.8"
SAGE_VERSION=2.8

gpl_prefix="GRIDSS-2.12.0"

snpeff_anno_db="GRCh38.manual.104"

#min_tumor_qual=100
#min_tumor_vaf=0.025
#max_tumor_only_vaf=0.975
#min_tumor_only_depth=10
#min_germline_depth=10
#max_germline_vaf=0.04
## germline_rel_AF is not an HMF setting
## set as AF.N/AF.T in R or FORMAT/AF[0]/FORMAT/AF[1] in bash
#max_germline_rel_AF=0.10
#max_germline_rel_raw_base_qual=0.04

GISTIC2_install_dir=~/tools/GISTIC2