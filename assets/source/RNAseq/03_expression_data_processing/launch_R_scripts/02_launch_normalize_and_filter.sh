#!/bin/bash

export LANG=en_US.UTF-8

RUN_DIR=${RUN_DIR}

. ${RUN_DIR}/config_files/common_config.sh
. ${HOME}/.R-4.1.0_setup.sh

echo "Launching normalize_and_filter R script for ${study_name}"
Rscript ${RNAseq_cg_pipeline_dir}/03_expression_data_processing/02_normalize_and_filter.R ${RUN_DIR}
