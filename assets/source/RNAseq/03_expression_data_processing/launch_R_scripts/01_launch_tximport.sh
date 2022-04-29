#!/bin/bash

RUN_DIR=${RUN_DIR}

. ${RUN_DIR}/config_files/common_config.sh
. ${HOME}/.R-4.1.0_setup.sh

echo "Launching tximport R script for ${study_name}"
Rscript ${RNAseq_cg_pipeline_dir}/03_expression_data_processing/01_tximport.R ${RUN_DIR}