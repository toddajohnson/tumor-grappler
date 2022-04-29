#!/bin/bash

export LANG=en_US.UTF-8

RUN_DIR=${RUN_DIR}
threads=${threads}
sample_index=${SGE_TASK_ID}

. ${RUN_DIR}/config_files/common_config.sh

module use /usr/local/package/modulefiles/
module load R/4.0.2

Rscript ${cg_pipeline_dir}/04_gridss_purple_linx/06B_extract_purple_germline_genomewide_vcf.R ${RUN_DIR} ${threads} ${sample_index}