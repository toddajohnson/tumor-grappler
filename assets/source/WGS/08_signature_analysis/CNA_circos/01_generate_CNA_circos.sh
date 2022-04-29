#!/bin/bash

run_dir="${RUN_DIR}"

. ${run_dir}/config_files/common_config.sh

#module use /usr/local/package/modulefiles/
#module load java/11

export PATH="/home/tjohnson/.local/bin:/home/tjohnson/tools/circos-0.69-9/bin:${PATH}"

COHORT_INFO_FILE="${run_dir}/config_files/CNA_circos_sample_info.csv"

## First line is header
let "LINENUM = ${SGE_TASK_ID} + 1"

COHORT_INFO=$(sed -n "${LINENUM}p" ${COHORT_INFO_FILE})

cancerType=$(echo ${COHORT_INFO} | awk -F "," '{print $1}')

cd ${OUTPUT_PATH}

conf_file=${run_dir}/result/CNA_circos/copyNumberSummary/${cancerType}.conf
output_dir=${run_dir}/result/CNA_circos/plots
output_fig_file=${cancerType}.png

circos -nosvg -debug_group text -debug_group textplace -conf ${conf_file} -outputdir ${output_dir} -outputfile ${output_fig_file}
