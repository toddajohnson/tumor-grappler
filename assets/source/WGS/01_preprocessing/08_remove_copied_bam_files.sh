#!/bin/bash

RUN_DIR=${RUN_DIR}

sample_index=${SGE_TASK_ID}
sample_type=${sample_type}

bam_dir="${RUN_DIR}/result/bam"
final_bam_dir="${RUN_DIR}/result/bam_final"
markdup_bam_dir="${RUN_DIR}/result/bam_markdup"
SAMPLE_INFO_FILE="${RUN_DIR}/config_files/tumor_normal_pair_info.csv"

mkdir -p ${final_bam_dir}

## First line is header
let "LINENUM = ${sample_index} + 1"
SAMPLE_INFO=$(sed -n "${LINENUM}p" ${SAMPLE_INFO_FILE})

if [[ ${sample_type} == "tumor" ]]; then
	sample_id=$(echo ${SAMPLE_INFO} | awk -F "," '{print $4}')
elif [[ ${sample_type} == "normal" ]]; then
	sample_id=$(echo ${SAMPLE_INFO} | awk -F "," '{print $6}')
fi

markdup_bam="${markdup_bam_dir}/${sample_id}.markdup.bam"
markdup_bai="${markdup_bam_dir}/${sample_id}.markdup.bam.bai"
markdup_md5="${markdup_bam_dir}/${sample_id}.markdup.bam.md5"
markdup_metrics="${markdup_bam_dir}/${sample_id}.markdup.metrics"

echo "Removing markdup bams to final bam directory for ${sample_id} at ${markdup_bam}"
rm ${markdup_bam}
rm ${markdup_bai}
rm ${markdup_md5} 
rm ${markdup_metrics}
