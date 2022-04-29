#!/bin/bash

RUN_DIR=${RUN_DIR}

if [[ -z ${use_mv} ]] ; then
	use_mv="false"
else
	use_mv=${use_mv}
fi

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

final_bam="${final_bam_dir}/${sample_id}.bam"
final_bai="${final_bam_dir}/${sample_id}.bam.bai"
final_md5="${final_bam_dir}/${sample_id}.bam.md5"
final_metrics="${final_bam_dir}/${sample_id}.metrics"
markdup_bam="${markdup_bam_dir}/${sample_id}.markdup.bam"
markdup_bai="${markdup_bam_dir}/${sample_id}.markdup.bam.bai"
markdup_md5="${markdup_bam_dir}/${sample_id}.markdup.bam.md5"
markdup_metrics="${markdup_bam_dir}/${sample_id}.markdup.metrics"

if [[ ! -f ${final_bam} ]]; then
	if [[ "${use_mv}" == "false" ]] ; then
		echo "Copying markdup bams to final bam directory for ${sample_id} at ${final_bam}"
		cp ${markdup_bam} ${final_bam}
		cp ${markdup_bai} ${final_bai}
		cp ${markdup_md5} ${final_md5}
		cp ${markdup_metrics} ${final_metrics}
	else
		echo "Moving markdup bams to final bam directory for ${sample_id} at ${final_bam}"
		mv ${markdup_bam} ${final_bam}
		mv ${markdup_bai} ${final_bai}
		mv ${markdup_md5} ${final_md5}
		mv ${markdup_metrics} ${final_metrics}
	fi
	echo "Removing intermediate sorted bam files at ${bam_dir}/${sample_id}_*.sorted.bam"
	rm ${bam_dir}/${sample_id}_*.sorted.bam
	rm ${bam_dir}/${sample_id}_*.sorted.bam.bai
else
	echo "Final bam file already exists"
fi
