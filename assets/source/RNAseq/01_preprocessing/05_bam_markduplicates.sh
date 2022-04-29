#!/bin/bash

RUN_DIR=${RUN_DIR}

if [[ -z ${load_biobambam2_from_tools} ]] ; then
	load_biobambam2_from_tools="false"
else
	load_biobambam2_from_tools=${load_biobambam2_from_tools}
fi

. ${RUN_DIR}/config_files/common_config.sh

sample_index=${SGE_TASK_ID}
BMD_THREAD_NUM=${THREADS}

## Corresponding file for WGS was much more complicated, but this first RNA-seq data comes from only one lane on one flowcell.
## May need to be modified if more complicated RNA-seq is used.

## STAR 2.7.9a, samtools and bcftools 1.12, cutadapt are in /home/tjohnson/.local
module use /usr/local/package/modulefiles/

if [[ "${load_biobambam2_from_tools}" == "false" ]] ; then
	export PATH="${HOME}/.local/bin:${PATH}"
	echo "Loading biobambam2-2.0.146 from module"
	module load biobambam/2.0.146
else
	echo "Adding biobambam2-2.0.182 to PATH"
	module load gcc/10.2.0
	export LD_LIBRARY_PATH="${HOME}/tools/libdeflate-1.7:${HOME}/tools/snappy-1.1.8/lib64:${HOME}/tools/io_lib-1.14.14/lib:${LD_LIBRARY_PATH}"
	export PATH="${HOME}/tools/biobambam2-2.0.182/bin:${HOME}/.local/bin:${PATH}"
fi

SAMPLE_INFO_FILE="${RUN_DIR}/config_files/sample_info.csv"

## First line is header
let "LINENUM = ${sample_index} + 1"
SAMPLE_INFO=$(sed -n "${LINENUM}p" ${SAMPLE_INFO_FILE})
subject_id=$(echo ${SAMPLE_INFO} | awk -F "," '{print $1}')
sample_id=$(echo ${SAMPLE_INFO} | awk -F "," '{print $2}')

bam_dir="${RUN_DIR}/result/STAR-${STAR_VERSION}/${subject_id}"
tmp_bam_dir="${RUN_DIR}/result/STAR-${STAR_VERSION}/${subject_id}/tmp"

mkdir -p ${bam_markdup} ${tmp_bam_dir}

input_file_prefix="${bam_dir}/${sample_id}"
markdup_input_paths="I=${input_file_prefix}Aligned.out.sorted.bam"

echo "Marking duplicates on ${sample_id}"

markdup_bam_path="${bam_dir}/${sample_id}.bam"
metrics_path="${bam_dir}/${sample_id}.metrics"
tmp_path="${tmp_bam_dir}/${sample_id}.markdup.tmp"

echo "Marking duplicated for ${markdup_input_paths}"

bammarkduplicates M=${metrics_path} tmpfile=${tmp_path} markthreads=${BMD_THREAD_NUM} rewritebam=1 rewritebamlevel=1 index=1 md5=1 ${markdup_input_paths} O=${markdup_bam_path}

echo "Finished marking duplicated into ${markdup_bam_path}"
