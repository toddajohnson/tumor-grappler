#!/bin/bash

RUN_DIR=${RUN_DIR}

if [[ -z ${load_biobambam2_from_tools} ]] ; then
	load_biobambam2_from_tools="false"
else
	load_biobambam2_from_tools=${load_biobambam2_from_tools}
fi

. ${RUN_DIR}/config_files/common_config.sh

blockmb=${BAMSORT_BLOCKMB}
bamsort_threads=${BAMSORT_THREADS}
bamsort_output_threads=${BAMSORT_OUTPUT_THREADS}

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


STAR_FASTQ_INFO_FILE="${RUN_DIR}/config_files/STAR_fastq_file_info.txt"

## First line is header
let "LINENUM = ${SGE_TASK_ID} + 1"

STAR_FASTQ_INFO=$(sed -n "${LINENUM}p" ${STAR_FASTQ_INFO_FILE})
SUBJECT_ID=$(echo "${STAR_FASTQ_INFO}" | awk -F '|' '{print $1}')
SAMPLE_ID=$(echo ${STAR_FASTQ_INFO} | awk -F "|" '{print $2}')

INPUT_PATH="${RUN_DIR}/result/STAR-${STAR_VERSION}/${SUBJECT_ID}"
OUTPUT_PATH="${RUN_DIR}/result/STAR-${STAR_VERSION}/${SUBJECT_ID}"

mkdir -p ${OUTPUT_PATH}/tmp

INPUT_FILE_PREFIX="${OUTPUT_PATH}/${SAMPLE_ID}"
TMP_FILE_PREFIX="${OUTPUT_PATH}/tmp/${SAMPLE_ID}"
OUTPUT_FILE_PREFIX="${OUTPUT_PATH}/${SAMPLE_ID}"

INPUT_BAM="${INPUT_FILE_PREFIX}Aligned.out.bam"

TMP_BAM_PATH="${TMP_FILE_PREFIX}.sorted.bam.tmp"
OUTPUT_BAM_PATH="${OUTPUT_FILE_PREFIX}Aligned.out.sorted.bam"
OUTPUT_BAM_INDEX_PATH="${OUTPUT_FILE_PREFIX}Aligned.out.sorted.bam.bai"

echo "Sorting bam files for ${SAMPLE_ID} RNA-seq data"

bamsort index=1 SO=coordinate level=1 \
	blockmb=${blockmb} \
	inputthreads=${bamsort_threads} outputthreads=${bamsort_output_threads} sortthreads=${bamsort_threads} \
	calmdnm=1 calmdnmrecompindentonly=1 calmdnmreference=${REFERENCE_GENOME_FA} \
	tmpfile=${TMP_BAM_PATH} \
	inputformat=sam \
	indexfilename=${OUTPUT_BAM_INDEX_PATH} \
	I=${INPUT_BAM} \
	O=${OUTPUT_BAM_PATH}

echo "Finished sorting ${SAMPLE_ID} RNA-seq data"
