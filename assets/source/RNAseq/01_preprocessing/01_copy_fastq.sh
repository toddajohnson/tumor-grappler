#!/bin/bash

RUN_DIR=${RUN_DIR}

. ${RUN_DIR}/config_files/common_config.sh

read_num=${READ_NUM}

FASTQ_DIR="${RUN_DIR}/fastq"
READPAIR_INFO_FILE="${RUN_DIR}/config_files/fastq_readpair_by_lane_info.csv"

## First line is header
let "LINENUM = ${SGE_TASK_ID} + 1"

READPAIR_INFO=$(sed -n "${LINENUM}p" ${READPAIR_INFO_FILE})
SUBJECT_ID=$(echo "${READPAIR_INFO}" | awk -F ',' '{print $1}')
SAMPLE_ID=$(echo ${READPAIR_INFO} | awk -F "," '{print $2}')
SAMPLE_NAME=$(echo ${READPAIR_INFO} | awk -F "," '{print $3}')
FLOWCELL_ID=$(echo ${READPAIR_INFO} | awk -F "," '{print $4}')
SAMPLE_TYPE=$(echo ${READPAIR_INFO} | awk -F "," '{print $5}')
COMPRESSION_TYPE=$(echo ${READPAIR_INFO} | awk -F "," '{print $6}')
FLOWCELL_LANE=$(echo ${READPAIR_INFO} | awk -F "," '{print $8}')
SEQ_RUN_DIR=$(echo ${READPAIR_INFO} | awk -F "," '{print $9}')
SAMPLE_FASTQ_DIR=$(echo ${READPAIR_INFO} | awk -F "," '{print $10}')
READ1_FASTQ_FILE=$(echo ${READPAIR_INFO} | awk -F "," '{print $11}')
READ2_FASTQ_FILE=$(echo ${READPAIR_INFO} | awk -F "," '{print $12}')

if [[ ${read_num} == 1 ]]; then
	FASTQ_FILE=${READ1_FASTQ_FILE}
elif [[ ${read_num} == 2 ]]; then
	FASTQ_FILE=${READ2_FASTQ_FILE}
fi

SOURCE_FASTQ_DIR="${SEQ_RUN_DIR}/${SAMPLE_FASTQ_DIR}"
FASTQ_PATH="${SOURCE_FASTQ_DIR}/${FASTQ_FILE}"

OUTPUT_PATH="${RUN_DIR}/fastq"

mkdir -p ${OUTPUT_PATH}

if [[ -f ${FASTQ_PATH} ]]; then
	echo "Copying fastq file for ${SAMPLE_ID} on ${FASTQ_PATH}"
	cp ${FASTQ_PATH} ${OUTPUT_PATH}
else
	echo "The fastq file for ${SAMPLE_ID} was not found at ${FASTQ_PATH}."
fi
