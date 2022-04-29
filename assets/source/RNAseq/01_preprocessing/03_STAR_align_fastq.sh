#!/bin/bash

RUN_DIR=${RUN_DIR}
threads=${threads}

. ${RUN_DIR}/config_files/common_config.sh

## STAR 2.7.9a, samtools and bcftools 1.12, cutadapt are in /home/tjohnson/.local
export PATH="/home/tjohnson/.local/bin:${PATH}"

#STAR_FASTQ_INFO_FILE="${RUN_DIR}/config_files/STAR_fastq_file_info.csv"
STAR_FASTQ_INFO_FILE="${RUN_DIR}/config_files/STAR_fastq_file_info.txt"

## First line is header
let "LINENUM = ${SGE_TASK_ID} + 1"

STAR_FASTQ_INFO=$(sed -n "${LINENUM}p" ${STAR_FASTQ_INFO_FILE})
SUBJECT_ID=$(echo "${STAR_FASTQ_INFO}" | awk -F '|' '{print $1}')
SAMPLE_ID=$(echo ${STAR_FASTQ_INFO} | awk -F "|" '{print $2}')
readFilesIn=$(echo ${STAR_FASTQ_INFO} | awk -F "|" '{print $3}')
outSAMattrRGline=$(echo ${STAR_FASTQ_INFO} | awk -F "|" '{print $4}')

#SUBJECT_ID=$(echo "${STAR_FASTQ_INFO}" | awk -F ',' '{print $1}')
#SAMPLE_ID=$(echo ${STAR_FASTQ_INFO} | awk -F "," '{print $2}')
#readFilesIn=$(echo ${STAR_FASTQ_INFO} | awk -F "," '{print $3}')
#outSAMattrRGline=$(echo ${STAR_FASTQ_INFO} | awk -F "," '{print $4}')

OUTPUT_PATH="${RUN_DIR}/result/STAR-${STAR_VERSION}/${SUBJECT_ID}"
mkdir -p ${OUTPUT_PATH}

OUTPUT_FILE_PREFIX="${OUTPUT_PATH}/${SAMPLE_ID}"

readFastqCmd="${readFilesIn} ${outSAMattrRGline}"

echo "Running genome mapping of ${SAMPLE_ID} RNA-seq data using STAR-${STAR_VERSION}"

STAR --runThreadN ${threads} \
	--genomeDir ${STAR_GENOME_INDICES_DIR} \
	${readFastqCmd} \
	--outFileNamePrefix ${OUTPUT_FILE_PREFIX} \
	--outSAMtype BAM Unsorted \
	--outSAMunmapped Within \
	--outBAMcompression 0 \
	--outSAMattributes All \
	--outFilterMultimapNmax 10 \
	--outFilterMismatchNmax 3 limitOutSJcollapsed 3000000 \
	--chimSegmentMin 10 \
	--chimOutType WithinBAM SoftClip  \
	--chimJunctionOverhangMin 10 \
	--chimSegmentReadGapMax 3 \
	--chimScoreMin 1 \
	--chimScoreDropMax 30 \
	--chimScoreJunctionNonGTAG 0 \
	--chimScoreSeparation 1 \
	--outFilterScoreMinOverLread 0.33 \
	--outFilterMatchNminOverLread 0.33 \
	--outFilterMatchNmin 35 \
	--alignSplicedMateMapLminOverLmate 0.33 \
	--alignSplicedMateMapLmin 35 \
	--alignSJstitchMismatchNmax 5 -1 5 5

echo "Finished genome mapping of ${SAMPLE_ID} RNA-seq data"
