#!/bin/bash

threads=${threads}
RUN_DIR=${RUN_DIR}

## samtools and bcftools 1.12, cutadapt are in /home/tjohnson/.local
export PATH="/home/tjohnson/.local/bin:${PATH}"

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
SAMPLE_FCID_LANE_INDEX=$(echo ${READPAIR_INFO} | awk -F "," '{print $13}')

mkdir -p ${RUN_DIR}/result/bam/stats

if [[ ${SAMPLE_FCID_LANE_INDEX} == 1 ]] ;  then
	FLOWCELL_LANE=${FLOWCELL_LANE}
elif [[ ${SAMPLE_FCID_LANE_INDEX} == 2 ]] ;  then
	FLOWCELL_LANE=${FLOWCELL_LANE}_${SAMPLE_FCID_LANE_INDEX}
else
	echo "The sample FCID lane index that was passed did not make sense. Exiting."
	exit 1
fi

input_bam_path="${RUN_DIR}/result/bam/${SAMPLE_ID}_${FLOWCELL_ID}_L00${FLOWCELL_LANE}.sorted.bam"
output_bam_stats="${RUN_DIR}/result/bam/stats/${SAMPLE_ID}_${FLOWCELL_ID}_L00${FLOWCELL_LANE}.sorted.bamstats"
	
if [[ -f ${output_bam_stats} ]]; then
	echo "Removing previous stats file: ${output_bam_stats}"
	rm ${output_bam_stats}
fi

if [[ -f ${input_bam_path} ]] ; then
	echo "Running samtools stats on ${input_bam_path}"
	samtools stats --threads ${threads} ${input_bam_path} > ${output_bam_stats}
else
	echo "There was no file ${input_bam_path}"
fi