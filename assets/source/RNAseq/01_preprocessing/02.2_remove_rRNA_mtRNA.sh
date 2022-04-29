#!/bin/bash

RUN_DIR=${RUN_DIR}
java_vmem=${java_vmem}

. ${RUN_DIR}/config_files/common_config.sh

## samtools and bcftools 1.12, cutadapt are in /home/tjohnson/.local
export PATH="/home/tjohnson/.local/bin:${PATH}"

module use /usr/local/package/modulefiles/
module load bwa/0.7.17

FASTQ_DIR="${RUN_DIR}/fastq"
READPAIR_INFO_FILE="${RUN_DIR}/config_files/fastq_readpair_by_lane_info.csv"

abundant_refs_dir="/home/tjohnson/reference/iGenomes/hg38/Homo_sapiens/UCSC/hg38/Sequence/AbundantSequences"

ABUNDANT_REF_LIST="${abundant_refs_dir}/chrM.fa,${abundant_refs_dir}/hum5SrDNA.fa,${abundant_refs_dir}/humRibosomal.fa"

bbmap_dir="/home/tjohnson/tools/BBMap/bbmap"

## First line is header
let "LINENUM = ${SGE_TASK_ID} + 1"

READPAIR_INFO=$(sed -n "${LINENUM}p" ${READPAIR_INFO_FILE})
SUBJECT_ID=$(echo "${READPAIR_INFO}" | awk -F ',' '{print $1}')
SAMPLE_ID=$(echo ${READPAIR_INFO} | awk -F "," '{print $2}')
SAMPLE_NAME=$(echo ${READPAIR_INFO} | awk -F "," '{print $3}')
FLOWCELL_ID=$(echo ${READPAIR_INFO} | awk -F "," '{print $4}')
SAMPLE_TYPE=$(echo ${READPAIR_INFO} | awk -F "," '{print $5}')
COMPRESSION_TYPE=$(echo ${READPAIR_INFO} | awk -F "," '{print $6}')
PLATFORM=$(echo ${READPAIR_INFO} | awk -F "," '{print $7}')
FLOWCELL_LANE=$(echo ${READPAIR_INFO} | awk -F "," '{print $8}')

echo "Running trim and map for ${SUBJECT_ID} and sample ${SAMPLE_ID}"

READGROUP_ID="${FLOWCELL_ID}.${SAMPLE_NAME}.${FLOWCELL_LANE}"
RG_PU=${READGROUP_ID}
READGROUP_STRING="@RG\\tID:${READGROUP_ID}\\tSM:${SAMPLE_ID}\\tLB:${SAMPLE_NAME}\\tPL:ILLUMINA\\tPU:${RG_PU}"

SUBJECT_FASTQ_DIR=${FASTQ_DIR}/${SUBJECT_ID}

input_r1_trimmed_fq="${SUBJECT_FASTQ_DIR}/${SAMPLE_ID}_${FLOWCELL_ID}_L00${FLOWCELL_LANE}_R1.trimmed.fastq"
input_r2_trimmed_fq="${SUBJECT_FASTQ_DIR}/${SAMPLE_ID}_${FLOWCELL_ID}_L00${FLOWCELL_LANE}_R2.trimmed.fastq"

output_r1_trimmed_m_fq="${SUBJECT_FASTQ_DIR}/${SAMPLE_ID}_${FLOWCELL_ID}_L00${FLOWCELL_LANE}_R1.trimmed.rRNA_mtRNA_mapped.fastq"
output_r2_trimmed_m_fq="${SUBJECT_FASTQ_DIR}/${SAMPLE_ID}_${FLOWCELL_ID}_L00${FLOWCELL_LANE}_R2.trimmed.rRNA_mtRNA_mapped.fastq"

output_r1_trimmed_um_fq="${SUBJECT_FASTQ_DIR}/${SAMPLE_ID}_${FLOWCELL_ID}_L00${FLOWCELL_LANE}_R1.trimmed.cleaned.fastq"
output_r2_trimmed_um_fq="${SUBJECT_FASTQ_DIR}/${SAMPLE_ID}_${FLOWCELL_ID}_L00${FLOWCELL_LANE}_R2.trimmed.cleaned.fastq"

statistics_file="${SUBJECT_FASTQ_DIR}/${SAMPLE_ID}_${FLOWCELL_ID}_L00${FLOWCELL_LANE}.cleaning_stats.txt"

echo "Running seal on ${input_r1_trimmed_fq} and ${input_r2_trimmed_fq}"
bash ${bbmap_dir}/seal.sh -Xmx${java_vmem} in=${input_r1_trimmed_fq} in2=${input_r2_trimmed_fq} ref=${ABUNDANT_REF_LIST} refstats=${statistics_file} out=${output_r1_trimmed_m_fq} out2=${output_r2_trimmed_m_fq} outu=${output_r1_trimmed_um_fq} outu2=${output_r2_trimmed_um_fq}


if [[ -s ${output_r1_trimmed_um_fq} ]] && [[ -s ${output_r2_trimmed_um_fq} ]]; then
	echo "Finished removing rRNA and mtRNA reads."
else
	echo "Read removal not performed.  Alignment aborted."
fi
