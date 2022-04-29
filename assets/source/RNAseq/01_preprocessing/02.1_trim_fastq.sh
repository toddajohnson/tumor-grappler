#!/bin/bash

RUN_DIR=${RUN_DIR}


if [[ -z ${use_temp_fastq} ]] ; then
	use_temp_fastq="false"
else
	use_temp_fastq=${use_temp_fastq}
fi

. ${RUN_DIR}/config_files/common_config.sh

bam_threads=${threads}

## samtools and bcftools 1.12, cutadapt are in /home/tjohnson/.local
export PATH="/home/tjohnson/.local/bin:${PATH}"

module use /usr/local/package/modulefiles/
module load bwa/0.7.17

FASTQ_DIR="${RUN_DIR}/fastq"
READPAIR_INFO_FILE="${RUN_DIR}/config_files/fastq_readpair_by_lane_info.csv"

base_ref_dir="/home/tjohnson/reference"

if [[ "${REFERENCE_GENOME}" == "38" ]] || [[ "${REFERENCE_GENOME}" == "HG38" ]] || [[ "${REFERENCE_GENOME}" == "GRCh38" ]]; then
	ref_dir="${base_ref_dir}/HMF/38"
	REFERENCE_GENOME_FA="${base_ref_dir}/HMF/38/refgenomes/Homo_sapiens.GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
elif [[ "${REFERENCE_GENOME}" == "37" ]] || [[ "${REFERENCE_GENOME}" == "HG37" ]] || [[ "${REFERENCE_GENOME}" == "GRCh38" ]] || [[ "${REFERENCE_GENOME}" == "hg19" ]]; then
	ref_dir="${base_ref_dir}/HMF/37"
	REFERENCE_GENOME_FA="${base_ref_dir}/HMF/37/refgenomes/Homo_sapiens.GRCh37/GRCh37.fa"
fi

# Cutadapt
CUTADAPT="cutadapt"
FWD_ADAPT=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
REV_ADAPT=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
CUTADAPT_OPTIONS="--minimum-length 1"


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
SEQ_RUN_DIR=$(echo ${READPAIR_INFO} | awk -F "," '{print $9}')
SAMPLE_FASTQ_DIR=$(echo ${READPAIR_INFO} | awk -F "," '{print $10}')
READ1_FASTQ_FILE=$(echo ${READPAIR_INFO} | awk -F "," '{print $11}')
READ2_FASTQ_FILE=$(echo ${READPAIR_INFO} | awk -F "," '{print $12}')


echo "Running trim and map for ${SUBJECT_ID} and sample ${SAMPLE_ID}"

READGROUP_ID="${FLOWCELL_ID}.${SAMPLE_NAME}.${FLOWCELL_LANE}"
RG_PU=${READGROUP_ID}
READGROUP_STRING="@RG\\tID:${READGROUP_ID}\\tSM:${SAMPLE_ID}\\tLB:${SAMPLE_NAME}\\tPL:ILLUMINA\\tPU:${RG_PU}"

SUBJECT_FASTQ_DIR=${FASTQ_DIR}/${SUBJECT_ID}
mkdir -p ${SUBJECT_FASTQ_DIR}

if [[ "${use_temp_fastq}" == "true" ]] ;  then
	TEMP_FASTQC_DIR="${RUN_DIR}/result/fastqc/tmp"
	
	TEMP_FASTQ_R1_PATH="${TEMP_FASTQC_DIR}/${SAMPLE_ID}_${FLOWCELL_ID}_L00${FLOWCELL_LANE}_R1.fastq"
	TEMP_FASTQ_R2_PATH="${TEMP_FASTQC_DIR}/${SAMPLE_ID}_${FLOWCELL_ID}_L00${FLOWCELL_LANE}_R2.fastq"
	
	FASTQ_R1_PATH="${TEMP_FASTQ_R1_PATH}"
	FASTQ_R2_PATH="${TEMP_FASTQ_R2_PATH}"
else
	SOURCE_FASTQ_DIR="${SEQ_RUN_DIR}/${SAMPLE_FASTQ_DIR}"
	
	FASTQ_R1_PATH="${SOURCE_FASTQ_DIR}/${READ1_FASTQ_FILE}"
	FASTQ_R2_PATH="${SOURCE_FASTQ_DIR}/${READ2_FASTQ_FILE}"
fi

output_r1_trimmed_fq="${SUBJECT_FASTQ_DIR}/${SAMPLE_ID}_${FLOWCELL_ID}_L00${FLOWCELL_LANE}_R1.trimmed.fastq"
output_r2_trimmed_fq="${SUBJECT_FASTQ_DIR}/${SAMPLE_ID}_${FLOWCELL_ID}_L00${FLOWCELL_LANE}_R2.trimmed.fastq"

if [[ ! -f ${output_r1_trimmed_fq} ]] || [[ ! -f ${output_r2_trimmed_fq} ]]; then
	if [[ -f ${FASTQ_R1_PATH} ]]; then
		if [[ "${PLATFORM}" == "Novaseq" ]] ; then
			echo "Running cutadapt with nextseq trimming on ${FASTQ_R1_PATH} and ${FASTQ_R2_PATH}"
			${CUTADAPT} -j ${bam_threads} -a ${FWD_ADAPT} -A ${REV_ADAPT} --nextseq-trim=20 ${CUTADAPT_OPTIONS} -o ${output_r1_trimmed_fq} -p ${output_r2_trimmed_fq} ${FASTQ_R1_PATH} ${FASTQ_R2_PATH}
		else
			echo "Running cutadapt on ${FASTQ_R1_PATH} and ${FASTQ_R2_PATH}"
			${CUTADAPT} -j ${bam_threads} -a ${FWD_ADAPT} -A ${REV_ADAPT} ${CUTADAPT_OPTIONS} -o ${output_r1_trimmed_fq} -p ${output_r2_trimmed_fq} ${FASTQ_R1_PATH} ${FASTQ_R2_PATH}
		fi
	fi
fi

if [[ -s ${output_r1_trimmed_fq} ]] && [[ -s ${output_r2_trimmed_fq} ]]; then
	echo "Finished adapter trimming."
else
	echo "Trimmed fastq files were not created.  Alignment aborted."
fi
