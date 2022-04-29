#!/bin/bash

RUN_DIR=${RUN_DIR}


if [[ -z ${use_temp_fastq} ]] ; then
	use_temp_fastq="false"
else
	use_temp_fastq=${use_temp_fastq}
fi

if [[ -z ${compressed_file} ]] ; then
	compressed_file="false"
else
	compressed_file=${compressed_file}
fi

if [[ -z ${load_biobambam2_from_tools} ]] ; then
	load_biobambam2_from_tools="false"
else
	load_biobambam2_from_tools=${load_biobambam2_from_tools}
fi

. ${RUN_DIR}/config_files/common_config.sh

blockmb=${BAMSORT_BLOCKMB}
bam_threads=${BAM_THREADS}
bamsort_threads=${BAMSORT_THREADS}
bamsort_output_threads=${BAMSORT_OUTPUT_THREADS}

## samtools and bcftools 1.12, cutadapt are in /home/tjohnson/.local
export PATH="/home/tjohnson/.local/bin:${PATH}"

module use /usr/local/package/modulefiles/
module load bwa/0.7.17

if [[ "${load_biobambam2_from_tools}" == "false" ]] ; then
	echo "Loading biobambam2-2.0.146 from module"
	module load biobambam/2.0.146
else
	echo "Adding biobambam2-2.0.182 to PATH"
	module load gcc/10.2.0
	export LD_LIBRARY_PATH="${HOME}/tools/libdeflate-1.7:${HOME}/tools/snappy-1.1.8/lib64:${HOME}/tools/io_lib-1.14.14/lib:${LD_LIBRARY_PATH}"
	export PATH="${HOME}/tools/biobambam2-2.0.182/bin:${HOME}/.local/bin:${PATH}"
fi

READPAIR_INFO_FILE="${RUN_DIR}/config_files/fastq_readpair_by_lane_info.csv"

base_ref_dir="/home/tjohnson/reference"

if [[ "${REFERENCE_GENOME}" == "38" ]] || [[ "${REFERENCE_GENOME}" == "HG38" ]] || [[ "${REFERENCE_GENOME}" == "GRCh38" ]]; then
	ref_dir="${base_ref_dir}/HMF/38"
	REFERENCE_GENOME_FA="${base_ref_dir}/HMF/38/refgenomes/Homo_sapiens.GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
elif [[ "${REFERENCE_GENOME}" == "37" ]] || [[ "${REFERENCE_GENOME}" == "HG37" ]] || [[ "${REFERENCE_GENOME}" == "hg19" ]]; then
	ref_dir="${base_ref_dir}/HMF/37"
	REFERENCE_GENOME_FA="${base_ref_dir}/HMF/37/refgenomes/Homo_sapiens.GRCh37/GRCh37.fa"
fi


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
SAMPLE_FCID_LANE_INDEX=$(echo ${READPAIR_INFO} | awk -F "," '{print $13}')

echo "Running trim and map for ${SUBJECT_ID} and sample ${SAMPLE_ID}"

READGROUP_ID="${FLOWCELL_ID}.${SAMPLE_NAME}.${FLOWCELL_LANE}"
RG_PU=${READGROUP_ID}
READGROUP_STRING="@RG\\tID:${READGROUP_ID}\\tSM:${SAMPLE_ID}\\tLB:${SAMPLE_NAME}\\tPL:ILLUMINA\\tPU:${RG_PU}"

if [[ "${use_temp_fastq}" == "true" ]] ;  then
	TEMP_FASTQC_DIR="${RUN_DIR}/result/fastqc/tmp"
	
	echo "Running analysis using ${TEMP_FASTQC_DIR}"
	
	if [[ ${SAMPLE_FCID_LANE_INDEX} == 1 ]] ;  then
		TEMP_FASTQ_R1_PATH="${TEMP_FASTQC_DIR}/${SAMPLE_ID}_${FLOWCELL_ID}_L00${FLOWCELL_LANE}_R1.fastq"
		TEMP_FASTQ_R2_PATH="${TEMP_FASTQC_DIR}/${SAMPLE_ID}_${FLOWCELL_ID}_L00${FLOWCELL_LANE}_R2.fastq"
	elif [[ ${SAMPLE_FCID_LANE_INDEX} == 2 ]] ;  then
		TEMP_FASTQ_R1_PATH="${TEMP_FASTQC_DIR}/${SAMPLE_ID}_${FLOWCELL_ID}_L00${FLOWCELL_LANE}_${SAMPLE_FCID_LANE_INDEX}_R1.fastq"
		TEMP_FASTQ_R2_PATH="${TEMP_FASTQC_DIR}/${SAMPLE_ID}_${FLOWCELL_ID}_L00${FLOWCELL_LANE}_${SAMPLE_FCID_LANE_INDEX}_R2.fastq"
	else
		echo "The sample FCID lane index that was passed did not make sense. Exiting."
		exit 1
	fi

	if [[ "${compressed_file}" == "true" ]] ;  then
		echo "Adding compressed file suffix"
		TEMP_FASTQ_R1_PATH="${TEMP_FASTQ_R1_PATH}.${COMPRESSION_TYPE}"
		TEMP_FASTQ_R2_PATH="${TEMP_FASTQ_R2_PATH}.${COMPRESSION_TYPE}"
	fi
	
	FASTQ_R1_PATH="${TEMP_FASTQ_R1_PATH}"
	FASTQ_R2_PATH="${TEMP_FASTQ_R2_PATH}"
else
	echo "Using original source fastq files"
	SOURCE_FASTQ_DIR="${SEQ_RUN_DIR}/${SAMPLE_FASTQ_DIR}"
	
	FASTQ_R1_PATH="${SOURCE_FASTQ_DIR}/${READ1_FASTQ_FILE}"
	FASTQ_R2_PATH="${SOURCE_FASTQ_DIR}/${READ2_FASTQ_FILE}"
fi

mkdir -p ${RUN_DIR}/result/bam/tmp

if [[ ${SAMPLE_FCID_LANE_INDEX} == 1 ]] ;  then
	TMP_BAM_PATH="${RUN_DIR}/result/bam/tmp/${SAMPLE_ID}_${FLOWCELL_ID}_L00${FLOWCELL_LANE}.sorted.bam.tmp"
	OUTPUT_BAM_PATH="${RUN_DIR}/result/bam/${SAMPLE_ID}_${FLOWCELL_ID}_L00${FLOWCELL_LANE}.sorted.bam"
	OUTPUT_BAM_INDEX_PATH="${RUN_DIR}/result/bam/${SAMPLE_ID}_${FLOWCELL_ID}_L00${FLOWCELL_LANE}.sorted.bam.bai"
elif [[ ${SAMPLE_FCID_LANE_INDEX} == 2 ]] ;  then
	TMP_BAM_PATH="${RUN_DIR}/result/bam/tmp/${SAMPLE_ID}_${FLOWCELL_ID}_L00${FLOWCELL_LANE}_${SAMPLE_FCID_LANE_INDEX}.sorted.bam.tmp"
	OUTPUT_BAM_PATH="${RUN_DIR}/result/bam/${SAMPLE_ID}_${FLOWCELL_ID}_L00${FLOWCELL_LANE}_${SAMPLE_FCID_LANE_INDEX}.sorted.bam"
	OUTPUT_BAM_INDEX_PATH="${RUN_DIR}/result/bam/${SAMPLE_ID}_${FLOWCELL_ID}_L00${FLOWCELL_LANE}_${SAMPLE_FCID_LANE_INDEX}.sorted.bam.bai"
else
	echo "The sample FCID lane index that was passed did not make sense. Exiting."
	exit 1
fi

if [[ -s ${FASTQ_R1_PATH} ]] && [[ -s ${FASTQ_R2_PATH} ]]; then
	echo "Running bwa mem and bamsort"
	echo "FASTQ1: ${FASTQ_R1_PATH}"
	echo "FASTQ2: ${FASTQ_R2_PATH}"

	if [[ "${compressed_file}" == "true" ]] && [[ ${COMPRESSION_TYPE} == "bz2" ]] ; then
		echo "Running bwa-mem with redirection of bzip2 processed fastq files"
		bwa mem -R ${READGROUP_STRING} -Y -t ${bam_threads} ${REFERENCE_GENOME_FA} <(bzip2 -dc ${FASTQ_R1_PATH}) <(bzip2 -dc ${FASTQ_R2_PATH}) | \
		bamsort index=1 SO=coordinate level=1 \
			blockmb=${blockmb} \
			inputthreads=${bamsort_threads} outputthreads=${bamsort_output_threads} sortthreads=${bamsort_threads} \
			calmdnm=1 calmdnmrecompindentonly=1 calmdnmreference=${REFERENCE_GENOME_FA} \
			tmpfile=${TMP_BAM_PATH} \
			inputformat=sam \
			indexfilename=${OUTPUT_BAM_INDEX_PATH} \
			O=${OUTPUT_BAM_PATH}
	else
		echo "Running bwa-mem with standard arguments"
		bwa mem -R ${READGROUP_STRING} -Y -t ${bam_threads} ${REFERENCE_GENOME_FA} ${FASTQ_R1_PATH} ${FASTQ_R2_PATH} | \
		bamsort index=1 SO=coordinate level=1 \
			blockmb=${blockmb} \
			inputthreads=${bamsort_threads} outputthreads=${bamsort_output_threads} sortthreads=${bamsort_threads} \
			calmdnm=1 calmdnmrecompindentonly=1 calmdnmreference=${REFERENCE_GENOME_FA} \
			tmpfile=${TMP_BAM_PATH} \
			inputformat=sam \
			indexfilename=${OUTPUT_BAM_INDEX_PATH} \
			O=${OUTPUT_BAM_PATH}
	fi
			
	if [[ -s ${OUTPUT_BAM_PATH} ]]; then
		echo "BAM created.  Removing intermediate fastq files."
		if [[ "${use_temp_fastq}" == "true" ]] ;  then
			echo "Removing temp fastq files from ${TEMP_FASTQ_R1_PATH} and ${TEMP_FASTQ_R2_PATH}"
			rm ${TEMP_FASTQ_R1_PATH}
			rm ${TEMP_FASTQ_R2_PATH}
		fi
	else
		echo "BAM file was not created. Keeping intermediate files"
	fi
fi

echo "Finished mapping.  Output to ${OUTPUT_BAM_PATH}"
