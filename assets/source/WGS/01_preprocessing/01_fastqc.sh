#!/bin/bash

RUN_DIR=${RUN_DIR}
read_num=${READ_NUM}
threads=${threads}

jvmheap="25G"

if [[ -z ${use_temp} ]] ; then
	${use_temp}="true"
else
	use_temp=${use_temp}
fi

if [[ -z ${decompress_file} ]] ; then
	decompress_file="true"
else
	decompress_file=${decompress_file}
fi

. ${RUN_DIR}/config_files/common_config.sh

## samtools and bcftools 1.12, cutadapt are in /home/tjohnson/.local
export PATH="/home/tjohnson/.local/bin:/home/tjohnson/tools/FastQC:${PATH}"

FASTQ_DIR="${RUN_DIR}/fastq"
READPAIR_INFO_FILE="${RUN_DIR}/config_files/fastq_readpair_by_lane_info.csv"

export _JAVA_OPTIONS="-Xmx${jvmheap}"

## First line is header
let "LINENUM = ${SGE_TASK_ID} + 1"

echo "Parsing sample information for row ${LINENUM}"

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
SAMPLE_FCID_LANE_INDEX=$(echo ${READPAIR_INFO} | awk -F "," '{print $13}')

echo "Parsed info for a file from ${SAMPLE_ID}"


if [[ ${SAMPLE_FCID_LANE_INDEX} == 1 ]] ;  then
	FLOWCELL_LANE=${FLOWCELL_LANE}
elif [[ ${SAMPLE_FCID_LANE_INDEX} == 2 ]] ;  then
	FLOWCELL_LANE=${FLOWCELL_LANE}_${SAMPLE_FCID_LANE_INDEX}
else
	echo "The sample FCID lane index that was passed did not make sense. Exiting."
	exit 1
fi

if [[ ${read_num} == 1 ]]; then
	FASTQ_FILE=${READ1_FASTQ_FILE}
	TEMP_FASTQ_FILE="${SAMPLE_ID}_${FLOWCELL_ID}_L00${FLOWCELL_LANE}_R1.fastq"
elif [[ ${read_num} == 2 ]]; then
	FASTQ_FILE=${READ2_FASTQ_FILE}
	TEMP_FASTQ_FILE="${SAMPLE_ID}_${FLOWCELL_ID}_L00${FLOWCELL_LANE}_R2.fastq"
fi

SOURCE_FASTQ_DIR="${SEQ_RUN_DIR}/${SAMPLE_FASTQ_DIR}"
FASTQ_PATH="${SOURCE_FASTQ_DIR}/${FASTQ_FILE}"

OUTPUT_PATH="${RUN_DIR}/result/fastqc"
TMP_PATH="${OUTPUT_PATH}/tmp"

mkdir -p ${OUTPUT_PATH} ${TMP_PATH}

echo "Output path: ${OUTPUT_PATH}"
echo "Fastq file path: ${FASTQ_PATH}"

fastqc_args="-o ${OUTPUT_PATH} -f fastq"

if [[ ! ${adapters_file} == "none" ]] ; then
	echo "Setting up FastQC for ${SAMPLE_ID} using ${adapters_file} adapters file"
	fastqc_args="${fastqc_args} -a ${adapters_file}"
fi

if [[ -f ${FASTQ_PATH} ]]; then
	echo "Testing if I should use the temp file (${use_temp})"
	if [[ "${use_temp}" == "true" ]] ; then
		temp_file=${TMP_PATH}/${TEMP_FASTQ_FILE}
	
		if [[ "${decompress_file}" == "true" ]]; then
			echo "Uncompressing fastq.bz/gz file for ${SAMPLE_ID} from ${FASTQ_PATH}"
			## temp_file used and removed in trim_and_map
			if [[ -s ${temp_file} ]] ; then
				echo "temp file is already present. Skipping file decompression."
			else
				if [[ "${COMPRESSION_TYPE}" == "bz2" ]] ; then
					pbunzip2 -dk -p${threads} -c ${FASTQ_PATH} > ${temp_file}
				elif [[ "${COMPRESSION_TYPE}" == "gz" ]] ; then
					gunzip -c ${FASTQ_PATH} > ${temp_file}
				else
					echo "Known compression type not found.  Aborting"
					exit
				fi
			fi

			echo "Running FastQC for ${SAMPLE_ID} on temp file created from ${FASTQ_PATH}"
			fastqc "${fastqc_args} ${temp_file}"
			echo "Report should be in ${OUTPUT_PATH}."
		else
			echo "Moving fastq.bz2/gz file for ${SAMPLE_ID} from ${FASTQ_PATH}"
			temp_file=${temp_file}.${COMPRESSION_TYPE}
			mv ${FASTQ_PATH} ${temp_file}
			
			fastqc_args="${fastqc_args} stdin:${TEMP_FASTQ_FILE}"

			if [[ "${COMPRESSION_TYPE}" == "bz2" ]] ; then
				echo "Running FastQC through pbunzip2 to stdin ${SAMPLE_ID} from ${temp_file}"
				pbunzip2 -dk -p${threads} -c ${temp_file} | fastqc ${fastqc_args}
				echo "Report should be in ${OUTPUT_PATH}."
			elif [[ "${COMPRESSION_TYPE}" == "gz" ]] ; then
				echo "Running FastQC through gunzip from stdin ${SAMPLE_ID} from ${temp_file}"
				gunzip -c ${temp_file} | fastqc ${fastqc_args}
				echo "Report should be in ${OUTPUT_PATH}."
			else
				echo "Known compression type not found.  Aborting"
				exit
			fi
		fi
	else
		fastqc_args="${fastqc_args} stdin:${TEMP_FASTQ_FILE}"

		if [[ "${COMPRESSION_TYPE}" == "bz2" ]] ; then
			echo "Running FastQC through pbunzip2 to stdin ${SAMPLE_ID} from ${FASTQ_PATH} labeled as ${TEMP_FASTQ_FILE}"
			pbunzip2 -dk -p${threads} -c ${FASTQ_PATH} | fastqc ${fastqc_args}
			echo "Report should be in ${OUTPUT_PATH}."
		elif [[ "${COMPRESSION_TYPE}" == "gz" ]] ; then
			echo "Running FastQC through gunzip from stdin ${SAMPLE_ID} from ${FASTQ_PATH} labeled as ${TEMP_FASTQ_FILE}"
			gunzip -c ${FASTQ_PATH} | fastqc ${fastqc_args}
			echo "Report should be in ${OUTPUT_PATH}."
		else
			echo "Known compression type not found.  Aborting"
			exit
		fi
	fi
else
	echo "The fastq file for ${SAMPLE_ID} was not found at ${FASTQ_PATH}."
fi

