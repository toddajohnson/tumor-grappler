#!/bin/bash

RUN_DIR=${RUN_DIR}
BMD_THREAD_NUM=${threads}

. ${RUN_DIR}/config_files/common_config.sh

if [[ -z ${load_biobambam2_from_tools} ]] ; then
	load_biobambam2_from_tools="false"
else
	load_biobambam2_from_tools=${load_biobambam2_from_tools}
fi

sample_index=${SGE_TASK_ID}
sample_type=${sample_type}

module use /usr/local/package/modulefiles/

if [[ "${load_biobambam2_from_tools}" == "false" ]] ; then
	echo "Loading biobambam2-2.0.146 from module"
	module load biobambam/2.0.146
else
	echo "Adding biobambam2-2.0.182 to PATH"
	module load gcc/10.2.0
	export LD_LIBRARY_PATH="${HOME}/tools/libdeflate-1.7:${HOME}/tools/snappy-1.1.8/lib64:${HOME}/tools/io_lib-1.14.14/lib:${LD_LIBRARY_PATH}"
	export PATH="${HOME}/tools/biobambam2-2.0.182/bin:${HOME}/.local/bin:${PATH}"
fi

bam_dir="${RUN_DIR}/result/bam"
bam_markdup="${RUN_DIR}/result/bam_markdup"
tmp_bam_dir="${RUN_DIR}/result/bam_markdup/tmp"
SAMPLE_INFO_FILE="${RUN_DIR}/config_files/tumor_normal_pair_info.csv"
READPAIR_INFO_FILE="${RUN_DIR}/config_files/fastq_readpair_info.csv"
READPAIR_BY_LANE_INFO_FILE="${RUN_DIR}/config_files/fastq_readpair_by_lane_info.csv"

mkdir -p ${bam_markdup} ${tmp_bam_dir}

markdup_input_paths=""

## First line is header
let "LINENUM = ${sample_index} + 1"
SAMPLE_INFO=$(sed -n "${LINENUM}p" ${SAMPLE_INFO_FILE})

if [[ ${sample_type} == "tumor" ]]; then
	sample_id=$(echo ${SAMPLE_INFO} | awk -F "," '{print $4}')
elif [[ ${sample_type} == "normal" ]]; then
	sample_id=$(echo ${SAMPLE_INFO} | awk -F "," '{print $6}')
	sample_id_index=$(echo ${SAMPLE_INFO} | awk -F "," '{print $10}')
	
	if [[ $sample_id_index == 0 ]]; then
		echo "Exiting.  There was no normal sample listed in the table (sample.id.N.index==0)."
		exit 1
	fi
	
	if [[ $sample_id_index > 1 ]]; then
		echo "Exiting. This script should only run bammarkduplicates on the first normal sample (sample.id.N.index==1) for a tumor that is listed in the table."
		exit 1
	fi
fi

READPAIR_BY_LANE_INFO_FILE_TMP=$(mktemp /tmp/readpair.info.XXXXXXXXX)
awk -F "," -v sample_id="${sample_id}" \
'$2 == sample_id' ${READPAIR_BY_LANE_INFO_FILE} > ${READPAIR_BY_LANE_INFO_FILE_TMP}
	
while IFS=$',' read -r _ _ _ flowcell_id _ _ _ flowcell_lane _ _ _ _ sample_flowcell_lane_index; do	
	if [[ ${sample_flowcell_lane_index} == 1 ]] ; then
		flowcell_lane=${flowcell_lane}
	elif [[ ${sample_flowcell_lane_index} == 2 ]] ; then
		flowcell_lane=${flowcell_lane}_${sample_flowcell_lane_index}
	else
		echo "The sample_flowcell_lane_index was out of range. Exiting."
		exit 1
	fi
	
	if [[ ${markdup_input_paths} == "" ]]; then
		markdup_input_paths="I=${bam_dir}/${sample_id}_${flowcell_id}_L00${flowcell_lane}.sorted.bam"
	else
		markdup_input_paths="${markdup_input_paths} I=${bam_dir}/${sample_id}_${flowcell_id}_L00${flowcell_lane}.sorted.bam"
	fi
done < ${READPAIR_BY_LANE_INFO_FILE_TMP}

rm ${READPAIR_BY_LANE_INFO_FILE_TMP}

echo "Marking duplicates on ${sample_id}"

markdup_bam_path="${bam_markdup}/${sample_id}.markdup.bam"
metrics_path="${bam_markdup}/${sample_id}.markdup.metrics"
tmp_path="${tmp_bam_dir}/${sample_id}.markdup.tmp"


bammarkduplicates M=${metrics_path} tmpfile=${tmp_path} markthreads=${BMD_THREAD_NUM} rewritebam=1 rewritebamlevel=1 index=1 md5=1 ${markdup_input_paths} O=${markdup_bam_path}

echo "Finished marking duplicates on ${sample_id}"
