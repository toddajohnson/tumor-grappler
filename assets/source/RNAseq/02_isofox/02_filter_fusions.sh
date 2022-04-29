#!/bin/bash

RUN_DIR=${RUN_DIR}

. ${RUN_DIR}/config_files/common_config.sh

#sample_index=${SGE_TASK_ID}
threads=${threads}
jvm_mem=${jvm_mem}

#module use /usr/local/package/modulefiles/
#module load java/8

#export _JAVA_OPTIONS="-Xmx${jvm_mem}"

if [[ "${REFERENCE_GENOME}" == "HG38" ]] || [[ "${REFERENCE_GENOME}" == "GRCh38" ]]; then
	REFERENCE_GENOME="38"
elif [[ "${REFERENCE_GENOME}" == "HG37" ]] || [[ "${REFERENCE_GENOME}" == "GRCh37" ]] || [[ "${REFERENCE_GENOME}" == "HG19" ]] || [[ "${REFERENCE_GENOME}" == "hg19" ]]; then
	REFERENCE_GENOME="37"
fi

#SAMPLE_INFO_FILE="${RUN_DIR}/config_files/sample_info.csv"

## First line is header
#let "LINENUM = ${sample_index} + 1"
#SAMPLE_INFO=$(sed -n "${LINENUM}p" ${SAMPLE_INFO_FILE})
#subject_id=$(echo ${SAMPLE_INFO} | awk -F "," '{print $1}')
#sample_id=$(echo ${SAMPLE_INFO} | awk -F "," '{print $2}')
#WGS_sample_id=$(echo ${SAMPLE_INFO} | awk -F "," '{print $5}')

#bam_dir="${RUN_DIR}/result/STAR-${STAR_VERSION}/${subject_id}"
#bam_path="${bam_dir}/${sample_id}.bam"

#sample_list="${RUN_DIR}/config_files/Isofox_cohort_sample_list.csv"
sample_list="${RUN_DIR}/config_files/Isofox_sample_list.csv"

output_root_path=${RUN_DIR}/result/isofox
output_path=${output_root_path}/${WGS_sample_id}
#mkdir -p ${output_path}

echo "Running isofox cohort generation and fusion filtering for ${study_name}."

#sample_list=$(mktemp /tmp/temp_XXXXX.txt)
#echo ${WGS_sample_id} > ${sample_list}

java -Xmx${jvm_mem} -cp ${ISOFOX_JAR} com.hartwig.hmftools.isofox.cohort.CohortAnalyser \
	-root_data_dir ${output_root_path} \
	-sample_data_file ${sample_list} \
	-use_sample_dir \
	-analyses "FUSION" \
	-fusion_write_filtered -fusion_write_combined \
	-known_fusion_file "${ref_dir}/dbs/knowledgebases/known_fusion_data.38_v3.csv" \
	-threads ${threads} \
	-output_dir ${output_root_path}

#-fusion_write_filtered -fusion_write_combined \
	
#	-analyses "FUSION" \
#	-fusion_gen_cohort -fusion_write_filtered -fusion_write_combined \
#	-output_dir ${output_root_path}
#	-gene_transcripts_dir ${ensembl_cache_dir} \

#-fusion_gen_cohort -fusion_write_filtered -fusion_write_combined \

#rm ${sample_list}

echo "Finished running isofox cohort generation and fusion filtering."
