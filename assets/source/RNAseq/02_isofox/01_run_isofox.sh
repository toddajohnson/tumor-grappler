#!/bin/bash

RUN_DIR=${RUN_DIR}

. ${RUN_DIR}/config_files/common_config.sh

sample_index=${SGE_TASK_ID}
threads=${threads}
jvm_mem=${jvm_mem}

if [[ "${REFERENCE_GENOME}" == "HG38" ]] || [[ "${REFERENCE_GENOME}" == "GRCh38" ]]; then
	REFERENCE_GENOME="38"
elif [[ "${REFERENCE_GENOME}" == "HG37" ]] || [[ "${REFERENCE_GENOME}" == "GRCh37" ]] || [[ "${REFERENCE_GENOME}" == "HG19" ]] || [[ "${REFERENCE_GENOME}" == "hg19" ]]; then
	REFERENCE_GENOME="37"
fi

SAMPLE_INFO_FILE="${RUN_DIR}/config_files/sample_info.csv"

## First line is header
let "LINENUM = ${sample_index} + 1"
SAMPLE_INFO=$(sed -n "${LINENUM}p" ${SAMPLE_INFO_FILE})
subject_id=$(echo ${SAMPLE_INFO} | awk -F "," '{print $1}')
sample_id=$(echo ${SAMPLE_INFO} | awk -F "," '{print $2}')
WGS_sample_id=$(echo ${SAMPLE_INFO} | awk -F "," '{print $5}')

bam_dir="${RUN_DIR}/result/STAR-${STAR_VERSION}/${subject_id}"
bam_path="${bam_dir}/${sample_id}.bam"

output_path=${RUN_DIR}/result/isofox/${WGS_sample_id}

mkdir -p ${output_path}

echo "Running isofox for ${sample_id}"

if [[ "${ISOFOX_VERSION}" == "1.4" ]] ;  then
	module use /usr/local/package/modulefiles/
	module load java/11
	functions_to_run="TRANSCRIPT_COUNTS;ALT_SPLICE_JUNCTIONS;FUSIONS"
	
	export _JAVA_OPTIONS="-Xmx${jvm_mem}"

	java -jar ${ISOFOX_JAR} \
		-sample ${WGS_sample_id} \
		-functions  ${functions_to_run}\
		-bam_file ${bam_path} \
		-ref_genome ${REFERENCE_GENOME_FA} \
		-ref_genome_version ${REFERENCE_GENOME} \
		-ensembl_data_dir ${ensembl_cache_dir} \
		-output_dir ${output_path} \
		-read_length ${read_length} \
		-exp_counts_file "${ref_dir}/dbs/isofox/read_${read_length}_exp_counts.csv" \
		-exp_gc_ratios_file "${ref_dir}/dbs/isofox/read_${read_length}_exp_gc_ratios.csv" \
		-long_frag_limit 550 \
		-known_fusion_file "${ref_dir}/dbs/knowledgebases/known_fusion_data.38_v3.csv" \
		-threads ${threads} 
else
	module use /usr/local/package/modulefiles/
	module load java/8

	functions_to_run="TRANSCRIPT_COUNTS;NOVEL_LOCATIONS;FUSIONS"
	
	export _JAVA_OPTIONS="-Xmx${jvm_mem}"

	java -jar ${ISOFOX_JAR} \
		-sample ${WGS_sample_id} \
		-functions  ${functions_to_run}\
		-bam_file ${bam_path} \
		-ref_genome ${REFERENCE_GENOME_FA} \
		-ref_genome_version ${REFERENCE_GENOME} \
		-ensembl_data_dir ${ensembl_cache_dir} \
		-excluded_gene_id_file ${excluded_genes_list} \
		-output_dir ${output_path} \
		-apply_calc_frag_lengths \
		-apply_exp_rates \
		-read_length ${read_length} \
		-exp_counts_file "${ref_dir}/dbs/isofox/read_${read_length}_exp_counts.csv" \
		-apply_gc_bias_adjust \
		-exp_gc_ratios_file "${ref_dir}/dbs/isofox/read_${read_length}_exp_gc_ratios.csv" \
		-long_frag_limit 550 \
		-apply_map_qual_adjust \
		-known_fusion_file "${ref_dir}/dbs/knowledgebases/known_fusion_data.38_v3.csv" \
		-threads ${threads} 
fi



# move below to above to turn on longer debug log	
#-log_debug \
	
echo "Finished running isofox for ${sample_id}"
