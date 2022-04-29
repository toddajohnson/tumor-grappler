#!/bin/bash

RUN_DIR=${RUN_DIR}

sage_run_info_file="${RUN_DIR}/config_files/sage_germline_run_info.tsv"

. ${RUN_DIR}/config_files/common_config.sh
. ${cg_pipeline_dir}/common_config_files/ref_data_config_20210826.sh

SAGE_THREAD_NUM=${threads}
SAGE_JVM_MEM=${jvm_mem}

module use /usr/local/package/modulefiles/
module load java/11
module load R/4.0.2

write_status() { # Before logging initialised
	echo "$(date): $1" 1>&2
}

## First line is header
let "LINENUM = ${SGE_TASK_ID} + 1"

sage_run_info=$(sed -n "${LINENUM}p" ${sage_run_info_file})

subject_id=$(echo "${sage_run_info}" | awk -F "\t" '{print $1}')
tumor_only_flag=$(echo "${sage_run_info}" | awk -F "\t" '{print $2}')
tumor_names=$(echo "${sage_run_info}" | awk -F "\t" '{print $3}')
tumor_bams=$(echo "${sage_run_info}" | awk -F "\t" '{print $4}')
reference_names=$(echo "${sage_run_info}" | awk -F "\t" '{print $5}')
reference_bams=$(echo "${sage_run_info}" | awk -F "\t" '{print $6}')

echo "Setting up SAGE ${call_target} variant calling run for ${subject_id}"

output_dir="${RUN_DIR}/result/mutations/${subject_id}/SAGE-${sage_version}/PON_mode"

log_prefix="${output_dir}/logs/sage_PON_mode_$(date +%Y%m%d_%H%M%S).$HOSTNAME.$$"

if [[ -d ${output_dir} ]]; then
	if [[ "${previous}" == "none" ]] ; then
		rm -R ${output_dir}
	else
		mv ${output_dir} ${output_dir}_${previous}
	fi
fi

echo "Creating output directory at ${output_dir}"
mkdir -p ${output_dir}
mkdir -p ${output_dir}/logs

write_status "Running SAGE for germline PON (tumor only mode) calling of ${subject_id} samples"

. ${cg_pipeline_dir}/common_config_files/germline_PON_config_20211005.sh

# for germline, tumor and reference names/bams are swapped
## from SageCommandBuilder.java 
#            arguments.add("-hard_filter_enabled true")
#                    .add("-soft_filter_enabled false")
#                    .add("-hard_min_tumor_qual 0")
#                    .add("-hard_min_tumor_raw_alt_support 3")
#                    .add("-hard_min_tumor_raw_base_quality 30");
## should also add -high_confidence_bed ${high_confidence_bed} \ but did not see at first glance
## Should not matter except for adding field in INFO, hotspots and panel_bed
java -Xms4G -Xmx${SAGE_JVM_MEM} -cp ${SAGE_JAR} com.hartwig.hmftools.sage.SageApplication \
    -threads ${SAGE_THREAD_NUM} \
    -tumor ${tumor_names} -tumor_bam ${tumor_bams} \
    -assembly ${ref_genome_version} \
    -ref_genome ${ref_genome} \
    -hotspots ${hotspots_vcf} \
    -panel_bed ${coding_panel_bed} \
    -hard_filter_enabled true \
    -soft_filter_enabled false \
    -hard_min_tumor_qual 0 \
    -hard_min_tumor_raw_alt_support 3 \
    -hard_min_tumor_raw_base_quality 30 \
    -out ${output_dir}/${subject_id}.sage.somatic.vcf.gz \
	2>&1 | tee ${log_prefix}.sage.log
