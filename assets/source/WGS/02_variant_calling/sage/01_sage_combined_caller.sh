#!/bin/bash

RUN_DIR=${RUN_DIR}

if [[ -z ${previous} ]] ; then
	previous="none"
else
	previous=${previous}
fi

if [[ -z ${germline_target} ]] ; then
	germline_target="panel"
else
	germline_target=${germline_target}
fi

if [[ -z ${sage_run_info_file} ]] ; then
	sage_run_info_file="${RUN_DIR}/config_files/sage_${call_target}_run_info.tsv"
else
	sage_run_info_file=${sage_run_info_file}
fi

. ${RUN_DIR}/config_files/common_config.sh
. ${cg_pipeline_dir}/common_config_files/ref_data_config_20211006.sh

SAGE_THREAD_NUM=${threads}
SAGE_JVM_MEM=${jvm_mem}

call_target=${call_target}

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

if [[ -d "${RUN_DIR}/result/mutations/${subject_id}/SAGE-${sage_version}_previous" ]]; then
	rm -R ${RUN_DIR}/result/mutations/${subject_id}/SAGE-${sage_version}_previous
fi

output_dir="${RUN_DIR}/result/mutations/${subject_id}/SAGE-${sage_version}/${call_target}"

if [[ "${germline_target}" == "genomewide" ]] ; then
	output_dir=${output_dir}_${germline_target}
fi

log_prefix="${output_dir}/logs/sage_${call_target}_$(date +%Y%m%d_%H%M%S).$HOSTNAME.$$"

if [[ -d "${RUN_DIR}/result/mutations/${subject_id}/SAGE-" ]]; then
	rm -R "${RUN_DIR}/result/mutations/${subject_id}/SAGE-"
fi

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

write_status "Running SAGE for ${call_target} calling of ${subject_id} samples"

if [[ "${tumor_only_flag}" == "FALSE" ]]; then
	if [[ "${call_target}" == "germline" ]]; then
		. ${cg_pipeline_dir}/common_config_files/germline_config_20210828.sh
	
		if [[ "${germline_target}" == "panel" ]] ; then
			write_status "Running SAGE in tumor/reference mode for gene panel germline variant calling for ${subject_id}"
					
			java -Xms4G -Xmx${SAGE_JVM_MEM} -cp ${SAGE_JAR} com.hartwig.hmftools.sage.SageApplication \
				-threads ${SAGE_THREAD_NUM} \
				-reference ${reference_names} -reference_bam ${reference_bams} \
				-tumor ${tumor_names} -tumor_bam ${tumor_bams} \
				-coverage_bed ${coverage_bed} \
				-panel_only \
				-hotspot_min_tumor_qual ${hotspot_min_tumor_qual} \
				-panel_min_tumor_qual ${panel_min_tumor_qual} \
				-hotspot_max_germline_vaf  ${hotspot_max_germline_vaf} \
				-hotspot_max_germline_rel_raw_base_qual ${hotspot_max_germline_rel_raw_base_qual} \
				-panel_max_germline_vaf ${panel_max_germline_vaf} \
				-panel_max_germline_rel_raw_base_qual ${panel_max_germline_rel_raw_base_qual} \
				-mnv_filter_enabled false \
				-assembly ${ref_genome_version} \
				-ref_genome ${ref_genome} \
				-hotspots ${hotspots_vcf} \
				-panel_bed ${coding_panel_bed} \
				-high_confidence_bed ${high_confidence_bed} \
				-out ${output_dir}/${subject_id}.${call_target}.sage.vcf.gz \
				2>&1 | tee ${log_prefix}.sage.log
		elif [[ "${germline_target}" == "genomewide" ]] ; then
			write_status "Running SAGE in tumor/reference mode for genomewide germline variant calling for ${subject_id}"
			output_vcf=${output_dir}/${subject_id}.${call_target}.sage.vcf.gz
			write_status "Output to ${output_vcf}; hotspot_min_tumor_qual = ${hotspot_min_tumor_qual}"
			
			java -Xms4G -Xmx${SAGE_JVM_MEM} -cp ${SAGE_JAR} com.hartwig.hmftools.sage.SageApplication \
				-threads ${SAGE_THREAD_NUM} \
				-reference ${reference_names} -reference_bam ${reference_bams} \
				-tumor ${tumor_names} -tumor_bam ${tumor_bams} \
				-mnv_filter_enabled false \
				-assembly ${ref_genome_version} \
				-ref_genome ${ref_genome} \
				-hotspots ${hotspots_vcf} \
				-panel_bed ${coding_panel_bed} \
				-high_confidence_bed ${high_confidence_bed} \
				-hotspot_min_tumor_qual ${hotspot_min_tumor_qual} \
				-panel_min_tumor_qual ${panel_min_tumor_qual} \
				-high_confidence_min_tumor_qual ${high_confidence_min_tumor_qual} \
				-low_confidence_min_tumor_qual ${low_confidence_min_tumor_qual} \
				-hotspot_min_tumor_vaf ${hotspot_min_tumor_vaf} \
				-panel_min_tumor_vaf ${panel_min_tumor_vaf} \
				-high_confidence_min_tumor_vaf ${high_confidence_min_tumor_vaf} \
				-low_confidence_min_tumor_vaf ${low_confidence_min_tumor_vaf} \
				-hotspot_min_germline_depth ${hotspot_min_germline_depth} \
				-panel_min_germline_depth ${panel_min_germline_depth} \
				-high_confidence_min_germline_depth ${high_confidence_min_germline_depth} \
				-low_confidence_min_germline_depth ${low_confidence_min_germline_depth} \
				-hotspot_min_germline_depth_allosome ${hotspot_min_germline_depth_allosome} \
				-panel_min_germline_depth_allosome ${panel_min_germline_depth_allosome} \
				-high_confidence_min_germline_depth_allosome ${high_confidence_min_germline_depth_allosome} \
				-low_confidence_min_germline_depth_allosome ${low_confidence_min_germline_depth_allosome} \
				-hotspot_max_germline_vaf ${hotspot_max_germline_vaf} \
				-panel_max_germline_vaf ${panel_max_germline_vaf} \
				-high_confidence_max_germline_vaf ${high_confidence_max_germline_vaf} \
				-low_confidence_max_germline_vaf ${low_confidence_max_germline_vaf} \
				-hotspot_max_germline_rel_raw_base_qual ${hotspot_max_germline_rel_raw_base_qual} \
				-panel_max_germline_rel_raw_base_qual ${panel_max_germline_rel_raw_base_qual} \
				-high_confidence_max_germline_rel_raw_base_qual ${high_confidence_max_germline_rel_raw_base_qual} \
				-low_confidence_max_germline_rel_raw_base_qual ${low_confidence_max_germline_rel_raw_base_qual} \
				-out ${output_vcf} \
				2>&1 | tee ${log_prefix}.sage.log
		fi
	elif [[ "${call_target}" == "somatic" ]]; then
		write_status "Running SAGE in tumor/reference mode for somatic variant calling for ${subject_id}"
		
		. ${cg_pipeline_dir}/common_config_files/somatic_config_20210828.sh
		
		java -Xms4G -Xmx${SAGE_JVM_MEM} -cp ${SAGE_JAR} com.hartwig.hmftools.sage.SageApplication \
			-threads ${SAGE_THREAD_NUM} \
			-reference ${reference_names} -reference_bam ${reference_bams} \
			-tumor ${tumor_names} -tumor_bam ${tumor_bams} \
			-assembly ${ref_genome_version} \
			-ref_genome ${ref_genome} \
			-hotspots ${hotspots_vcf} \
			-panel_bed ${coding_panel_bed} \
			-high_confidence_bed ${high_confidence_bed} \
			-hotspot_min_tumor_qual ${hotspot_min_tumor_qual} \
			-panel_min_tumor_qual ${panel_min_tumor_qual} \
			-high_confidence_min_tumor_qual ${high_confidence_min_tumor_qual} \
			-low_confidence_min_tumor_qual ${low_confidence_min_tumor_qual} \
			-hotspot_min_tumor_vaf ${hotspot_min_tumor_vaf} \
			-panel_min_tumor_vaf ${panel_min_tumor_vaf} \
			-high_confidence_min_tumor_vaf ${high_confidence_min_tumor_vaf} \
			-low_confidence_min_tumor_vaf ${low_confidence_min_tumor_vaf} \
			-hotspot_min_germline_depth ${hotspot_min_germline_depth} \
			-panel_min_germline_depth ${panel_min_germline_depth} \
			-high_confidence_min_germline_depth ${high_confidence_min_germline_depth} \
			-low_confidence_min_germline_depth ${low_confidence_min_germline_depth} \
			-hotspot_min_germline_depth_allosome ${hotspot_min_germline_depth_allosome} \
			-panel_min_germline_depth_allosome ${panel_min_germline_depth_allosome} \
			-high_confidence_min_germline_depth_allosome ${high_confidence_min_germline_depth_allosome} \
			-low_confidence_min_germline_depth_allosome ${low_confidence_min_germline_depth_allosome} \
			-hotspot_max_germline_vaf ${hotspot_max_germline_vaf} \
			-panel_max_germline_vaf ${panel_max_germline_vaf} \
			-high_confidence_max_germline_vaf ${high_confidence_max_germline_vaf} \
			-low_confidence_max_germline_vaf ${low_confidence_max_germline_vaf} \
			-hotspot_max_germline_rel_raw_base_qual ${hotspot_max_germline_rel_raw_base_qual} \
			-panel_max_germline_rel_raw_base_qual ${panel_max_germline_rel_raw_base_qual} \
			-high_confidence_max_germline_rel_raw_base_qual ${high_confidence_max_germline_rel_raw_base_qual} \
			-low_confidence_max_germline_rel_raw_base_qual ${low_confidence_max_germline_rel_raw_base_qual} \
			-out ${output_dir}/${subject_id}.${call_target}.sage.vcf.gz \
			2>&1 | tee ${log_prefix}.sage.log
	fi
elif [[ "${tumor_only_flag}" == "TRUE" ]]; then
	if [[ "${call_target}" == "somatic" ]]; then
		write_status "Running SAGE in tumor only mode for somatic variant calling for ${subject_id}"
		
		. ${cg_pipeline_dir}/common_config_files/somatic_tumor_only_config_20211005.sh
		
		java -Xms4G -Xmx${SAGE_JVM_MEM} -cp ${SAGE_JAR} com.hartwig.hmftools.sage.SageApplication \
		    -threads ${SAGE_THREAD_NUM} \
		    -tumor ${tumor_names} -tumor_bam ${tumor_bams} \
		    -assembly ${ref_genome_version} \
		    -ref_genome ${ref_genome} \
		    -hotspots ${hotspots_vcf} \
		    -panel_bed ${coding_panel_bed} \
		    -high_confidence_bed ${high_confidence_bed} \
		    -hotspot_min_tumor_qual ${hotspot_min_tumor_qual} \
			-panel_min_tumor_qual ${panel_min_tumor_qual} \
			-high_confidence_min_tumor_qual ${high_confidence_min_tumor_qual} \
			-low_confidence_min_tumor_qual ${low_confidence_min_tumor_qual} \
			-hotspot_min_tumor_vaf ${hotspot_min_tumor_vaf} \
			-panel_min_tumor_vaf ${panel_min_tumor_vaf} \
			-high_confidence_min_tumor_vaf ${high_confidence_min_tumor_vaf} \
			-low_confidence_min_tumor_vaf ${low_confidence_min_tumor_vaf} \
		    -out ${output_dir}/${subject_id}.${call_target}.sage.vcf.gz \
			2>&1 | tee ${log_prefix}.sage.log
	else
		write_status "Tumor only mode should only be for somatic variant calling"
	fi
fi
