#!/bin/bash

RUN_DIR=${RUN_DIR}
JVM_MEM=${jvm_mem}
call_target=${call_target}
threads=${threads}

if [[ -z ${use_existing_annotated} ]] ; then
	use_existing_annotated="false"
else
	use_existing_annotated=${use_existing_annotated}
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

if [[ -z ${refilter_vcf_by_vaf} ]] ; then
	refilter_vcf_by_vaf="false"
else
	refilter_vcf_by_vaf=${refilter_vcf_by_vaf}
fi

if [[ -z ${refilter_vcf_by_AF_MAF} ]] ; then
	refilter_vcf_by_AF_MAF="false"
else
	refilter_vcf_by_AF_MAF=${refilter_vcf_by_AF_MAF}
fi

if [[ -z ${refilter_by_chrM} ]] ; then
	refilter_by_chrM="false"
else
	refilter_by_chrM=${refilter_by_chrM}
fi

if [[ -z ${somatic_filter_string} ]] ; then
	somatic_filter_string="filtered"
else
	somatic_filter_string=${somatic_filter_string}
fi


. ${RUN_DIR}/config_files/common_config.sh
. ${cg_pipeline_dir}/common_config_files/ref_data_config_20211006.sh

## samtools and bcftools 1.12, cutadapt are in /home/tjohnson/.local
export PATH="/home/tjohnson/tools/annovar:/home/tjohnson/.local/bin:${PATH}"
module use /usr/local/package/modulefiles/
module load java/12
module load R/4.0.2

vcf_af_calc_R="${cg_pipeline_dir}/common_scripts/vcf_calculate_AF_from_AN_AC.R"

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

base_output_dir="${RUN_DIR}/result/mutations/${subject_id}/SAGE-${sage_version}/${call_target}"

write_status "Running SAGE for ${call_target} calling of ${subject_id} samples"

write_status "Loaded ${call_target} run information."
write_status "Tumor names: ${tumor_names}.  Reference names: ${reference_names}"
write_status "Tumor names: ${tumor_bams}."
write_status "Reference names: ${reference_bams}"

output_dir="${base_output_dir}"

if [[ "${germline_target}" == "genomewide" ]] ; then
	output_dir=${output_dir}_${germline_target}
fi

tmp_output_dir="${output_dir}/tmp"
mkdir -p "${tmp_output_dir}"


## probably need to launch this after any contained variables are instantiated
. ${cg_pipeline_dir}/common_scripts/vcf_processing.sh

cd ${output_dir}

if [[ -d ${output_dir}/${call_target} ]] ; then
	write_status "Removing old misplaced snpEff directory at ${output_dir}/${call_target}"
	rm -R ${output_dir}/${call_target}
fi

if [[ "${call_target}" == "germline" ]]; then
	. ${cg_pipeline_dir}/common_config_files/germline_config_20210828.sh

	if [[ "${use_existing_annotated}" == "false" ]] || [[ ! -s "${output_dir}/${subject_id}.${call_target}.sage.filtered.vcf.gz" ]] ; then	
		remove_previous_filtered_output
	
		write_status "Running post-processing of SAGE in tumor/reference mode germline variant calls for ${subject_id}"
	
		sample_names="${tumor_names},${reference_names}"
	
		write_status "Filtering and annotating mappability"
	
		bcftools view -f PASS -s ${sample_names} -Ob -o "${tmp_output_dir}/sage.bcf" ${output_dir}/${subject_id}.${call_target}.sage.vcf.gz
		bcftools index "${tmp_output_dir}/sage.bcf"
	
		annotate_mappability \
			"${tmp_output_dir}/sage.bcf" \
			"${tmp_output_dir}/sage.mappability.annotated.bcf"
		annotate_clinvar \
			"${tmp_output_dir}/sage.mappability.annotated.bcf" \
			"${tmp_output_dir}/sage.mappability_clinvar.annotated.bcf"
		annotate_blacklisted \
			"${tmp_output_dir}/sage.mappability_clinvar.annotated.bcf" \
			"${tmp_output_dir}/sage.mappability_clinvar_blacklist.annotated.vcf.gz"
		# how to deal with MNVs?
		annotate_ALFA_dbSNP \
			"${tmp_output_dir}/sage.mappability_clinvar_blacklist.annotated.vcf.gz" \
			"${tmp_output_dir}/sage.mappability_clinvar_blacklist_ALFA_dbSNP.annotated.vcf.gz"	
		filter_PASS_variants \
			"${tmp_output_dir}/sage.mappability_clinvar_blacklist_ALFA_dbSNP.annotated.vcf.gz" \
			"${output_dir}/${subject_id}.${call_target}.sage.filtered.vcf.gz" \
			"filtered"
	fi
	
	if [[ "${refilter_vcf_by_AF_MAF}" == "true" ]] ; then
		if [[ -s ${output_dir}/${subject_id}.${call_target}.sage.AF_MAF_filtered.vcf.gz ]] ; then
			rm ${output_dir}/${subject_id}.${call_target}.sage.AF_MAF_filtered.vcf.gz ${output_dir}/${subject_id}.${call_target}.sage.AF_MAF_filtered.vcf.gz.tbi
		fi
		filter_on_AF_MAF \
			${output_dir}/${subject_id}.${call_target}.sage.filtered.vcf.gz \
			${output_dir}/${subject_id}.${call_target}.sage.with_AF_MAF_filter.vcf.gz \
			${output_dir}/${subject_id}.${call_target}.sage.AF_MAF_filtered.vcf.gz \
			"AF_MAF_filtered"
	fi
	
	if [[ "${refilter_by_chrM}" == "true" ]] ; then
		if [[ -s ${output_dir}/${subject_id}.${call_target}.sage.chrM_filtered.vcf.gz ]] ; then
			rm ${output_dir}/${subject_id}.${call_target}.sage.chrM_filtered.vcf.gz ${output_dir}/${subject_id}.${call_target}.sage.chrM_filtered.vcf.gz.tbi
		fi
		
		filter_on_chrM \
			${output_dir}/${subject_id}.${call_target}.sage.filtered.vcf.gz \
			${output_dir}/${subject_id}.${call_target}.sage.chrM_filtered.vcf.gz \
			"chrM_filtered"
	fi
	
	rm ${tmp_output_dir}/*
elif [[ "${call_target}" == "somatic" ]]; then
	write_status "Running post-processing of SAGE in tumor/reference mode somatic variant calls for ${subject_id}"

	. ${cg_pipeline_dir}/common_config_files/somatic_config_20210828.sh

	if [[ "${tumor_only_flag}" == "TRUE" ]]; then
		sample_names="${tumor_names}"
	elif [[ "${tumor_only_flag}" == "FALSE" ]]; then
		sample_names="${tumor_names},${reference_names}"
	fi

	if [[ "${use_existing_annotated}" == "false" ]] || [[ ! -s "${output_dir}/${subject_id}.${call_target}.sage.PON_filtered.hg38_multianno.vcf.gz" ]] ; then
		write_status "Selecting unannotated vcf to include tumor and reference samples"
		bcftools view -s ${sample_names} -Ob -o "${tmp_output_dir}/sage.bcf" ${output_dir}/${subject_id}.${call_target}.sage.vcf.gz
		bcftools index "${tmp_output_dir}/sage.bcf"

		remove_previous_annotation_output

		write_status "Annotating somatic calls using PON: ${germlinePon}"
	
		annotate_PON_mappability \
			"${tmp_output_dir}/sage.bcf" \
			"${tmp_output_dir}/sage.PON_mappability.annotated.bcf"
		annotate_clinvar \
			"${tmp_output_dir}/sage.PON_mappability.annotated.bcf" \
			"${tmp_output_dir}/sage.PON_mappability_clinvar.annotated.bcf"
		remove_PON_filtered \
			"${tmp_output_dir}/sage.PON_mappability_clinvar.annotated.bcf" \
			"${output_dir}/${subject_id}.${call_target}.sage.PON_filtered.vcf.gz"
		annotate_ALFA_dbSNP \
			"${output_dir}/${subject_id}.${call_target}.sage.PON_filtered.vcf.gz" \
			"${tmp_output_dir}/sage.PON_filtered.ALFA_dbSNP.annotated.vcf.gz"
		annotate_with_annovar \
			"${tmp_output_dir}/sage.PON_filtered.ALFA_dbSNP.annotated.vcf.gz" \
			"${output_dir}/${subject_id}.${call_target}.sage.PON_filtered"
	fi
	
	if [[ "${refilter_vcf_by_vaf}" == "false" ]] ; then
		remove_previous_filtered_output

		filter_PASS_variants \
			${output_dir}/${subject_id}.${call_target}.sage.PON_filtered.hg38_multianno.vcf.gz \
			${output_dir}/${subject_id}.${call_target}.sage.weak_filtered.vcf.gz \
			"weak_filtered"

		write_status "Filtering for new relaxed cutoffs"
		
		filter_on_min_tumor_qual \
			${output_dir}/${subject_id}.${call_target}.sage.weak_filtered.vcf.gz \
			${output_dir}/${subject_id}.${call_target}.sage.filtered.vcf.gz \
			"min_tumor_qual" \
			"filtered"

		filter_on_AF_MAF \
			${output_dir}/${subject_id}.${call_target}.sage.filtered.vcf.gz \
			${output_dir}/${subject_id}.${call_target}.sage.with_AF_MAF_filter.vcf.gz \
			${output_dir}/${subject_id}.${call_target}.sage.AF_MAF_filtered.vcf.gz \
			"AF_MAF_filtered"
	elif [[ "${refilter_vcf_by_vaf}" == "true" ]] ; then
		## This should only be used when a subject with multiple tumors needs to have the annotated vcf split by tumor
		## so each refiltering run should have only one sample in tumor_names
		write_status "Selecting annotated-filtered vcf to include current tumor and reference sample"
		bcftools view -s ${sample_names} -Ob -o "${tmp_output_dir}/sage.${tumor_names}.bcf" ${output_dir}/${subject_id}.${call_target}.sage.${somatic_filter_string}.vcf.gz
		bcftools index "${tmp_output_dir}/sage.${tumor_names}.bcf"
	
		filter_on_tumor_vaf \
			"${tmp_output_dir}/sage.${tumor_names}.bcf" \
			${output_dir}/${tumor_names}.${call_target}.sage.${somatic_filter_string}.vcf.gz \
			"min_tumor_vaf" \
			"${tumor_names}_${somatic_filter_string}"
		
		rm ${tmp_output_dir}/sage.${tumor_names}.bcf ${tmp_output_dir}/sage.${tumor_names}.bcf.csi
	fi

	rm ${tmp_output_dir}/*
fi
write_status "Finished processing ${subject_id} ${call_target} variant calls"