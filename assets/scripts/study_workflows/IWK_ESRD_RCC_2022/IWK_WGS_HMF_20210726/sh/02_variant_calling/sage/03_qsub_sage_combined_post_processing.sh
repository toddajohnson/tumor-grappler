#!/bin/bash

export LANG=en_US.UTF-8

working_dir="$PWD"

dir_name=${working_dir}
dir_base=$(basename ${dir_name})

while [[ ! "${dir_base}" == "sh" ]]; do
	dir_name=$(dirname ${dir_name})
	dir_base=$(basename ${dir_name})
done

dir_name=$(dirname ${dir_name})

. ${dir_name}/config_files/common_config.sh

base_logdir="${RUN_DIR}/log/02_variant_calling"

threads=12
vmem=4G
jvm_mem="12g"

mkdir -p ${base_logdir}


#call_target="germline"
#germline_target="panel"
#logdir="${base_logdir}/sage/${call_target}/03_sage_combined_${call_target}_post_processing"
#mkdir -p ${logdir}
#use_existing_annotated="false"
#refilter_vcf_by_AF_MAF="true"
#logname="${call_target}_${germline_target}"
#qsub -t 1-37 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},jvm_mem=${jvm_mem},call_target=${call_target},germline_target=${germline_target},use_existing_annotated=${use_existing_annotated},refilter_vcf_by_AF_MAF=${refilter_vcf_by_AF_MAF} ${cg_pipeline_dir}/02_variant_calling/sage/03_sage_combined_post_processing.sh

#call_target="germline"
#germline_target="genomewide"
#logdir="${base_logdir}/sage/${call_target}/03_sage_combined_${call_target}_post_processing"
#mkdir -p ${logdir}
#use_existing_annotated="false"
#refilter_vcf_by_AF_MAF="true"
#logname="${call_target}_${germline_target}_refilter"
#qsub -t 1 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},jvm_mem=${jvm_mem},call_target=${call_target},germline_target=${germline_target},use_existing_annotated=${use_existing_annotated},refilter_vcf_by_AF_MAF=${refilter_vcf_by_AF_MAF} ${cg_pipeline_dir}/02_variant_calling/sage/03_sage_combined_post_processing.sh
#qsub -t 1-37 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},jvm_mem=${jvm_mem},call_target=${call_target},germline_target=${germline_target},use_existing_annotated=${use_existing_annotated},refilter_vcf_by_AF_MAF=${refilter_vcf_by_AF_MAF} ${cg_pipeline_dir}/02_variant_calling/sage/03_sage_combined_post_processing.sh


call_target="somatic"
logdir="${base_logdir}/sage/${call_target}/03_sage_combined_${call_target}_post_processing"
mkdir -p ${logdir}
logname="${call_target}"
#sage_run_info_file="${RUN_DIR}/config_files/sage_${call_target}_run_info.tsv"
#use_existing_annotated="false"
#use_existing_annotated="true"
#qsub -t 1-37 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},jvm_mem=${jvm_mem},call_target=${call_target},use_existing_annotated=${use_existing_annotated},sage_run_info_file=${sage_run_info_file} ${cg_pipeline_dir}/02_variant_calling/sage/03_sage_combined_post_processing.sh

## Run calling as joint and then run annotation jointly, then filter that annotated file separately
#sage_run_info_file="${RUN_DIR}/config_files/sage_${call_target}_run_info.tsv"
#use_existing_annotated="false"
#refilter_vcf_by_vaf="false"
#qsub -t 14 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},jvm_mem=${jvm_mem},call_target=${call_target},use_existing_annotated=${use_existing_annotated},sage_run_info_file=${sage_run_info_file},refilter_vcf_by_vaf=${refilter_vcf_by_vaf} ${cg_pipeline_dir}/02_variant_calling/sage/03_sage_combined_post_processing.sh

#sage_run_info_file="${RUN_DIR}/config_files/sage_${call_target}_run_info_single.tsv"
#use_existing_annotated="true"
#refilter_vcf_by_vaf="true"
#qsub -t 14-15 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},jvm_mem=${jvm_mem},call_target=${call_target},use_existing_annotated=${use_existing_annotated},sage_run_info_file=${sage_run_info_file},refilter_vcf_by_vaf=${refilter_vcf_by_vaf} ${cg_pipeline_dir}/02_variant_calling/sage/03_sage_combined_post_processing.sh


somatic_filter_string="weak_filtered"
sage_run_info_file="${RUN_DIR}/config_files/sage_${call_target}_run_info_single.tsv"
use_existing_annotated="true"
refilter_vcf_by_vaf="true"
qsub -t 14-15 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},jvm_mem=${jvm_mem},call_target=${call_target},use_existing_annotated=${use_existing_annotated},sage_run_info_file=${sage_run_info_file},refilter_vcf_by_vaf=${refilter_vcf_by_vaf},somatic_filter_string=${somatic_filter_string} ${cg_pipeline_dir}/02_variant_calling/sage/03_sage_combined_post_processing.sh
