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

vmem="16G"
jvmheap="12G"

gpl_prefix="GRIDSS-2.12.0"

run_dir="${RUN_DIR}"

base_logdir="${RUN_DIR}/log/05_signature_analysis/maftools"
logdir="${base_logdir}/01_vep_vcf"

mkdir -p ${base_logdir} ${logdir}

#run_type="somatic"
#logname="log_01_vep_${run_type}_vcf"
#qsub -t 1-35 -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v run_dir=${run_dir},run_type=${run_type},gpl_prefix=${gpl_prefix} ${cg_pipeline_dir}/06_post_processing/01_vep_vcf.sh
#qsub -t 8 -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v run_dir=${run_dir},run_type=${run_type},gpl_prefix=${gpl_prefix} ${cg_pipeline_dir}/06_post_processing/01_vep_vcf.sh
#qsub -t 16 -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v run_dir=${run_dir},run_type=${run_type},gpl_prefix=${gpl_prefix} ${cg_pipeline_dir}/06_post_processing/01_vep_vcf.sh
#qsub -t 33 -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v run_dir=${run_dir},run_type=${run_type},gpl_prefix=${gpl_prefix} ${cg_pipeline_dir}/06_post_processing/01_vep_vcf.sh
#qsub -t 14-15 -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v run_dir=${run_dir},run_type=${run_type},gpl_prefix=${gpl_prefix} ${cg_pipeline_dir}/06_post_processing/01_vep_vcf.sh


#run_type="germline"
#logname="log_01_vep_${run_type}_vcf"
#qsub -t 1-35 -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v run_dir=${run_dir},run_type=${run_type},gpl_prefix=${gpl_prefix},jvmheap=${jvmheap} ${cg_pipeline_dir}/06_post_processing/01_vep_vcf.sh

run_type="SVs"
logname="log_01_vep_${run_type}_vcf"
qsub -t 18 -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v run_dir=${run_dir},run_type=${run_type},gpl_prefix=${gpl_prefix},jvmheap=${jvmheap} ${cg_pipeline_dir}/06_post_processing/01_vep_vcf.sh
qsub -t 20 -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v run_dir=${run_dir},run_type=${run_type},gpl_prefix=${gpl_prefix},jvmheap=${jvmheap} ${cg_pipeline_dir}/06_post_processing/01_vep_vcf.sh
