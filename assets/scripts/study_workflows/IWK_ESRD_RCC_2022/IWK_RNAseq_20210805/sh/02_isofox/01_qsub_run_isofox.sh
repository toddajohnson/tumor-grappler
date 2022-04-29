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

threads=12
vmem=4G
#threads=1
#vmem=48G

jvm_mem=25g

base_logdir="${RUN_DIR}/log/02_isofox"
logdir="${base_logdir}"

mkdir -p ${base_logdir}
mkdir -p ${logdir}

logname="log_01_run_isofox"
qsub -t 1-26 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},jvm_mem=${jvm_mem} ${RNAseq_cg_pipeline_dir}/02_isofox/01_run_isofox.sh

#qsub -t 1 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},jvm_mem=${jvm_mem} ${RNAseq_cg_pipeline_dir}/02_isofox/01_run_isofox.sh
#qsub -t 2-22 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},jvm_mem=${jvm_mem} ${RNAseq_cg_pipeline_dir}/02_isofox/01_run_isofox.sh
#qsub -t 23 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},jvm_mem=${jvm_mem} ${RNAseq_cg_pipeline_dir}/02_isofox/01_run_isofox.sh
#qsub -t 24-26 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},jvm_mem=${jvm_mem} ${RNAseq_cg_pipeline_dir}/02_isofox/01_run_isofox.sh
