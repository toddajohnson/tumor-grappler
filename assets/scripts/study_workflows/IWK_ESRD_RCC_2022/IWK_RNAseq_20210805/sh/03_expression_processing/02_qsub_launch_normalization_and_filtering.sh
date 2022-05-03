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

base_logdir="${RUN_DIR}/log"
logdir="${base_logdir}/03_expression_data_processing"

threads=6
vmem="5G"

mkdir -p ${base_logdir}
mkdir -p ${logdir}

logname="log_02_normalize_and_filter"
qsub -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR} ${RNAseq_cg_pipeline_dir}/03_expression_data_processing/launch_R_scripts/02_launch_normalize_and_filter.sh
