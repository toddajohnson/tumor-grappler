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

gpl_prefix=$(basename "$PWD")

vmem="8G"

base_logdir="${RUN_DIR}/log/03_gridss_purple_linx/${gpl_prefix}"
logdir="${base_logdir}/05_post_pipeline"

mkdir -p ${base_logdir}
mkdir -p ${logdir}

logname="log_05_merge_figures"

qsub -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR} ${cg_pipeline_dir}/04_gridss_purple_linx/05_combine_figures.sh
