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

base_logdir="${RUN_DIR}/log/00_Preparation"
logdir="${base_logdir}"

mkdir -p ${base_logdir}
mkdir -p ${logdir}

threads=18
jvm_mem=31G

excluded_genes_list=${excluded_genes_list}

logname="create_expected_trans_count_file"
qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},jvm_mem=${jvm_mem},excluded_genes_list=${excluded_genes_list} ${RNAseq_cg_pipeline_dir}/00_Preparation/04_create_isofox_expected_trans_counts.sh
