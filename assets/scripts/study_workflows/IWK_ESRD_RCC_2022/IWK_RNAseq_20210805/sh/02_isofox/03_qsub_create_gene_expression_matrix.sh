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

logname="log_03_create_gene_expression_matrix"
qsub -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},jvm_mem=${jvm_mem} ${RNAseq_cg_pipeline_dir}/02_isofox/04_create_gene_expression_matrix.sh
