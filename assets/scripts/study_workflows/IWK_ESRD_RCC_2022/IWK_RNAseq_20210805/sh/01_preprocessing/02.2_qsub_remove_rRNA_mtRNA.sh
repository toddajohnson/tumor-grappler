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

base_logdir="${RUN_DIR}/log/01_preprocessing"
logdir="${base_logdir}/02.2_remove_rRNA_mtRNA"

mkdir -p ${base_logdir}
mkdir -p ${logdir}

use_temp="true"

threads=12
#threads=1
vmem=5G
java_vmem="48g"

logname="remove_rRNA_mtRNA"
#qsub -t 1 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},java_vmem=${java_vmem} ${RNAseq_cg_pipeline_dir}/01_preprocessing/02.2_remove_rRNA_mtRNA.sh
qsub -t 2-27 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},java_vmem=${java_vmem} ${RNAseq_cg_pipeline_dir}/01_preprocessing/02.2_remove_rRNA_mtRNA.sh
