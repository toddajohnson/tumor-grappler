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
logdir="${base_logdir}/01_fastqc"

mkdir -p ${base_logdir}
mkdir -p ${logdir}

use_temp="true"

#threads=8
threads=1

read_num=1
logname="fastqc_R${read_num}"
qsub -t 1-27 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${RNAseq_cg_pipeline_dir}/01_preprocessing/01_copy_fastq.sh

read_num=2
logname="fastqc_R${read_num}"
qsub -t 1-27 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${RNAseq_cg_pipeline_dir}/01_preprocessing/01_copy_fastq.sh
