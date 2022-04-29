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

threads=18
vmem=4G

base_logdir="${RUN_DIR}/log/01_preprocessing"
logdir="${base_logdir}/03_STAR_align_fastq"

mkdir -p ${base_logdir}
mkdir -p ${logdir}

logname="log_03_STAR_align"
#qsub -t 1-26 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads} ${RNAseq_cg_pipeline_dir}/01_preprocessing/03B_STAR_align_fastq.sh
qsub -t 26 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads} ${RNAseq_cg_pipeline_dir}/01_preprocessing/03B_STAR_align_fastq.sh
