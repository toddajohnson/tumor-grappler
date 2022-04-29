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


base_logdir="${RUN_DIR}/log/00.1_QDNAseq"
logdir="${base_logdir}/01_downsample_bams"

vmem=5G

mkdir -p ${base_logdir}
mkdir -p ${logdir}

logname="downsample_bam"

qsub -tc 10 -t 1-75 -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR} ${cg_pipeline_dir}/00.1_QDNAseq/01_downsample_bam.sh
