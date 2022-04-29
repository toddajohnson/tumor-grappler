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

logname="setup_ref"
qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},OVERHANG_BP=125,threads=${threads} ${RNAseq_cg_pipeline_dir}/00_Preparation/01_setup_reference_genome.sh
