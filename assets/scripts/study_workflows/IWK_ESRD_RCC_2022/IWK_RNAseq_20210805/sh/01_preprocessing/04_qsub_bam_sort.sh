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

bamsort_blockmb=32768
bamsort_threads=9
bamsort_output_threads=18

base_logdir="${RUN_DIR}/log/01_preprocessing"
logdir="${base_logdir}/04_bam_sort"

mkdir -p ${base_logdir}
mkdir -p ${logdir}

logname="log_04_bam_sort"
#qsub -t 1-26 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},BAMSORT_BLOCKMB=${bamsort_blockmb},BAMSORT_THREADS=${bamsort_threads},BAMSORT_OUTPUT_THREADS=${bamsort_output_threads},load_biobambam2_from_tools="true" ${RNAseq_cg_pipeline_dir}/01_preprocessing/04B_bam_sort.sh
qsub -t 26 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},BAMSORT_BLOCKMB=${bamsort_blockmb},BAMSORT_THREADS=${bamsort_threads},BAMSORT_OUTPUT_THREADS=${bamsort_output_threads},load_biobambam2_from_tools="true" ${RNAseq_cg_pipeline_dir}/01_preprocessing/04B_bam_sort.sh
