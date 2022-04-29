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

module use /usr/local/package/modulefiles/
module load R/4.0.2

gpl_prefix=$(basename "$PWD")

base_logdir="${RUN_DIR}/log/05_signature_analysis/SigProfiler"
logdir="${base_logdir}"

threads=6

mkdir -p ${base_logdir}
mkdir -p ${logdir}

logname="log_00_copy_somatic_vcf"

qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/local/package/r/4.0.2/bin/Rscript ${cg_pipeline_dir}/08_signature_analysis/SigProfiler/00_copy_somatic_vcf.R
