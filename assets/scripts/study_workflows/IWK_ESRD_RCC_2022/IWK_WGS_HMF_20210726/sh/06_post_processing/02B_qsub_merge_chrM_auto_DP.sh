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

base_logdir="${RUN_DIR}/log/06_post_processing"
logdir="${base_logdir}/extract_chrM_auto_depth"

threads=6
vmem="4G"

mkdir -p ${base_logdir}
mkdir -p ${logdir}

logname="log_02B_merge_chrM_auto_depth"

qsub -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/local/package/r/4.0.2/bin/Rscript ${cg_pipeline_dir}/06_post_processing/02B_merge_chrM_auto_DP_for_germline_genomewide_vcf.R
