#!/bin/bash

export LANG=en_US.UTF-8

. ~/.R-4.1.0_setup.sh

working_dir="$PWD"

dir_name=${working_dir}
dir_base=$(basename ${dir_name})

while [[ ! "${dir_base}" == "sh" ]]; do
	dir_name=$(dirname ${dir_name})
	dir_base=$(basename ${dir_name})
done

dir_name=$(dirname ${dir_name})

. ${dir_name}/config_files/common_config.sh

gpl_prefix=$(basename "$PWD")

logdir="${RUN_DIR}/log/06_post_processing"

vmem="64G"

mkdir -p ${logdir}

logname="log_01_extract_somatic_chrM"

qsub -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S ~/.local/bin/Rscript ${cg_pipeline_dir}/06_post_processing/03_extract_somatic_chrM_from_vcf.R
