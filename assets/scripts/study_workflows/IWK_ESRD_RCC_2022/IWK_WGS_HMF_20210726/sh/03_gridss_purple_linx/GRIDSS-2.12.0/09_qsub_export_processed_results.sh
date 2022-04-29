#!/bin/bash

. ../../../config_files/common_config.sh

module use /usr/local/package/modulefiles/
module load R/4.0.2

gpl_prefix=$(basename "$PWD")

base_logdir="${RUN_DIR}/log/03_gridss_purple_linx/${gpl_prefix}"
logdir="${base_logdir}/05_post_pipeline"

vmem="16G"

mkdir -p ${base_logdir}
mkdir -p ${logdir}

logname="log_09_export_process_R_output_ver.20220210"

qsub -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/local/package/r/4.0.2/bin/Rscript ${cg_pipeline_dir}/04_gridss_purple_linx/08_export_processed_results.R
