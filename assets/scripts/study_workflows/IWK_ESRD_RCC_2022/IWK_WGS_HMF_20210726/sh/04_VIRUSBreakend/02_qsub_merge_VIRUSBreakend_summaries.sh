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

base_logdir="${RUN_DIR}/log/03_VIRUSBreakend"
logdir="${base_logdir}/01_run_virusbreakend"

mkdir -p ${base_logdir}
mkdir -p ${logdir}

logname="merge_summaries"

slotct=4
vmem=4G

qsub -pe def_slot ${slotct} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/local/package/r/4.0.2/bin/Rscript  ${cg_pipeline_dir}/03_VIRUSBreakend/02_merge_VIRUSBreakend_summaries.R