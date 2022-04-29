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

threads=12
vmem=4G

jvm_mem="36g"

base_logdir="${RUN_DIR}/log/02_variant_calling"

logdir="${base_logdir}/sage/germline_PON_mode"
mkdir -p ${logdir}

logname="log_06_germline_PON_mode"
qsub -t 1-37 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},jvm_mem=${jvm_mem},previous="none" ${cg_pipeline_dir}/02_variant_calling/sage/05_sage_germline_PON_caller.sh
