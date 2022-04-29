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

base_logdir="${RUN_DIR}/log/02_variant_calling"

threads=12
vmem=4G
jvm_mem="12g"

mkdir -p ${base_logdir}


call_target="somatic"
logdir="${base_logdir}/sage/${call_target}/03_sage_combined_${call_target}_post_processing"
mkdir -p ${logdir}
logname="${call_target}"

qsub -t 1-37 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},jvm_mem=${jvm_mem},call_target=${call_target},use_existing_annotated=${use_existing_annotated},sage_run_info_file=${sage_run_info_file} ${cg_pipeline_dir}/02_variant_calling/sage/03_sage_combined_post_processing.sh


