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

jvm_mem="36g"

mkdir -p ${base_logdir}

#call_target="germline"
#logdir="${base_logdir}/sage/${call_target}/01_sage_combined_${call_target}_calls"
#mkdir -p ${logdir}
#logname="${call_target}"

#qsub -t 1-37 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},jvm_mem=${jvm_mem},call_target=${call_target} ${cg_pipeline_dir}/02_variant_calling/sage/01_sage_combined_caller.sh

call_target="germline"
logdir="${base_logdir}/sage/${call_target}/01_sage_combined_${call_target}_calls"
germline_target="genomewide"
logname="${call_target}_${germline_target}"
qsub -t 1-37 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},jvm_mem=${jvm_mem},call_target=${call_target},germline_target=${germline_target},previous="none" ${cg_pipeline_dir}/02_variant_calling/sage/01_sage_combined_caller.sh


#call_target="somatic"
#logdir="${base_logdir}/sage/${call_target}/01_sage_combined_${call_target}_calls"
#mkdir -p ${logdir}
#logname="${call_target}"

#qsub -t 1-13 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},jvm_mem=${jvm_mem},call_target=${call_target},previous="none",min_tumor_qual=70 ${cg_pipeline_dir}/02_variant_calling/sage/01_sage_combined_caller.sh
#qsub -t 15-37 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},jvm_mem=${jvm_mem},call_target=${call_target},previous="none",min_tumor_qual=70 ${cg_pipeline_dir}/02_variant_calling/sage/01_sage_combined_caller.sh

#single tumor runs of IWK017 subject tumors
#sage_run_info_file="${RUN_DIR}/config_files/sage_${call_target}_run_info_single.tsv"
#qsub -t 14-15 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},jvm_mem=${jvm_mem},call_target=${call_target},previous="none",min_tumor_qual=70,sage_run_info_file=${sage_run_info_file} ${cg_pipeline_dir}/02_variant_calling/sage/01_sage_combined_caller.sh

#sage_run_info_file="${RUN_DIR}/config_files/sage_${call_target}_run_info.tsv"
#qsub -t 14 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},jvm_mem=${jvm_mem},call_target=${call_target},previous="none",min_tumor_qual=70,sage_run_info_file=${sage_run_info_file} ${cg_pipeline_dir}/02_variant_calling/sage/01_sage_combined_caller.sh
