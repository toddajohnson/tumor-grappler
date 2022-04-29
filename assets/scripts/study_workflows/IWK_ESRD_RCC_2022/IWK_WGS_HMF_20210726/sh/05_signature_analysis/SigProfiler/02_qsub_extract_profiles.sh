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
vmem="16G"
log_dir="${RUN_DIR}/log/05_signature_analysis/SigProfiler"


working_dir="${RUN_DIR}/result/SigProfiler/${study_name}_mutation_signature_analysis"

mut_type="SBS"
sig_type="SBS96"
log_name="log_02_extract_${mut_type}_profile"
qsub -pe def_slot ${threads} -l s_vmem=${vmem} -N ${log_name} -o ${log_dir} -wd ${working_dir} -q '!mjobs_rerun.q' -j y -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},study=${study_name},mut_type=${mut_type},sig_type=${sig_type} ${cg_pipeline_dir}/08_signature_analysis/SigProfiler/launch_02_extract_profiles.sh

#mut_type="SBS"
#sig_type="SBS288"
#log_name="log_02_extract_${mut_type}_profile"
#qsub -pe def_slot ${threads} -l s_vmem=${vmem} -N ${log_name} -o ${log_dir} -wd ${working_dir} -q '!mjobs_rerun.q' -j y -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},study=${study_name},mut_type=${mut_type},sig_type=${sig_type} ${cg_pipeline_dir}/08_signature_analysis/SigProfiler/launch_02_extract_profiles.sh
#
#mut_type="SBS"
#sig_type="SBS1536"
#log_name="log_02_extract_${mut_type}_profile"
#qsub -pe def_slot ${threads} -l s_vmem=${vmem} -N ${log_name} -o ${log_dir} -wd ${working_dir} -q '!mjobs_rerun.q' -j y -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},study=${study_name},mut_type=${mut_type},sig_type=${sig_type} ${cg_pipeline_dir}/08_signature_analysis/SigProfiler/launch_02_extract_profiles.sh

mut_type="DBS"
sig_type="DBS78"
log_name="log_02_extract_${mut_type}_profile"
qsub -pe def_slot ${threads} -l s_vmem=${vmem} -N ${log_name} -o ${log_dir} -wd ${working_dir} -q '!mjobs_rerun.q' -j y -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},study=${study_name},mut_type=${mut_type},sig_type=${sig_type} ${cg_pipeline_dir}/08_signature_analysis/SigProfiler/launch_02_extract_profiles.sh

mut_type="ID"
sig_type="ID83"
log_name="log_02_extract_${mut_type}_profile"
qsub -pe def_slot ${threads} -l s_vmem=${vmem} -N ${log_name} -o ${log_dir} -wd ${working_dir} -q '!mjobs_rerun.q' -j y -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},study=${study_name},mut_type=${mut_type},sig_type=${sig_type} ${cg_pipeline_dir}/08_signature_analysis/SigProfiler/launch_02_extract_profiles.sh


min_segment_length=0
mut_type="CN_gt_${min_segment_length}bp"
sig_type="ASCAT.CNV"

log_name="log_02_extract_${study_name}_resegmented_${mut_type}_profile"
working_dir="${RUN_DIR}/result/SigProfiler/${study_name}_resegmented_${mut_type}_signature_analysis"
qsub -pe def_slot ${threads} -l s_vmem=${vmem} -N ${log_name} -o ${log_dir} -wd ${working_dir} -q '!mjobs_rerun.q' -j y -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},study=${study_name}_resegmented,mut_type=${mut_type},sig_type=${sig_type} ${cg_pipeline_dir}/08_signature_analysis/SigProfiler/launch_02_extract_profiles.sh
