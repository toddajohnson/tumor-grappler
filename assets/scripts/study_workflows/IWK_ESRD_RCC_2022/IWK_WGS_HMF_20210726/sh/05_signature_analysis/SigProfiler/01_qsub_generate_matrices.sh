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

threads=8

log_dir="${RUN_DIR}/log/05_signature_analysis/SigProfiler"

ref_genome='GRCh38'

min_segment_length=0
mut_type="CN_gt_${min_segment_length}bp"
working_dir="${RUN_DIR}/result/SigProfiler/${study_name}_${mut_type}_signature_analysis"
mkdir -p ${working_dir}

log_name="log_01_generate_${mut_type}_matrix"
qsub -pe def_slot ${threads} -N ${log_name} -o ${log_dir} -wd ${working_dir} -q '!mjobs_rerun.q' -j y -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},study=${study_name}_resegmented,ref_genome=${ref_genome},mut_type=${mut_type} ${cg_pipeline_dir}/08_signature_analysis/SigProfiler/launch_01_generate_matrices.sh

mut_type="mutation"
working_dir="${RUN_DIR}/result/SigProfiler/${study_name}_${mut_type}_signature_analysis"
log_name="log_01_generate_${mut_type}_matrix"
qsub -pe def_slot ${threads} -N ${log_name} -o ${log_dir} -wd ${working_dir} -q '!mjobs_rerun.q' -j y -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},study=${study_name},ref_genome=${ref_genome},mut_type=${mut_type} ${cg_pipeline_dir}/08_signature_analysis/SigProfiler/launch_01_generate_matrices.sh
