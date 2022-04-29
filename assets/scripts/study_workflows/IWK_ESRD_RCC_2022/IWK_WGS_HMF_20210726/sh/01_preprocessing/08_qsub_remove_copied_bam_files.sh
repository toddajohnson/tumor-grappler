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


base_logdir="${RUN_DIR}/log/01_preprocessing"
logdir="${base_logdir}/08_bam_removal"

mkdir -p ${base_logdir}
mkdir -p ${logdir}


sample_type="tumor"
logname="remove_${sample_type}_bam_files"
#qsub -t 1-10 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/08_remove_copied_bam_files.sh
#qsub -t 11-21 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/08_remove_copied_bam_files.sh
#qsub -t 22-23 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/08_remove_copied_bam_files.sh
#qsub -t 24-25 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/08_remove_copied_bam_files.sh
#qsub -t 27-30 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/08_remove_copied_bam_files.sh
qsub -t 26 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/08_remove_copied_bam_files.sh
qsub -t 31-38 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/08_remove_copied_bam_files.sh

sample_type="normal"
logname="remove_${sample_type}_bam_files"
#qsub -t 1-10 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/08_remove_copied_bam_files.sh
#qsub -t 11-14 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/08_remove_copied_bam_files.sh
#qsub -t 16-21 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/08_remove_copied_bam_files.sh
#qsub -t 22-23 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/08_remove_copied_bam_files.sh
#qsub -t 24-25 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/08_remove_copied_bam_files.sh
#qsub -t 26 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/08_remove_copied_bam_files.sh
#qsub -t 28-30 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/08_remove_copied_bam_files.sh
qsub -t 27 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/08_remove_copied_bam_files.sh
qsub -t 31-38 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/08_remove_copied_bam_files.sh
