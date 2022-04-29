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
logdir="${base_logdir}/07_bam_finalization"

mkdir -p ${base_logdir}
mkdir -p ${logdir}

## If the size of each final bam file is above about 100GB, then it may
## be better to set striping on the bam_final directory and then run this script.
## cd ${RUN_DIR}/result
## mkdir bam_final
## lfs setstripe -c 8 bam_final

sample_type="tumor"
logname="finalize_${sample_type}_bam_files"
#qsub -t 1-10 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/07_finalize_bam_files.sh
#qsub -t 11-21 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/07_finalize_bam_files.sh
#qsub -t 22 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/07_finalize_bam_files.sh
#qsub -t 23 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/07_finalize_bam_files.sh
#qsub -t 24 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/07_finalize_bam_files.sh
#qsub -t 25 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/07_finalize_bam_files.sh
#qsub -t 27 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/07_finalize_bam_files.sh
#qsub -t 28 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/07_finalize_bam_files.sh
#qsub -t 29 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/07_finalize_bam_files.sh
#qsub -t 30 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/07_finalize_bam_files.sh
qsub -t 26 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/07_finalize_bam_files.sh
qsub -t 31-38 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/07_finalize_bam_files.sh

sample_type="normal"
logname="finalize_${sample_type}_bam_files"
#qsub -t 1-10 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/07_finalize_bam_files.sh
#qsub -t 11-14 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/07_finalize_bam_files.sh
#qsub -t 16-21 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/07_finalize_bam_files.sh
#qsub -t 22 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/07_finalize_bam_files.sh
#qsub -t 23 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/07_finalize_bam_files.sh
#qsub -t 24 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/07_finalize_bam_files.sh
#qsub -t 25 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/07_finalize_bam_files.sh
#qsub -t 26 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/07_finalize_bam_files.sh
#qsub -t 28 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/07_finalize_bam_files.sh
#qsub -t 29 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/07_finalize_bam_files.sh
#qsub -t 30 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/07_finalize_bam_files.sh
qsub -t 27 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/07_finalize_bam_files.sh
qsub -t 31-38 -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/07_finalize_bam_files.sh
