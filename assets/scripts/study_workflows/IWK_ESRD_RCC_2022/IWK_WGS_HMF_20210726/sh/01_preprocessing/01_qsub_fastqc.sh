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
logdir="${base_logdir}/01_fastqc"

mkdir -p ${base_logdir}
mkdir -p ${logdir}

use_temp="true"

#threads=8
threads=1

read_num=1
logname="fastqc_R${read_num}"
#qsub -t 1 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -tc 5 -t 2-35 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -tc 10 -t 36-140 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -tc 10 -t 141-175 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -t 149 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -tc 10 -t 176-340 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -tc 10 -t 341-510 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -tc 10 -t 341-671 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -tc 10 -t 341-451 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -tc 10 -t 452-561 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -tc 10 -t 562-671 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -t 580 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -t 619 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -t 645-646 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -t 651-656 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -t 657-671 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -t 672-966 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -t 672-820 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -t 821-966 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -t 967-1229 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
qsub -t 1146-1229 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh

read_num=2
logname="fastqc_R${read_num}"
#qsub -t 1 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -tc 5 -t 2-35 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -tc 10 -t 36-140 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -tc 10 -t 141-175 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -t 149 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -tc 10 -t 176-340 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -tc 10 -t 341-510 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -tc 10 -t 341-671 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -tc 10 -t 341-451 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -tc 10 -t 452-561 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -tc 10 -t 562-671 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -t 571 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -t 646-648 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -t 653-658 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -t 659-671 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -t 672-966 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -t 672-820 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -t 821-966 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
#qsub -t 967-1229 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
qsub -t 1146-1229 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},READ_NUM=${read_num},use_temp=${use_temp},threads=${threads} ${cg_pipeline_dir}/01_preprocessing/01_fastqc.sh
