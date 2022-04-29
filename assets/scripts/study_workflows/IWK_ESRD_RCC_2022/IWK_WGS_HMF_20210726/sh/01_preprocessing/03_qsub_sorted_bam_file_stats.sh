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

vmem=4G
threads=4

base_logdir="${RUN_DIR}/log/01_preprocessing"
logdir="${base_logdir}/03_bam_file_stats"

mkdir -p ${base_logdir}
mkdir -p ${logdir}

logname="sorted_bam_file_stats"


#qsub -t 1-35 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v threads=${threads},RUN_DIR=${RUN_DIR} ${cg_pipeline_dir}/01_preprocessing/03_sorted_bam_file_stats.sh
#qsub -t 36-70 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v threads=${threads},RUN_DIR=${RUN_DIR} ${cg_pipeline_dir}/01_preprocessing/03_sorted_bam_file_stats.sh
#qsub -t 71-105 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v threads=${threads},RUN_DIR=${RUN_DIR} ${cg_pipeline_dir}/01_preprocessing/03_sorted_bam_file_stats.sh
#qsub -t 106-140 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v threads=${threads},RUN_DIR=${RUN_DIR} ${cg_pipeline_dir}/01_preprocessing/03_sorted_bam_file_stats.sh
#qsub -tc 10 -t 141-175 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v threads=${threads},RUN_DIR=${RUN_DIR} ${cg_pipeline_dir}/01_preprocessing/03_sorted_bam_file_stats.sh
#qsub -tc 10 -t 176-245 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v threads=${threads},RUN_DIR=${RUN_DIR} ${cg_pipeline_dir}/01_preprocessing/03_sorted_bam_file_stats.sh
#qsub -tc 10 -t 246-340 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v threads=${threads},RUN_DIR=${RUN_DIR} ${cg_pipeline_dir}/01_preprocessing/03_sorted_bam_file_stats.sh
#qsub -tc 10 -t 341-570 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v threads=${threads},RUN_DIR=${RUN_DIR} ${cg_pipeline_dir}/01_preprocessing/03_sorted_bam_file_stats.sh
#qsub -tc 10 -t 571-647 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v threads=${threads},RUN_DIR=${RUN_DIR} ${cg_pipeline_dir}/01_preprocessing/03_sorted_bam_file_stats.sh
#qsub -t 648-671 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v threads=${threads},RUN_DIR=${RUN_DIR} ${cg_pipeline_dir}/01_preprocessing/03_sorted_bam_file_stats.sh
#qsub -tc 20 -t 672-820 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v threads=${threads},RUN_DIR=${RUN_DIR} ${cg_pipeline_dir}/01_preprocessing/03_sorted_bam_file_stats.sh
#qsub -t 821-853 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v threads=${threads},RUN_DIR=${RUN_DIR} ${cg_pipeline_dir}/01_preprocessing/03_sorted_bam_file_stats.sh
#qsub -t 854-893 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v threads=${threads},RUN_DIR=${RUN_DIR} ${cg_pipeline_dir}/01_preprocessing/03_sorted_bam_file_stats.sh
#qsub -t 894-908 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v threads=${threads},RUN_DIR=${RUN_DIR} ${cg_pipeline_dir}/01_preprocessing/03_sorted_bam_file_stats.sh
#qsub -t 909-966 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v threads=${threads},RUN_DIR=${RUN_DIR} ${cg_pipeline_dir}/01_preprocessing/03_sorted_bam_file_stats.sh
##qsub -t 967-1229 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v threads=${threads},RUN_DIR=${RUN_DIR} ${cg_pipeline_dir}/01_preprocessing/03_sorted_bam_file_stats.sh
qsub -t 833 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v threads=${threads},RUN_DIR=${RUN_DIR} ${cg_pipeline_dir}/01_preprocessing/03_sorted_bam_file_stats.sh
qsub -t 841 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v threads=${threads},RUN_DIR=${RUN_DIR} ${cg_pipeline_dir}/01_preprocessing/03_sorted_bam_file_stats.sh
