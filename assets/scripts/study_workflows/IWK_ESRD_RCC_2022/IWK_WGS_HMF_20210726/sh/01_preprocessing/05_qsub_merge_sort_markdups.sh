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
logdir="${base_logdir}/05_merge_sort_markdups"

#threads=1
threads=12
vmem=4G

mkdir -p ${base_logdir}
mkdir -p ${logdir}


logname="bammarkdups_T"
sample_type="tumor"
#qsub -t 1 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/05_merge_sort_markdups.sh
#qsub -t 2 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/05_merge_sort_markdups.sh
#qsub -t 3 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/05_merge_sort_markdups.sh
#qsub -t 4 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/05_merge_sort_markdups.sh
#qsub -t 5 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/05_merge_sort_markdups.sh
#qsub -t 6-7 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/05_merge_sort_markdups.sh
#qsub -t 8-10 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/05_merge_sort_markdups.sh
#qsub -t 11-21 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/05_merge_sort_markdups.sh
#qsub -t 22-30 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/05_merge_sort_markdups.sh
qsub -t 26 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/05_merge_sort_markdups.sh
qsub -t 31-38 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/05_merge_sort_markdups.sh

####
## If there are multiple tumor samples with one normal sample, then care should be taken to run this script against that normal just once.
## Since several instances of this script running against the same sample at the same time will likely result in bad output.
## I have added some logic to the tumor_normal_pair_info.csv file and 05_merge_sort_markdups.sh to exit if this tries to run against
## other than the first instance (row) listing a normal sample.
## Check the log that those tumor's normal index >1 exit
####
logname="bammarkdups_R"
sample_type="normal"
#qsub -t 1 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/05_merge_sort_markdups.sh
#qsub -t 2 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/05_merge_sort_markdups.sh
#qsub -t 3 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/05_merge_sort_markdups.sh
#qsub -t 4 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/05_merge_sort_markdups.sh
#qsub -t 5 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/05_merge_sort_markdups.sh
#qsub -t 6-7 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/05_merge_sort_markdups.sh
#sub -t 8-10 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/05_merge_sort_markdups.sh
#qsub -t 11-14 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/05_merge_sort_markdups.sh
#qsub -t 16-21 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/05_merge_sort_markdups.sh
#qsub -t 22-30 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/05_merge_sort_markdups.sh
qsub -t 27 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/05_merge_sort_markdups.sh
qsub -t 31-38 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},threads=${threads},sample_type=${sample_type} ${cg_pipeline_dir}/01_preprocessing/05_merge_sort_markdups.sh
