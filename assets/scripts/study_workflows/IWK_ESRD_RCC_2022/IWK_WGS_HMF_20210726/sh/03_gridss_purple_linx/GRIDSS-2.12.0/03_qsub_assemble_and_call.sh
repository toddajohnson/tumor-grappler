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

run_dir=${RUN_DIR}

gridss_version="2.12.0"

. ${cg_pipeline_dir}/04_gridss_purple_linx/gridss_common.sh

job_nodes=1
def_slot=12
s_vmem=8
thread_num=8
jvm_mem=31g

run_job_indexes="1"
base_log_filename=""
logdir=""
subject_id=""
joint_sample_labels=""
joint_bams_list=""

gridss_version="2.12.0"
use_gridss_suffix="false"
#steps=c("all", "setupreference", "preprocess", "assemble", "call")
steps="all"

base_shdir="${cg_pipeline_dir}/04_gridss_purple_linx"
base_logdir="${run_dir}/log/03_gridss_purple_linx/GRIDSS-${gridss_version}"
GRIDSS_JOINT_sh="01_GRIDSS_joint.sh"

mkdir -p ${base_logdir}

JOINT_CALLING_INFO_FILE="${run_dir}/config_files/joint_tumor_normal_calling_info.tsv"

submit_gridss_script
