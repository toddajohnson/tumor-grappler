#!/bin/bash
#$ -pe def_slot 12
#$ -l s_vmem=4G
#$ -N "log_01_filter_resegment_CN_segments"
#$ -o /home/tjohnson/workspace/runs/IWK_WGS_HMF_20210726/log/05_signature_analysis/00_Preparation
#$ -q '!mjobs_rerun.q'
#$ -j y
#$ -cwd
#$ -S /bin/bash

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

. $HOME/.R-4.1.0_setup.sh

#Rscript ${cg_pipeline_dir}/08_signature_analysis/00_Preparation/filter_resegment_copy_number_segments.R 0
Rscript ${cg_pipeline_dir}/08_signature_analysis/00_Preparation/filter_resegment_copy_number_segments.ver_20220131.R 0
