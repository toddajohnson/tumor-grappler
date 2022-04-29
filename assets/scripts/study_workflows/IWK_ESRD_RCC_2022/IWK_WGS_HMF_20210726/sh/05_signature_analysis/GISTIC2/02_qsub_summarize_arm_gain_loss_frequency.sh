#!/bin/bash
#$ -pe def_slot 12
#$ -l s_vmem=4G
#$ -N "log_02_summarize_gain_loss"
#$ -o /home/tjohnson/workspace/runs/IWK_WGS_HMF_20210726/log/05_signature_analysis/GISTIC2
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
#Rscript ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/02_summarize_arm_gain_loss_frequency.R "0.1" "20220201" "AllSamples,ccRCC,non-ccRCC"
#Rscript ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/02_summarize_arm_gain_loss_frequency.R "0.3" "20220201" "AllSamples,ccRCC,non-ccRCC"

#Rscript ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/02_summarize_arm_gain_loss_frequency.R "0.1" "20220202_median" "AllSamples,ccRCC,non-ccRCC"
#Rscript ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/02_summarize_arm_gain_loss_frequency.R "0.3" "20220202_median" "AllSamples,ccRCC,non-ccRCC"
#Rscript ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/02_summarize_arm_gain_loss_frequency.R "0.1" "20220202_none" "AllSamples,ccRCC,non-ccRCC"

Rscript ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/02_summarize_arm_gain_loss_frequency.R "0.1" "20220202_median" "AllSamples,ccRCC,non-ccRCC,ACD-RCC,pRCC,other-RCC"
Rscript ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/02_summarize_arm_gain_loss_frequency.R "0.3" "20220202_median" "AllSamples,ccRCC,non-ccRCC,ACD-RCC,pRCC,other-RCC"