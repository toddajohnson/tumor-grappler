#!/bin/bash

RUN_DIR=${RUN_DIR}
threads=${threads}
study=${study}
mut_type=${mut_type}
sig_type=${sig_type}

. ${RUN_DIR}/config_files/common_config.sh

export LANG=en_US.UTF-8

echo "Activating minoconda environment"
eval  "$(~/miniconda3/bin/conda shell.bash hook)"
source $(conda info --base)/etc/profile.d/conda.sh
conda activate sigprofiler

echo "Running ${mut_type} signature profiling for ${study}"
python ${cg_pipeline_dir}/08_signature_analysis/SigProfiler/02_extract_profiles.py ${RUN_DIR} ${study} ${mut_type} ${sig_type} "results/${mut_type}/${sig_type}" ${threads}
echo "Finished ${mut_type} signature profiling"