#!/bin/bash

export LANG=en_US.UTF-8

RUN_DIR=${RUN_DIR}
study=${study}
ref_genome=${ref_genome}
mut_type=${mut_type}

. ${RUN_DIR}/config_files/common_config.sh

echo "Activating minoconda environment"
eval  "$(~/miniconda3/bin/conda shell.bash hook)"
source $(conda info --base)/etc/profile.d/conda.sh
conda activate getz_signature_analyzer

dir_name=$(pwd)

echo "Running matrix generation for ${mut_type} in $dir_name"
python ${cg_pipeline_dir}/08_signature_analysis/SigProfiler/01_generate_matrices.py ${RUN_DIR} ${study} ${ref_genome} ${mut_type}
echo "Finished matrix generation"