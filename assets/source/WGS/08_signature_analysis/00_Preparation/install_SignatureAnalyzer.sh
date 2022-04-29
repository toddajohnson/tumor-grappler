#!/bin/bash

## first time
eval  "$(~/miniconda3/bin/conda shell.bash hook)"
source $(conda info --base)/etc/profile.d/conda.sh
conda create --name getz_signature_analyzer python=3.9
conda activate getz_signature_analyzer
#conda install pip

pip install twobitreader
conda install pandas



conda update -n base -c defaults conda

cd ~/tools/SignatureAnalyzer/getz/getzlab-SignatureAnalyzer
pip3 install -e .

pip install torch===1.10.0 torchvision===0.5.0 -f https://download.pytorch.org/whl/torch_stable.html

module use /usr/local/package/modulefiles/
module load gcc/9.3.0

pip install signatureanalyzer

conda deactivate