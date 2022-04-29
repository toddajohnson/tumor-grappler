#!/bin/bash

## first time
eval  "$(~/miniconda3/bin/conda shell.bash hook)"
source $(conda info --base)/etc/profile.d/conda.sh
conda create --name sigprofiler
conda activate sigprofiler
conda install pip

pip install SigProfilerMatrixGenerator
 
python
>>> from SigProfilerMatrixGenerator import install as genInstall
>>> genInstall.install('GRCh38', rsync=False, bash=True)

conda install -c conda-forge mscorefonts
cd /home/tjohnson/miniconda3/envs/sigprofiler/fonts
cp *.ttf ../lib/python3.9/site-packages/matplotlib/mpl-data/fonts/ttf/

cd ~/.cache

conda deactivate



## after previous env creation
eval  "$(~/miniconda3/bin/conda shell.bash hook)"
source $(conda info --base)/etc/profile.d/conda.sh
conda activate sigprofiler
##	which pip
##	~/miniconda3/envs/sigprofiler/bin/pip

## upgraded to version 1.20 on 9/28/2021
pip install --upgrade SigProfilerMatrixGenerator

## some bug fixes made since last install, but release files are old
## installed 11/09/2021
cd ~/tools
mkdir SigProfiler
cd SigProfiler
git clone https://github.com/AlexandrovLab/SigProfilerExtractor.git

pip install SigProfilerExtractor/

## 2022/01/12 versions installed
## pip list 
#SigProfilerExtractor       1.1.3
#SigProfilerMatrixGenerator 1.2.0
#sigProfilerPlotting        1.1.17