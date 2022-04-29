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

threads=18
vmem=4G

base_logdir="${RUN_DIR}/log/00.1_QDNAseq"
logdir="${base_logdir}"

mkdir -p ${base_logdir} ${logdir}

logname="log_01_run_QDNAseq"

qsub -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/local/package/r/4.0.2/bin/Rscript ${cg_pipeline_dir}/00.1_QDNAseq/02_run_QDNAseq.R

## combine figures after run
## cd ../../result/QDNAseq_ACE
##	pdfunite qdnaseq_output/*segments.pdf IWK_aegments.pdf
##	pdfunite qdnaseq_output/*calls.pdf IWK_calls.pdf