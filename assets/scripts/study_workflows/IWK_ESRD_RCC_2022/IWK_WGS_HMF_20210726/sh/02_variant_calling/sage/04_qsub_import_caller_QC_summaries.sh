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

threads=12
vmem=4G
jvm_mem="12g"

base_logdir="${RUN_DIR}/log/02_variant_calling"
logdir="${base_logdir}/sage"

mkdir -p ${base_logdir}
mkdir -p ${logdir}

logname="log_04_import_called_QC_summaries"

#previous_output=${RUN_DIR}/result_summaries/sage_variant_snpEff_summaries.Rdata
#previous_renamed_output=${RUN_DIR}/result_summaries/sage_variant_snpEff_summaries_20210804.Rdata

#mv $previous_output $previous_renamed_output

qsub -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/local/package/r/4.0.2/bin/Rscript ${cg_pipeline_dir}/02_variant_calling/sage/04_import_caller_QC_summaries.R
