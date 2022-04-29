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

#module use /usr/local/package/modulefiles/
#module load R/4.0.2

gpl_prefix=$(basename "$PWD")

base_logdir="${RUN_DIR}/log/05_signature_analysis/CNA_circos"
logdir="${base_logdir}"

threads=4

mkdir -p ${base_logdir}
mkdir -p ${logdir}

logname="log_01_build_CNA_circos_plot"

qsub -t 1 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /bin/bash -v RUN_DIR=${RUN_DIR} ${cg_pipeline_dir}/08_signature_analysis/CNA_circos/01_generate_CNA_circos.sh
qsub -t 2 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /bin/bash -v RUN_DIR=${RUN_DIR} ${cg_pipeline_dir}/08_signature_analysis/CNA_circos/01_generate_CNA_circos.sh
qsub -t 3 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /bin/bash -v RUN_DIR=${RUN_DIR} ${cg_pipeline_dir}/08_signature_analysis/CNA_circos/01_generate_CNA_circos.sh
qsub -t 4 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /bin/bash -v RUN_DIR=${RUN_DIR} ${cg_pipeline_dir}/08_signature_analysis/CNA_circos/01_generate_CNA_circos.sh
qsub -t 5 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /bin/bash -v RUN_DIR=${RUN_DIR} ${cg_pipeline_dir}/08_signature_analysis/CNA_circos/01_generate_CNA_circos.sh
qsub -t 6 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /bin/bash -v RUN_DIR=${RUN_DIR} ${cg_pipeline_dir}/08_signature_analysis/CNA_circos/01_generate_CNA_circos.sh
