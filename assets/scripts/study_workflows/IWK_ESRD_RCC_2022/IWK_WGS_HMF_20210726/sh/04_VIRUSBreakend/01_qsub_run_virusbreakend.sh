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

base_logdir="${RUN_DIR}/log/03_VIRUSBreakend"
logdir="${base_logdir}/01_run_virusbreakend_to_GPL"

mkdir -p ${base_logdir}
mkdir -p ${logdir}

threads=8
slotct=18
vmem="8G"
#vmem="64G"
jvmheap="31G"

gridss_version="2.12.0"
reference_genome_version="38"

logname="vbe"

#qsub -t 1 -pe def_slot ${slotct} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v threads=${threads},jvmheap=${jvmheap},gridss_version=${gridss_version},reference_genome_version=${reference_genome_version},RUN_DIR=${RUN_DIR},use_archive_run_dir="false" ${cg_pipeline_dir}/03_VIRUSBreakend/01_VIRUSBreakend_to_GPL_dir.sh
qsub -t 2-75 -pe def_slot ${slotct} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v threads=${threads},jvmheap=${jvmheap},gridss_version=${gridss_version},reference_genome_version=${reference_genome_version},RUN_DIR=${RUN_DIR},use_archive_run_dir="false" ${cg_pipeline_dir}/03_VIRUSBreakend/01_VIRUSBreakend_to_GPL_dir.sh

