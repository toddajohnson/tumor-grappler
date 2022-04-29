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

run_dir=${RUN_DIR}

gpl_prefix=$(basename "$PWD")

base_logdir="${run_dir}/log/03_gridss_purple_linx/${gpl_prefix}"
logdir="${base_logdir}/03_gripss_purple_linx"

## Make sure the threads * vmem > 40G
#threads=8
#vmem="8G"
threads=18
vmem="7G"

#PGL_JVM_MEM="31G"
PGL_JVM_MEM="96G"

mkdir -p ${base_logdir}
mkdir -p ${logdir}


#threads=8
#vmem="8G"

#PGL_JVM_MEM="48G"


#germline_target="panel"
germline_target="genomewide"

write_neo_epitopes="true"
somatic_filter_string="filtered"

snvvcf_split_by_tumor="false"

## try forcing maximum purity lower in range
#target_subject_id='IWK042'
#target_sample_id='IWK042_T'
#
#min_ploidy=1.0
#max_ploidy=8
#min_purity=0.08
#max_purity=1.0
#from_purple_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/purple/${target_sample_id}"
#to_purple_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/purple/${target_sample_id}_ploidy_${min_ploidy}_${max_ploidy}_purity_${min_purity}_${max_purity}"
#mv ${from_purple_dir} ${to_purple_dir}
#
#from_linx_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/linx/${target_sample_id}"
#to_linx_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/linx/${target_sample_id}_ploidy_${min_ploidy}_${max_ploidy}_purity_${min_purity}_${max_purity}"
#mv ${from_linx_dir} ${to_linx_dir}
#
#min_ploidy=1.0
#max_ploidy=2.5
#min_purity=0.50
#max_purity=0.55
#logname="GPL_20211216_germline_${germline_target}"
#qsub -t 35 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${run_dir},PGL_THREAD_NUM=${threads},PGL_JVM_MEM=${PGL_JVM_MEM},ignore_bams="true",somatic_filter_string=${somatic_filter_string},snvvcf_split_by_tumor=${snvvcf_split_by_tumor},germline_target=${germline_target},write_neo_epitopes=${write_neo_epitopes},min_ploidy=${min_ploidy},max_ploidy=${max_ploidy},min_purity=${min_purity},max_purity=${max_purity} ${cg_pipeline_dir}/04_gridss_purple_linx/04_gripss_purple_linx.sh
#
#
#target_subject_id='IWK010'
#target_sample_id='IWK010_T'
#
#min_ploidy=1.0
#max_ploidy=8
#min_purity=0.08
#max_purity=1.0
#from_purple_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/purple/${target_sample_id}"
#to_purple_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/purple/${target_sample_id}_ploidy_${min_ploidy}_${max_ploidy}_purity_${min_purity}_${max_purity}"
#mv ${from_purple_dir} ${to_purple_dir}
#
#from_linx_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/linx/${target_sample_id}"
#to_linx_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/linx/${target_sample_id}_ploidy_${min_ploidy}_${max_ploidy}_purity_${min_purity}_${max_purity}"
#mv ${from_linx_dir} ${to_linx_dir}
#
#min_ploidy=1.0
#max_ploidy=2.5
#min_purity=0.35
#max_purity=0.50
#logname="GPL_20211216_germline_${germline_target}"
#qsub -t 8 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${run_dir},PGL_THREAD_NUM=${threads},PGL_JVM_MEM=${PGL_JVM_MEM},ignore_bams="true",somatic_filter_string=${somatic_filter_string},snvvcf_split_by_tumor=${snvvcf_split_by_tumor},germline_target=${germline_target},write_neo_epitopes=${write_neo_epitopes},min_ploidy=${min_ploidy},max_ploidy=${max_ploidy},min_purity=${min_purity},max_purity=${max_purity} ${cg_pipeline_dir}/04_gridss_purple_linx/04_gripss_purple_linx.sh
#
#
# deleted genes from original settings
# are all on chrY and increased purity reduced chr16 gain
target_subject_id='IWK018'
target_sample_id='IWK018_T'
#
#min_ploidy=1.0
#max_ploidy=8
#min_purity=0.08
#max_purity=1.0
#

min_ploidy=1.0
max_ploidy=8
min_purity=0.08
max_purity=1.0

from_purple_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/purple/${target_sample_id}"
to_purple_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/purple/${target_sample_id}_ploidy_${min_ploidy}_${max_ploidy}_purity_${min_purity}_${max_purity}"
mv ${from_purple_dir} ${to_purple_dir}

from_linx_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/linx/${target_sample_id}"
to_linx_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/linx/${target_sample_id}_ploidy_${min_ploidy}_${max_ploidy}_purity_${min_purity}_${max_purity}"
mv ${from_linx_dir} ${to_linx_dir}


min_ploidy=1.0
max_ploidy=8
min_purity=0.30
max_purity=0.31

to_purple_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/purple/${target_sample_id}"
from_purple_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/purple/${target_sample_id}_ploidy_${min_ploidy}_${max_ploidy}_purity_${min_purity}_${max_purity}"
mv ${from_purple_dir} ${to_purple_dir}

to_linx_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/linx/${target_sample_id}"
from_linx_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/linx/${target_sample_id}_ploidy_${min_ploidy}_${max_ploidy}_purity_${min_purity}_${max_purity}"
mv ${from_linx_dir} ${to_linx_dir}

#logname="GPL_20211216_germline_${germline_target}"
#qsub -t 16 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${run_dir},PGL_THREAD_NUM=${threads},PGL_JVM_MEM=${PGL_JVM_MEM},ignore_bams="true",somatic_filter_string=${somatic_filter_string},snvvcf_split_by_tumor=${snvvcf_split_by_tumor},germline_target=${germline_target},write_neo_epitopes=${write_neo_epitopes},min_ploidy=${min_ploidy},max_ploidy=${max_ploidy},min_purity=${min_purity},max_purity=${max_purity} ${cg_pipeline_dir}/04_gridss_purple_linx/04_gripss_purple_linx.sh


#target_subject_id='IWK020'
#target_sample_id='IWK020_T'

#min_ploidy=1.0
#max_ploidy=8
#min_purity=0.08
#max_purity=1.0

#min_ploidy=1.8
#max_ploidy=2.2
#min_purity=0.35
#max_purity=0.45

#min_ploidy=1.0
#max_ploidy=8.0
#min_purity=0.75
#max_purity=0.80

#min_ploidy=1.8
#max_ploidy=2.2
#min_purity=0.77
#max_purity=0.79

#min_ploidy=3.0
#max_ploidy=3.20
#min_purity=0.77
#max_purity=0.79

#min_ploidy=2.8
#max_ploidy=3.2
#min_purity=0.65
#max_purity=0.68
#
#from_purple_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/purple/${target_sample_id}"
#to_purple_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/purple/${target_sample_id}_ploidy_${min_ploidy}_${max_ploidy}_purity_${min_purity}_${max_purity}"
#mv ${from_purple_dir} ${to_purple_dir}
#
#from_linx_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/linx/${target_sample_id}"
#to_linx_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/linx/${target_sample_id}_ploidy_${min_ploidy}_${max_ploidy}_purity_${min_purity}_${max_purity}"
#mv ${from_linx_dir} ${to_linx_dir}
#
# Renamed folder and moved original back into place
#min_ploidy=1.5
#max_ploidy=1.9
#min_purity=0.75
#max_purity=0.78
#logname="GPL_20220414_germline_${germline_target}"
#qsub -t 18 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${run_dir},PGL_THREAD_NUM=${threads},PGL_JVM_MEM=${PGL_JVM_MEM},ignore_bams="true",somatic_filter_string=${somatic_filter_string},snvvcf_split_by_tumor=${snvvcf_split_by_tumor},germline_target=${germline_target},write_neo_epitopes=${write_neo_epitopes},min_ploidy=${min_ploidy},max_ploidy=${max_ploidy},min_purity=${min_purity},max_purity=${max_purity} ${cg_pipeline_dir}/04_gridss_purple_linx/04_gripss_purple_linx.sh


# Changing settings does not really benefit sample
#
#target_subject_id='IWK011'
#target_sample_id='IWK011_T'

#min_ploidy=1.0
#max_ploidy=8
#min_purity=0.08
#max_purity=1.0

#min_ploidy=0.9
#max_ploidy=1.2
#min_purity=0.07
#max_purity=0.12

#min_ploidy=1.8
#max_ploidy=2.5
#min_purity=0.07
#max_purity=0.20

#min_ploidy=1.8
#max_ploidy=2.5
#min_purity=0.07
#max_purity=0.08

#min_ploidy=3.8
#max_ploidy=4.5
#min_purity=0.17
#max_purity=0.25
#
#from_purple_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/purple/${target_sample_id}"
#to_purple_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/purple/${target_sample_id}_ploidy_${min_ploidy}_${max_ploidy}_purity_${min_purity}_${max_purity}"
#mv ${from_purple_dir} ${to_purple_dir}
#
#from_linx_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/linx/${target_sample_id}"
#to_linx_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/linx/${target_sample_id}_ploidy_${min_ploidy}_${max_ploidy}_purity_${min_purity}_${max_purity}"
#mv ${from_linx_dir} ${to_linx_dir}
#
#min_ploidy=1.0
#max_ploidy=8
#min_purity=0.08
#max_purity=1.0
#
#to_purple_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/purple/${target_sample_id}"
#from_purple_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/purple/${target_sample_id}_ploidy_${min_ploidy}_${max_ploidy}_purity_${min_purity}_${max_purity}"
#mv ${from_purple_dir} ${to_purple_dir}
#
#to_linx_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/linx/${target_sample_id}"
#from_linx_dir="${RUN_DIR}/result/${gpl_prefix}/${target_subject_id}/linx/${target_sample_id}_ploidy_${min_ploidy}_${max_ploidy}_purity_${min_purity}_${max_purity}"
#mv ${from_linx_dir} ${to_linx_dir}

#logname="GPL_20220422_germline_${germline_target}"
#qsub -t 9 -pe def_slot ${threads} -l s_vmem=${vmem} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${run_dir},PGL_THREAD_NUM=${threads},PGL_JVM_MEM=${PGL_JVM_MEM},ignore_bams="true",somatic_filter_string=${somatic_filter_string},snvvcf_split_by_tumor=${snvvcf_split_by_tumor},germline_target=${germline_target},write_neo_epitopes=${write_neo_epitopes},min_ploidy=${min_ploidy},max_ploidy=${max_ploidy},min_purity=${min_purity},max_purity=${max_purity} ${cg_pipeline_dir}/04_gridss_purple_linx/04_gripss_purple_linx.sh
