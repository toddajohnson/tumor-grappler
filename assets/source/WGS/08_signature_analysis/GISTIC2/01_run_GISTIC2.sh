#!/bin/bash

RUN_DIR=${RUN_DIR}
output_dir=${output_dir}
segment_file=${segment_file}

if [[ -z ${use_GDC_parameters} ]] ; then
	use_GDC_parameters="false"
else
	use_GDC_parameters=${use_GDC_parameters}
fi

if [[ -z ${conf_threshold} ]] ; then
	conf_threshold=0.75
else
	conf_threshold=${conf_threshold}
fi

if [[ -z ${broad_length} ]] ; then
	broad_length=0.98
else
	broad_length=${broad_length}
fi

if [[ ${use_GDC_parameters} == "false" ]] ; then
	amp_threshold=${amp_threshold}
	del_threshold=${del_threshold}
	sample_centering=${sample_centering}
	broad_analysis=${broad_analysis}
fi

. ${RUN_DIR}/config_files/common_config.sh

export XAPPLRESDIR=${GISTIC2_install_dir}/MATLAB_Compiler_Runtime/v83/X11/app-defaults:$XAPPLRESDIR
export LD_LIBRARY_PATH=${GISTIC2_install_dir}/MATLAB_Compiler_Runtime/v83/runtime/glnxa64:${GISTIC2_install_dir}/MATLAB_Compiler_Runtime/v83/bin/glnxa64:${GISTIC2_install_dir}/MATLAB_Compiler_Runtime/v83/sys/os/glnxa64:${LD_LIBRARY_PATH}

ref_genome=${GISTIC2_install_dir}/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat

GISTIC2=${GISTIC2_install_dir}/gistic2

mkdir -p ${output_dir}

if [[ ${use_GDC_parameters} == "false" ]] ; then
	${GISTIC2} \
	-b ${output_dir} \
	-seg ${segment_file} \
	-refgene ${ref_genome} \
	-scent ${sample_centering} \
	-broad ${broad_analysis} \
	-brlen ${broad_length} \
	-conf ${conf_threshold} \
	-ta ${amp_threshold} \
	-td ${del_threshold} \
	-savegene 1 \
	-rx 1 \
	-genegistic 0 \
	-v 30 
#	-twoside 1 \
#	-rx 0 
else
	${GISTIC2} \
	-b ${output_dir} \
	-seg ${segment_file} \
	-refgene ${ref_genome} \
	-ta 0.1 \
	-armpeel 1 \
	-brlen 0.7 \
	-cap 1.5 \
	-conf 0.99 \
	-td 0.1 \
	-genegistic 1 \
	-gcm extreme \
	-js 4 \
	-maxseg 2000 \
	-qvt 0.25 \
	-rx 0 \
	-savegene 1 \
	-broad 1
fi