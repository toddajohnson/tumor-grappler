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

base_logdir="${RUN_DIR}/log/05_signature_analysis"
logdir="${base_logdir}/GISTIC2"

mkdir -p ${base_logdir}
mkdir -p ${logdir}

threads=12

##analysis_date="20211221"

## Default CN segment output changed to use log2(CN)-log2(2) instead of ploidy
#sample_centering="median"
#amp_threshold=0.1
#del_threshold=0.1
#broad_analysis=1

#annotation='AllSamples'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh
#
#annotation='ACD_RCC'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh
#
#annotation='clear_cell'
#segment_file="${RUN_DIR}/result/GISTIC2}/CN_segments/merged_CN_segs.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh
#
#annotation='papillary'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh
#
#annotation='clear_papillary'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh
#
#annotation='chromophobe'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh

#analysis_date="20220131"
#
# Default CN segment output changed to use log2(CN)-log2(2) instead of ploidy
#sample_centering="median"
#amp_threshold=0.1
#del_threshold=0.1
#broad_analysis=1
#
#annotation='ccRCC'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh
#
#annotation='non-ccRCC'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh
#
#
#analysis_date="20220131"
#
# Default CN segment output changed to use log2(CN)-log2(2) instead of ploidy
#sample_centering="median"
#amp_threshold=0.3
#del_threshold=0.3
#broad_analysis=1
#
#annotation='ccRCC'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh
#
#annotation='non-ccRCC'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh


#analysis_date="20220201"
#
# Default CN segment output changed to use log2(CN)-log2(2) instead of ploidy
#sample_centering="median"
#amp_threshold=0.1
#del_threshold=0.1
#broad_analysis=1
#
#annotation='AllSamples'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh

#
#annotation='ccRCC'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh
#
#annotation='non-ccRCC'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh


#analysis_date="20220201"
#
# Default CN segment output changed to use log2(CN)-log2(2) instead of ploidy
#sample_centering="median"
#amp_threshold=0.3
#del_threshold=0.3
#broad_analysis=1
#
#annotation='AllSamples'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh

#annotation='ccRCC'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh
#
#annotation='non-ccRCC'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh


#analysis_date="20220202"
#
# Default CN segment output changed to use log2(CN)-log2(2) instead of ploidy
#sample_centering="median"
#amp_threshold=0.1
#del_threshold=0.1
#broad_analysis=1
#
#analysis_date="${analysis_date}_${sample_centering}"
#
#annotation='AllSamples'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh
#
#annotation='ccRCC'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh
#
#annotation='non-ccRCC'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh
#
# Default CN segment output changed to use log2(CN)-log2(2) instead of ploidy
#amp_threshold=0.3
#del_threshold=0.3
#broad_analysis=1
#
#annotation='AllSamples'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh
#
#annotation='ccRCC'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh
#
#annotation='non-ccRCC'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh


#analysis_date="20220202"
# Default CN segment output changed to use log2(CN)-log2(2) instead of ploidy
#sample_centering="none"
#amp_threshold=0.1
#del_threshold=0.1
#broad_analysis=1
#analysis_date="${analysis_date}_${sample_centering}"
#
#annotation='AllSamples'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh
#
#annotation='ccRCC'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh
#
#annotation='non-ccRCC'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh


#analysis_date="20220202"
#
# Default CN segment output changed to use log2(CN)-log2(2) instead of ploidy
#sample_centering="median"
#amp_threshold=0.1
#del_threshold=0.1
#broad_analysis=1
#analysis_date="${analysis_date}_${sample_centering}"
#
#annotation='ACD-RCC'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh
#
#annotation='pRCC'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh
#
#annotation='other-RCC'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh


#analysis_date="20220202"
#
# Default CN segment output changed to use log2(CN)-log2(2) instead of ploidy
#sample_centering="median"
#amp_threshold=0.3
#del_threshold=0.3
#broad_analysis=1
#analysis_date="${analysis_date}_${sample_centering}"
#
#annotation='ACD-RCC'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh
#
#annotation='pRCC'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh

#annotation='other-RCC'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh


## Default CN segment output changed to use log2(CN)-log2(adj_diploid)
# and stronger threshold based on examining CN segments
#analysis_date="20220203"
#sample_centering="median"
#amp_threshold=0.3
#del_threshold=0.3
#broad_analysis=1
#analysis_date="${analysis_date}_${sample_centering}"
#
#annotation='AllSamples'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh
#
#analysis_date="20220203"
#sample_centering="none"
#amp_threshold=0.3
#del_threshold=0.3
#broad_analysis=1
#analysis_date="${analysis_date}_${sample_centering}"
#
#annotation='AllSamples'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh


#analysis_date="20220203"
#sample_centering="none"
#amp_threshold=0.3
#del_threshold=0.3
#broad_analysis=1
#analysis_date="${analysis_date}_${sample_centering}"
#
#annotation='non-ccRCC'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh

#analysis_date="20220204"
#sample_centering="none"
#amp_threshold=0.3
#del_threshold=0.3
#broad_analysis=1
#analysis_date="${analysis_date}_${sample_centering}_noChrX"
#
#annotation='AllSamples'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh


#analysis_date="20220204"
#sample_centering="none"
#amp_threshold=0.3
#del_threshold=0.3
#broad_analysis=1
#analysis_date="${analysis_date}_${sample_centering}_noChrX"
#
#annotation='non-ccRCC'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh

#analysis_date="20220205"
#sample_centering="none"
#amp_threshold=0.3
#del_threshold=0.3
#broad_analysis=1
#broad_length=0.5
#analysis_date="${analysis_date}_${sample_centering}_noChrX"
#
#annotation='AllSamples'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold},broad_length=${broad_length} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh
#
#analysis_date="20220205"
#sample_centering="none"
#amp_threshold=0.3
#del_threshold=0.3
#broad_analysis=1
#broad_length=0.5
#analysis_date="${analysis_date}_${sample_centering}_noChrX"
#
#annotation='non-ccRCC'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold},broad_length=${broad_length} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh

#analysis_date="20220205"
#sample_centering="none"
#amp_threshold=0.3
#del_threshold=0.3
#broad_analysis=1
#broad_length=0.98
#analysis_date="${analysis_date}_${sample_centering}_noChrX_brlen_${broad_length}"
#
#annotation='AllSamples'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold},broad_length=${broad_length} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh
#
#analysis_date="20220205"
#sample_centering="none"
#amp_threshold=0.3
#del_threshold=0.3
#broad_analysis=1
#broad_length=0.98
#analysis_date="${analysis_date}_${sample_centering}_noChrX_brlen_${broad_length}"
#
#annotation='non-ccRCC'
#segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
#output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
#qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold},broad_length=${broad_length} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh

analysis_date="20220205"
sample_centering="none"
amp_threshold=0.3
del_threshold=0.3
broad_analysis=1
broad_length=0.50
analysis_date="${analysis_date}_${sample_centering}_noChrX_brlen_${broad_length}"

annotation='AllSamples'
segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold},broad_length=${broad_length} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh

analysis_date="20220205"
sample_centering="none"
amp_threshold=0.3
del_threshold=0.3
broad_analysis=1
broad_length=0.50
analysis_date="${analysis_date}_${sample_centering}_noChrX_brlen_${broad_length}"

annotation='non-ccRCC'
segment_file="${RUN_DIR}/result/GISTIC2/CN_segments/merged_CN_segs.ver_20220131.${annotation}.txt"
output_dir="${RUN_DIR}/result/GISTIC2/${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
logname="log_01_run_GISTIC2_${analysis_date}_${amp_threshold}_${del_threshold}_${annotation}"
qsub -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -wd ${GISTIC2_install_dir} -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},output_dir=${output_dir},segment_file=${segment_file},sample_centering=${sample_centering},broad_analysis=${broad_analysis},amp_threshold=${amp_threshold},del_threshold=${del_threshold},broad_length=${broad_length} ${cg_pipeline_dir}/08_signature_analysis/GISTIC2/01_run_GISTIC2.sh
