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
#threads=12

# blockmb for 0.75*(72/2)*1024
# current run of small files use ~ 20GB, BQA uses 8-9GB ram
# to be safer, set to 32*1024
#bamsort_blockmb=41472
bamsort_blockmb=32768
bam_threads=12
bamsort_threads=6
bamsort_output_threads=6

base_logdir="${RUN_DIR}/log/01_preprocessing"
logdir="${base_logdir}/02_trim_and_map"

mkdir -p ${base_logdir} ${logdir}

logname="trim_and_map_fastq"

#qsub -tc 2 -t 1-35 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},use_temp_fastq="true",load_biobambam2_from_tools="true",BAMSORT_BLOCKMB=${bamsort_blockmb},BAM_THREADS=${bam_threads},BAMSORT_THREADS=${bamsort_threads},BAMSORT_OUTPUT_THREADS=${bamsort_output_threads} ${cg_pipeline_dir}/01_preprocessing/02_trim_and_map_nocutadapt.sh
#qsub -t 1 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},use_temp_fastq="true",load_biobambam2_from_tools="true",BAMSORT_BLOCKMB=${bamsort_blockmb},BAM_THREADS=${bam_threads},BAMSORT_THREADS=${bamsort_threads},BAMSORT_OUTPUT_THREADS=${bamsort_output_threads} ${cg_pipeline_dir}/01_preprocessing/02_trim_and_map_nocutadapt.sh
#qsub -t 2-35 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},use_temp_fastq="true",load_biobambam2_from_tools="true",BAMSORT_BLOCKMB=${bamsort_blockmb},BAM_THREADS=${bam_threads},BAMSORT_THREADS=${bamsort_threads},BAMSORT_OUTPUT_THREADS=${bamsort_output_threads} ${cg_pipeline_dir}/01_preprocessing/02_trim_and_map_nocutadapt.sh
#qsub -t 36-140 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},use_temp_fastq="true",load_biobambam2_from_tools="true",BAMSORT_BLOCKMB=${bamsort_blockmb},BAM_THREADS=${bam_threads},BAMSORT_THREADS=${bamsort_threads},BAMSORT_OUTPUT_THREADS=${bamsort_output_threads} ${cg_pipeline_dir}/01_preprocessing/02_trim_and_map_nocutadapt.sh
#qsub -t 141-175 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},use_temp_fastq="true",load_biobambam2_from_tools="true",BAMSORT_BLOCKMB=${bamsort_blockmb},BAM_THREADS=${bam_threads},BAMSORT_THREADS=${bamsort_threads},BAMSORT_OUTPUT_THREADS=${bamsort_output_threads} ${cg_pipeline_dir}/01_preprocessing/02_trim_and_map_nocutadapt.sh
#qsub -t 176-253 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},use_temp_fastq="true",load_biobambam2_from_tools="true",BAMSORT_BLOCKMB=${bamsort_blockmb},BAM_THREADS=${bam_threads},BAMSORT_THREADS=${bamsort_threads},BAMSORT_OUTPUT_THREADS=${bamsort_output_threads} ${cg_pipeline_dir}/01_preprocessing/02_trim_and_map_nocutadapt.sh
#qsub -t 254-340 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},use_temp_fastq="true",load_biobambam2_from_tools="true",BAMSORT_BLOCKMB=${bamsort_blockmb},BAM_THREADS=${bam_threads},BAMSORT_THREADS=${bamsort_threads},BAMSORT_OUTPUT_THREADS=${bamsort_output_threads} ${cg_pipeline_dir}/01_preprocessing/02_trim_and_map_nocutadapt.sh
#qsub -t 341-451 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},use_temp_fastq="true",load_biobambam2_from_tools="true",BAMSORT_BLOCKMB=${bamsort_blockmb},BAM_THREADS=${bam_threads},BAMSORT_THREADS=${bamsort_threads},BAMSORT_OUTPUT_THREADS=${bamsort_output_threads} ${cg_pipeline_dir}/01_preprocessing/02_trim_and_map_nocutadapt.sh
#qsub -t 452-502 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},use_temp_fastq="true",load_biobambam2_from_tools="true",BAMSORT_BLOCKMB=${bamsort_blockmb},BAM_THREADS=${bam_threads},BAMSORT_THREADS=${bamsort_threads},BAMSORT_OUTPUT_THREADS=${bamsort_output_threads} ${cg_pipeline_dir}/01_preprocessing/02_trim_and_map_nocutadapt.sh
#qsub -t 503-561 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},use_temp_fastq="true",load_biobambam2_from_tools="true",BAMSORT_BLOCKMB=${bamsort_blockmb},BAM_THREADS=${bam_threads},BAMSORT_THREADS=${bamsort_threads},BAMSORT_OUTPUT_THREADS=${bamsort_output_threads} ${cg_pipeline_dir}/01_preprocessing/02_trim_and_map_nocutadapt.sh
#qsub -t 562-671 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},use_temp_fastq="true",load_biobambam2_from_tools="true",BAMSORT_BLOCKMB=${bamsort_blockmb},BAM_THREADS=${bam_threads},BAMSORT_THREADS=${bamsort_threads},BAMSORT_OUTPUT_THREADS=${bamsort_output_threads} ${cg_pipeline_dir}/01_preprocessing/02_trim_and_map_nocutadapt.sh
#qsub -t 562-570 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},use_temp_fastq="true",load_biobambam2_from_tools="true",BAMSORT_BLOCKMB=${bamsort_blockmb},BAM_THREADS=${bam_threads},BAMSORT_THREADS=${bamsort_threads},BAMSORT_OUTPUT_THREADS=${bamsort_output_threads} ${cg_pipeline_dir}/01_preprocessing/02_trim_and_map_nocutadapt.sh
#qsub -t 572-579 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},use_temp_fastq="true",load_biobambam2_from_tools="true",BAMSORT_BLOCKMB=${bamsort_blockmb},BAM_THREADS=${bam_threads},BAMSORT_THREADS=${bamsort_threads},BAMSORT_OUTPUT_THREADS=${bamsort_output_threads} ${cg_pipeline_dir}/01_preprocessing/02_trim_and_map_nocutadapt.sh
#qsub -t 581-618 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},use_temp_fastq="true",load_biobambam2_from_tools="true",BAMSORT_BLOCKMB=${bamsort_blockmb},BAM_THREADS=${bam_threads},BAMSORT_THREADS=${bamsort_threads},BAMSORT_OUTPUT_THREADS=${bamsort_output_threads} ${cg_pipeline_dir}/01_preprocessing/02_trim_and_map_nocutadapt.sh
#qsub -t 620-644 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},use_temp_fastq="true",load_biobambam2_from_tools="true",BAMSORT_BLOCKMB=${bamsort_blockmb},BAM_THREADS=${bam_threads},BAMSORT_THREADS=${bamsort_threads},BAMSORT_OUTPUT_THREADS=${bamsort_output_threads} ${cg_pipeline_dir}/01_preprocessing/02_trim_and_map_nocutadapt.sh

#qsub -t 571 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},use_temp_fastq="true",load_biobambam2_from_tools="true",BAMSORT_BLOCKMB=${bamsort_blockmb},BAM_THREADS=${bam_threads},BAMSORT_THREADS=${bamsort_threads},BAMSORT_OUTPUT_THREADS=${bamsort_output_threads} ${cg_pipeline_dir}/01_preprocessing/02_trim_and_map_nocutadapt.sh
#qsub -t 580 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},use_temp_fastq="true",load_biobambam2_from_tools="true",BAMSORT_BLOCKMB=${bamsort_blockmb},BAM_THREADS=${bam_threads},BAMSORT_THREADS=${bamsort_threads},BAMSORT_OUTPUT_THREADS=${bamsort_output_threads} ${cg_pipeline_dir}/01_preprocessing/02_trim_and_map_nocutadapt.sh
#qsub -t 619 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},use_temp_fastq="true",load_biobambam2_from_tools="true",BAMSORT_BLOCKMB=${bamsort_blockmb},BAM_THREADS=${bam_threads},BAMSORT_THREADS=${bamsort_threads},BAMSORT_OUTPUT_THREADS=${bamsort_output_threads} ${cg_pipeline_dir}/01_preprocessing/02_trim_and_map_nocutadapt.sh
#qsub -t 645-656 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},use_temp_fastq="true",load_biobambam2_from_tools="true",BAMSORT_BLOCKMB=${bamsort_blockmb},BAM_THREADS=${bam_threads},BAMSORT_THREADS=${bamsort_threads},BAMSORT_OUTPUT_THREADS=${bamsort_output_threads} ${cg_pipeline_dir}/01_preprocessing/02_trim_and_map_nocutadapt.sh
#qsub -t 657-671 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},use_temp_fastq="true",load_biobambam2_from_tools="true",BAMSORT_BLOCKMB=${bamsort_blockmb},BAM_THREADS=${bam_threads},BAMSORT_THREADS=${bamsort_threads},BAMSORT_OUTPUT_THREADS=${bamsort_output_threads} ${cg_pipeline_dir}/01_preprocessing/02_trim_and_map_nocutadapt.sh
#qsub -tc 10 -t 672-820 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},use_temp_fastq="true",load_biobambam2_from_tools="true",BAMSORT_BLOCKMB=${bamsort_blockmb},BAM_THREADS=${bam_threads},BAMSORT_THREADS=${bamsort_threads},BAMSORT_OUTPUT_THREADS=${bamsort_output_threads} ${cg_pipeline_dir}/01_preprocessing/02_trim_and_map_nocutadapt.sh
#qsub -tc 10 -t 821-966 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},use_temp_fastq="true",load_biobambam2_from_tools="true",BAMSORT_BLOCKMB=${bamsort_blockmb},BAM_THREADS=${bam_threads},BAMSORT_THREADS=${bamsort_threads},BAMSORT_OUTPUT_THREADS=${bamsort_output_threads} ${cg_pipeline_dir}/01_preprocessing/02_trim_and_map_nocutadapt.sh
#qsub -tc 10 -t 967-1229 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},use_temp_fastq="true",load_biobambam2_from_tools="true",BAMSORT_BLOCKMB=${bamsort_blockmb},BAM_THREADS=${bam_threads},BAMSORT_THREADS=${bamsort_threads},BAMSORT_OUTPUT_THREADS=${bamsort_output_threads} ${cg_pipeline_dir}/01_preprocessing/02_trim_and_map_nocutadapt.sh
qsub -t 833 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},use_temp_fastq="true",load_biobambam2_from_tools="true",BAMSORT_BLOCKMB=${bamsort_blockmb},BAM_THREADS=${bam_threads},BAMSORT_THREADS=${bamsort_threads},BAMSORT_OUTPUT_THREADS=${bamsort_output_threads} ${cg_pipeline_dir}/01_preprocessing/02_trim_and_map_nocutadapt.sh
qsub -t 841 -pe def_slot ${threads} -N ${logname} -o ${logdir} -q '!mjobs_rerun.q' -j y -cwd -S /usr/bin/bash -v RUN_DIR=${RUN_DIR},use_temp_fastq="true",load_biobambam2_from_tools="true",BAMSORT_BLOCKMB=${bamsort_blockmb},BAM_THREADS=${bam_threads},BAMSORT_THREADS=${bamsort_threads},BAMSORT_OUTPUT_THREADS=${bamsort_output_threads} ${cg_pipeline_dir}/01_preprocessing/02_trim_and_map_nocutadapt.sh
