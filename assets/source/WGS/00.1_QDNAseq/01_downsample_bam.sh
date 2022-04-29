#!/bin/bash

RUN_DIR=${RUN_DIR}

. ${RUN_DIR}/config_files/common_config.sh

# After checking in step 07_ markdup bams are consistence with fastq and earlier files
# Move to final bam folder
sample_index=${SGE_TASK_ID}

## samtools and bcftools 1.12, cutadapt are in /home/tjohnson/.local
export PATH="/home/tjohnson/.local/bin:${PATH}"

bam_dir="${RUN_DIR}/result/bam_final"
downsampled_bam_dir="${RUN_DIR}/result/downsampled_bam"
SAMPLE_INFO_FILE="${RUN_DIR}/config_files/sample_final_bam_read_cts.csv"

mkdir -p ${downsampled_bam_dir}

## LINENUM reflects line not subject index;  needs to changes as more samples added
## First line is header
let "LINENUM = ${sample_index} + 1"
SAMPLE_INFO=$(sed -n "${LINENUM}p" ${SAMPLE_INFO_FILE})
subject_id=$(echo ${SAMPLE_INFO} | awk -F "," '{print $1}')
sample_id=$(echo ${SAMPLE_INFO} | awk -F "," '{print $2}')
read_pair_ct=$(echo ${SAMPLE_INFO} | awk -F "," '{print $3}')

echo "Getting ready to downsample ${sample_id} bam file."

input_bam_file="${bam_dir}/${sample_id}.bam"
downsampled_bam_file="${downsampled_bam_dir}/${sample_id}.downsampled.bam"

target_read_pair_ct=10000000
sampling_seed=569
## 1540 is: read unmapped (0x4); read fails platform/vendor quality checks (0x200); read is PCR or optical duplicate (0x400)
filter_bits=1540
mapq=30

frac=$(bc -l <<< "scale=4; ${target_read_pair_ct}/${read_pair_ct}")

echo "The target read pair ct. was ${frac} of the file read pair ct."

if (( $(echo "${frac} < 1" | bc -l) )); then
	echo "Downsampling bam for ${sample_id} to about ${target_read_pair_ct} of ${read_pair_ct} read pairs"
	samtools view -h -q 30 -F ${filter_bits} ${input_bam_file} | samtools view -s ${sampling_seed}${frac} -bo ${downsampled_bam_file}
elif (( $(echo "${frac} >= 1" | bc -l) )); then
	echo "Bam file ${input_bam_file} for ${sample_id} has ${read_pair_ct} read pairs.  That is less the target of ${target_read_pair_ct}."
	echo "Filtering bam file for MAPQ of ${mapq} and filter on ${filter_bits}"
	samtools view -h -q ${mapq} -F ${filter_bits} ${input_bam_file} | samtools view -s ${sampling_seed}${frac} -bo ${downsampled_bam_file}
fi

if [[ -s ${downsampled_bam_file} ]]; then
	echo "Indexing bam file"
	samtools index ${downsampled_bam_file}
fi
