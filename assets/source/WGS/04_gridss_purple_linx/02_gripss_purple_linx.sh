#!/bin/bash

GPL_PREFIX=$(basename "$PWD")

run_dir="${RUN_DIR}"
PGL_THREAD_NUM="${PGL_THREAD_NUM}"
PGL_JVM_MEM="${PGL_JVM_MEM}"

#run_set=${run_set}

if [[ -z ${snvvcf_split_by_tumor} ]] ; then
	snvvcf_split_by_tumor="false"
else
	snvvcf_split_by_tumor=${snvvcf_split_by_tumor}
fi

if [[ -z ${ignore_bams} ]] ; then
	ignore_bams="false"
else
	ignore_bams=${ignore_bams}
fi

if [[ -z ${markdup_bams} ]] ; then
	markdup_bams="false"
else
	markdup_bams=${markdup_bams}
fi

if [[ -z ${germline_target} ]] ; then
	germline_target="panel"
else
	germline_target=${germline_target}
fi

if [[ -z ${somatic_filter_string} ]] ; then
	somatic_filter_string="filtered"
else
	somatic_filter_string=${somatic_filter_string}
fi

if [[ -z ${germline_filter_string} ]] ; then
	germline_filter_string="filtered"
#	germlne_filter_string="AF_MAF_filtered"
else
	germline_filter_string=${germline_filter_string}
fi

if [[ -z ${write_neo_epitopes} ]] ; then
	write_neo_epitopes="false"
else
	write_neo_epitopes=${write_neo_epitopes}
fi

if [[ -z ${min_ploidy} ]] ; then
	min_ploidy=1
else
	min_ploidy=${min_ploidy}
fi

if [[ -z ${max_ploidy} ]] ; then
	max_ploidy=8
else
	max_ploidy=${max_ploidy}
fi

if [[ -z ${min_purity} ]] ; then
	min_purity=0.08
else
	min_purity=${min_purity}
fi

if [[ -z ${max_purity} ]] ; then
	max_purity=1
else
	max_purity=${max_purity}
fi

. ${run_dir}/config_files/common_config.sh
#. ${cg_pipeline_dir}/common_config_files/ref_data_config_20210826.sh
. ${cg_pipeline_dir}/common_config_files/ref_data_config_20211006.sh

module use /usr/local/package/modulefiles/
#module load java/8
module load java/11
module load bwa/0.7.17

export PATH="/home/tjohnson/.local/bin:/home/tjohnson/.local/RepeatMasker:/home/tjohnson/tools/kraken2:/home/tjohnson/tools/circos-0.69-9/bin:${PATH}"

module load R/4.0.2

#export _JAVA_OPTIONS="-Xmx${PGL_JVM_MEM}"

#REFERENCE_GENOME_VERSION="HG38"
REFERENCE_GENOME_VERSION="38"
REF_DATA_PATH="/home/tjohnson/reference/HMF/38"

## Need to determine how to pass VIRUSBREAKENDDB to the called script
VIRUSBREAKENDDB="/home/tjohnson/reference/VIRUSBreakend/virusbreakenddb"

TUMOR_NORMAL_PAIR_INFO_FILE="${RUN_DIR}/config_files/tumor_normal_pair_info_GPL.csv"

## First line is header
let "LINENUM = ${SGE_TASK_ID} + 1"

TUMOR_NORMAL_PAIR_INFO=$(sed -n "${LINENUM}p" ${TUMOR_NORMAL_PAIR_INFO_FILE})

SUBJECT_ID=$(echo ${TUMOR_NORMAL_PAIR_INFO} | awk -F "," '{print $1}')
TUMOR_SAMPLE_ID=$(echo ${TUMOR_NORMAL_PAIR_INFO} | awk -F "," '{print $4}')
NORMAL_SAMPLE_ID=$(echo ${TUMOR_NORMAL_PAIR_INFO} | awk -F "," '{print $6}')

GRIDSS_OUTPUT_DIR=${RUN_DIR}/result/${GPL_PREFIX}
OUTPUT_PATH="${GRIDSS_OUTPUT_DIR}/${SUBJECT_ID}"

cd ${OUTPUT_PATH}

if [[ "${markdup_bams}" == "true" ]] ; then
	TUMOR_BAM_PATH="${RUN_DIR}/result/bam_markdup/${TUMOR_SAMPLE_ID}.markdup.bam"
else
	TUMOR_BAM_PATH="${RUN_DIR}/result/bam_final/${TUMOR_SAMPLE_ID}.bam"
fi

if [[ "${snvvcf_split_by_tumor}" == "true" ]]; then
	MERGED_ANNOTATED_SOMATIC_VCF=${RUN_DIR}/result/mutations/${SUBJECT_ID}/SAGE-${SAGE_VERSION}/somatic/${TUMOR_SAMPLE_ID}.somatic.sage.${somatic_filter_string}.vcf.gz
else
	MERGED_ANNOTATED_SOMATIC_VCF=${RUN_DIR}/result/mutations/${SUBJECT_ID}/SAGE-${SAGE_VERSION}/somatic/${SUBJECT_ID}.somatic.sage.${somatic_filter_string}.vcf.gz
fi

if [[ "${NORMAL_SAMPLE_ID}" == "none" ]] ; then
	NORMAL_BAM_PATH=""
	NORMAL_SAMPLE_ID=""
	purple_args=""
else
	if [[ "${markdup_bams}" == "true" ]] ; then
		NORMAL_BAM_PATH="${RUN_DIR}/result/bam_markdup/${NORMAL_SAMPLE_ID}.markdup.bam"
	else
		NORMAL_BAM_PATH="${RUN_DIR}/result/bam_final/${NORMAL_SAMPLE_ID}.bam"
	fi
	
	if [[ "${germline_target}" == "panel" ]] ; then
		MERGED_ANNOTATED_GERMLINE_VCF=${RUN_DIR}/result/mutations/${SUBJECT_ID}/SAGE-${SAGE_VERSION}/germline/${SUBJECT_ID}.germline.sage.${germline_filter_string}.vcf.gz
		purple_args="-germline_vcf ${MERGED_ANNOTATED_GERMLINE_VCF}"
	elif [[ "${germline_target}" == "genomewide" ]] ; then
		MERGED_ANNOTATED_GERMLINE_VCF=${RUN_DIR}/result/mutations/${SUBJECT_ID}/SAGE-${SAGE_VERSION}/germline_${germline_target}/${SUBJECT_ID}.germline.sage.${germline_filter_string}.vcf.gz
		purple_args="-germline_vcf ${MERGED_ANNOTATED_GERMLINE_VCF}"
	fi
fi

purple_args="-min_ploidy ${min_ploidy} -max_ploidy ${max_ploidy} -min_purity ${min_purity} -max_purity ${max_purity} ${purple_args}"

#if [[ ${run_set} == 1 ]]; then
#	. ${cg_pipeline_dir}/04_gridss_purple_linx/04_gripss_purple_linx_run_settings_1.sh
#elif [[ ${run_set} == 2 ]]; then
#	. ${cg_pipeline_dir}/04_gridss_purple_linx/04_gripss_purple_linx_run_settings_2.sh
#fi

## One line and some logic added to gripps-purple-linx.sh to
## quickly allow for just running purple and linx with no bam files present
if [[ "${ignore_bams}" == "true" ]] ; then
	export IGNORE_BAMS="true"
else
	export IGNORE_BAMS="false"
fi


gripps_purple_linx_cmd=${cg_pipeline_dir}/04_gridss_purple_linx/gripss-purple-linx.sh

echo "Running gripss-purple-linx using ${gripps_purple_linx_cmd}"

echo "Somatic VCF: ${MERGED_ANNOTATED_SOMATIC_VCF}"
echo "Germline VCF: ${MERGED_ANNOTATED_GERMLINE_VCF}"

if [[ "${write_neo_epitopes}" == "true" ]] ; then
	echo "Setting arguments to write neo-epitope file."
	linx_args="-write_neo_epitopes"

	if [[ -f "${RNAseq_RUN_DIR}/config_files/Isofox_sample_list.csv" ]] ; then
		echo "Isofox sample list file exits.  Checking if current sample is in the list."
		## If Linx fails, it may be because of misformatted isofox.combined_fusions.csv; check columns and header
		isofox_sample_check=$(grep ${TUMOR_SAMPLE_ID} "${RNAseq_RUN_DIR}/config_files/Isofox_sample_list.csv" | wc -l)
	
		if [[ "${isofox_sample_check}" == "1" ]] ; then
			echo "Sample found.  Adding RNA fusion files to Linx arguments."
			rna_fusions_file="${RNAseq_RUN_DIR}/result/isofox/isofox.combined_fusions.csv"
			rna_file_source="ISOFOX"
			linx_args="${linx_args} -rna_fusions_file ${rna_fusions_file} -rna_file_source ${rna_file_source}"
		fi
	fi
else
	linx_args=""
fi

if [[ "$GRIPSS_VERSION" == "1.11" ]] ; then
	gripss_args="gripss.somatic"
elif [[ "$GRIPSS_VERSION" == "v2.0_beta" ]] || [["$GRIPSS_VERSION" == "v2.0"]]; then
	gripss_args="gripss"
fi

${gripps_purple_linx_cmd} \
	-n "${NORMAL_BAM_PATH}" \
	-t "${TUMOR_BAM_PATH}" \
	--sample "${SUBJECT_ID}" \
	--normal_sample "${NORMAL_SAMPLE_ID}" \
	--tumour_sample "${TUMOR_SAMPLE_ID}" \
	--snvvcf "${MERGED_ANNOTATED_SOMATIC_VCF}" \
	--purple_args "${purple_args}" \
	--ref_genome_version "${REFERENCE_GENOME_VERSION}" \
	--ref_dir "${REF_DATA_PATH}" \
	--install_dir "${install_dir}" \
	--threads "${PGL_THREAD_NUM}" \
	--jvmheap "${PGL_JVM_MEM}" \
	--output_dir "${OUTPUT_PATH}" \
	--linx_args "${linx_args}" \
	--gripss_args "$gripss_args"
