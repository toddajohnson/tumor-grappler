#!/bin/bash

run_dir=${run_dir}
gpl_prefix=${gpl_prefix}
JVM_MEM=${jvmheap}

study_dir=$(basename ${run_dir})

. ${run_dir}/config_files/common_config.sh

write_status() { # Before logging initialised
	echo "$(date): $1" 1>&2
}

## liftOver, STAR 2.7.9a, samtools and bcftools 1.12, cutadapt are in /home/tjohnson/.local
export PATH="/home/tjohnson/.local/bin:${PATH}"
eval $(perl -I${HOME}/perl5/lib/perl5 -Mlocal::lib)

VCF2MAF="/home/tjohnson/tools/vcf2maf-1.6.21/vcf2maf.pl"

#export _JAVA_OPTIONS="-Xmx${jvmheap}"

GATK="/home/tjohnson/tools/gatk-4.1.9.0/gatk"

vep_path="/home/tjohnson/tools/ensembl-vep"
ref_path="/home/tjohnson/reference/HMF/38"
ref_fa="${ref_path}/refgenomes/Homo_sapiens.GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

tumor_normal_pair_info_file="${run_dir}/result_summaries/${gpl_prefix}/vep_annotation_sample_info.tsv"

## First line is header
let "LINENUM = ${SGE_TASK_ID} + 1"

tumor_normal_pair_info=$(sed -n "${LINENUM}p" ${tumor_normal_pair_info_file})

subject_id=$(echo "${tumor_normal_pair_info}" | awk -F "\t" '{print $1}')
tumor_id=$(echo "${tumor_normal_pair_info}" | awk -F "\t" '{print $2}')
normal_id=$(echo "${tumor_normal_pair_info}" | awk -F "\t" '{print $3}')

mkdir -p "${run_dir}/result/${gpl_prefix}/${subject_id}/purple/${tumor_id}/tmp"

if [[ "${run_type}" == "somatic" ]] ; then
	input_gzvcf="${run_dir}/result/${gpl_prefix}/${subject_id}/purple/${tumor_id}/${tumor_id}.purple.somatic.vcf.gz"
	output_maf="${run_dir}/result/${gpl_prefix}/${subject_id}/purple/${tumor_id}/${tumor_id}.purple.somatic.vep.maf"
elif [[ "${run_type}" == "germline" ]] ; then
	. ${cg_pipeline_dir}/common_config_files/germline_config_20210828.sh
	. ${cg_pipeline_dir}/common_config_files/ref_data_config_20210826.sh
	. ${cg_pipeline_dir}/common_scripts/vcf_processing.sh

	original_gzvcf="${run_dir}/result/${gpl_prefix}/${subject_id}/purple/${tumor_id}/${tumor_id}.purple.germline.vcf.gz"
	temp_original_vcf="${run_dir}/result/${gpl_prefix}/${subject_id}/purple/${tumor_id}/tmp/${tumor_id}.purple.germline.vcf"
	temp_newheader_gzvcf="${run_dir}/result/${gpl_prefix}/${subject_id}/purple/${tumor_id}/tmp/${tumor_id}.purple.germline.vcf.gz"
		
	filter_gzvcf="${run_dir}/result/${gpl_prefix}/${subject_id}/purple/${tumor_id}/${tumor_id}.purple.germline.with_AF_MAF_filter.vcf.gz"
	input_gzvcf="${run_dir}/result/${gpl_prefix}/${subject_id}/purple/${tumor_id}/${tumor_id}.purple.germline.AF_MAF_filtered.vcf.gz"
	output_maf="${run_dir}/result/${gpl_prefix}/${subject_id}/purple/${tumor_id}/${tumor_id}.purple.germline.vep.maf"
	
	${GATK} --java-options "-Xmx${JVM_MEM}" FixVcfHeader -I ${original_gzvcf} -O ${temp_newheader_gzvcf}

	filter_on_AF_MAF \
		${temp_newheader_gzvcf} \
		${filter_gzvcf} \
		${input_gzvcf} \
		"AF_MAF_filtered"
elif [[ "${run_type}" == "SVs" ]] ; then
	input_gzvcf="${run_dir}/result/${gpl_prefix}/merged.sv.vcf.gz"
	output_maf="${run_dir}/result/${gpl_prefix}/${subject_id}/purple/${tumor_id}/${tumor_id}.linx.sv.vep.maf"
fi

if [[ -f ${output_maf} ]] ; then
	rm ${output_maf}
fi

input_vcf=$(mktemp /tmp/${tumor_id}.XXXXXX.vcf)
gunzip -c ${input_gzvcf} > ${input_vcf}

echo "Running vcf2maf on ${run_type} vcf for ${tumor_id}"

if [[ "${run_type}" == "SVs" ]] ; then
	perl ${VCF2MAF} --vep-path ${vep_path} --ref-fasta ${ref_fa} --species homo_sapiens --ncbi-build GRCh38 --tumor-id ${tumor_id} --normal-id ${normal_id} --input-vcf ${input_vcf} --output-maf ${output_maf}
else
	perl ${VCF2MAF} --vep-path ${vep_path} --ref-fasta ${ref_fa} --species homo_sapiens --ncbi-build GRCh38 --tumor-id ${tumor_id} --normal-id ${normal_id} --retain-info REPORTED,PURPLE_AF,PURPLE_VCN, --input-vcf ${input_vcf} --output-maf ${output_maf}
fi

rm ${input_vcf}

echo "Finished vcf2maf"