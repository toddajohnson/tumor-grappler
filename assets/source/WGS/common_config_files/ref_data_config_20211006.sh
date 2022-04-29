#!/bin/bash

ref_data_path="/home/tjohnson/reference/HMF/38"
ref_genome_version="hg38"
#need to figure out where, SAGE, GRIDSS, PURPLE used hg38(SAGE) vs other versions
#ref_genome_version="38"
ref_genome="${ref_data_path}/refgenomes/Homo_sapiens.GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"


dbsnp_dir="/home/tjohnson/reference/dbsnp/human_9606_bb151_GRCh38p7"
dbsnp_common_vcf="${dbsnp_dir}/common_all_20180418.vcf.gz"

alfa_dir="/home/tjohnson/reference/alfa/20210106"
alfa_vcf="${alfa_dir}/freq.vcf.gz"
alfa_AF_MAF_txt=${alfa_dir}/merged_alfa_AF.txt.gz

alfa_AF_MAF_annotation_header_vcf="${alfa_dir}/alfa_annotation_header.vcf"

dbSNP_dir="/home/tjohnson/reference/dbSNP/vcf/20210521"
dbSNP_vcf="${dbSNP_dir}/GCF_000001405.39.gz"
dbSNP_decomposed_bcf="${dbSNP_dir}/GCF_000001405.39.decomposed.bcf"
#dbSNP_extracted_bcf="${dbSNP_dir}/GCF_000001405.39.FREQ_extracted.bcf"

ucsc_to_ncbi_chr_conv_file="${cg_pipeline_dir}/common_config_files/ucsc_to_ncbi_chr_name_conv.txt"
ncbi_to_ucsc_chr_conv_file="${cg_pipeline_dir}/common_config_files/ncbi_to_ucsc_chr_name_conv.txt"

ucsc_to_NC_accession_conv_file="${cg_pipeline_dir}/common_config_files/ucsc_chr_name_to_NC_accession_conv.txt"
NC_accession_to_ucsc_conv_file="${cg_pipeline_dir}/common_config_files/NC_accession_to_ucsc_chr_name_conv.txt"

alfa_pop_id_to_pop_name_conv="${cg_pipeline_dir}/common_config_files/alfa_pop_id_to_pop_name_conv.txt"

humandb_path="/home/tjohnson/tools/annovar/humandb/"
snpeff_tool_dir="/home/tjohnson/tools/snpEff/"

snpeff_anno_db="GRCh38.manual.104"

clinvar_vcf="/home/tjohnson/tools/snpEff/clinvar/GRCh38/clinvar.with.chr.vcf.gz"

sage_resource_dir="${ref_data_path}/dbs/sage"

mappability_bed="${ref_data_path}/dbs/mappability/out_150.mappability.38.bed.gz"
mappability_hdr="${ref_data_path}/dbs/mappability/mappability.hdr"
known_blacklist_bed="${sage_resource_dir}/38/KnownBlacklist.germline.38.bed.gz"
known_blacklist_vcf="${sage_resource_dir}/38/KnownBlacklist.germline.38.vcf.gz"
germlinePon="${sage_resource_dir}/38/SageGermlinePon.JP.149x.38.vcf.gz"

if [[ "${call_target}" == "germline" ]]; then
	hotspots_vcf="${sage_resource_dir}/38/KnownHotspots.germline.38.vcf.gz"
	coverage_bed="${sage_resource_dir}/38/CoverageCodingPanel.germline.38.bed"
	coding_panel_bed="${sage_resource_dir}/38/ActionableCodingPanel.germline.38.bed"
elif [[ "${call_target}" == "somatic" ]]; then
	hotspots_vcf="${sage_resource_dir}/38/KnownHotspots.somatic.38.vcf.gz"
	coding_panel_bed="${sage_resource_dir}/38/ActionableCodingPanel.somatic.38.bed"
else
	hotspots_vcf="${sage_resource_dir}/38/KnownHotspots.germline.38.vcf.gz"
	coverage_bed="${sage_resource_dir}/38/CoverageCodingPanel.germline.38.bed"
	coding_panel_bed="${sage_resource_dir}/38/ActionableCodingPanel.germline.38.bed"
fi

high_confidence_bed="${sage_resource_dir}/38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed"
