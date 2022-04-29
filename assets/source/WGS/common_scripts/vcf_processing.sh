#germlinePon=''
#mappability_bed=''
#mappability_hdr=''
#clinvar_vcf=''
#alfa_vcf=''
#ucsc_to_NC_accession_conv_file=''
#NC_accession_to_ucsc_conv_file=''
#alfa_pop_id_to_pop_name_conv_file=''
#snpeff_tool_dir=''
#humandb_path=''
#threads=''
#JVM_MEM=''
#tmp_output_dir=''

remove_previous_annotation_output(){
	if [[ -f ${output_dir}/${subject_id}.${call_target}.sage.PON_filtered.vcf.gz ]]; then
		write_status "Removing previous PON filtered file"
		rm ${output_dir}/${subject_id}.${call_target}.sage.PON_filtered.vcf.gz ${output_dir}/${subject_id}.${call_target}.sage.PON_filtered.vcf.gz.tbi
	fi

	if [[ -f ${output_dir}/${subject_id}.${call_target}.sage.PON_filtered.hg38_multianno.vcf.gz ]]; then
		write_status "Removing previous PON annovar annotations file"
		rm ${output_dir}/${subject_id}.${call_target}.sage.PON_filtered.hg38_multianno.vcf.gz ${output_dir}/${subject_id}.${call_target}.sage.PON_filtered.hg38_multianno.vcf.gz.tbi
		rm ${output_dir}/${subject_id}.${call_target}.sage.PON_filtered.avinput ${output_dir}/${subject_id}.${call_target}.sage.PON_filtered.hg38_multianno.txt
	fi

	if [[ -f ${output_dir}/${subject_id}.${call_target}.sage.annotated.vcf.gz ]]; then
		write_status "Removing previous annotated file"
		rm ${output_dir}/${subject_id}.${call_target}.sage.annotated.vcf.gz ${output_dir}/${subject_id}.${call_target}.sage.annotated.vcf.gz.tbi
	fi
	
	if [[ -f ${output_dir}/${subject_id}.${call_target}.sage.annotated.chrM.vcf.gz ]]; then
		write_status "Removing previous annotated chrM file"
		rm ${output_dir}/${subject_id}.${call_target}.sage.annotated.chrM.vcf.gz ${output_dir}/${subject_id}.${call_target}.sage.annotated.chrM.vcf.gz.tbi
	fi
}

remove_previous_filtered_output(){
	if [[ -f ${output_dir}/${subject_id}.${call_target}.sage.filtered.premod.vcf.gz ]]; then
		write_status "Removing previous pre-mod filtered file"
		rm ${output_dir}/${subject_id}.${call_target}.sage.filtered.premod.vcf.gz ${output_dir}/${subject_id}.${call_target}.sage.filtered.premod.vcf.gz.tbi
	fi
	
	if [[ -f ${output_dir}/${subject_id}.${call_target}.sage.weak_filtered.vcf.gz ]]; then
		write_status "Removing previous weak_filtered file"
		rm ${output_dir}/${subject_id}.${call_target}.sage.weak_filtered.vcf.gz ${output_dir}/${subject_id}.${call_target}.sage.weak_filtered.vcf.gz.tbi
		rm ${output_dir}/snpEff_summary_${snpeff_anno_db}.weak_filtered.html
		rm ${output_dir}/snpEff_summary_${snpeff_anno_db}.weak_filtered.genes.txt
	fi
	
	if [[ -f ${output_dir}/${subject_id}.${call_target}.sage.filtered.vcf.gz ]]; then
		write_status "Removing previous filtered file"
		rm ${output_dir}/${subject_id}.${call_target}.sage.filtered.vcf.gz ${output_dir}/${subject_id}.${call_target}.sage.filtered.vcf.gz.tbi
		rm ${output_dir}/snpEff_summary_${snpeff_anno_db}.filtered.html
		rm ${output_dir}/snpEff_summary_${snpeff_anno_db}.filtered.genes.txt
		
		if [[ -f ${output_dir}/snpEff_summary_${snpeff_anno_db}.html ]]; then
			rm ${output_dir}/snpEff_summary_${snpeff_anno_db}.html
			rm ${output_dir}/snpEff_summary_${snpeff_anno_db}.genes.txt
		fi
	fi
	
	if [[ -f ${output_dir}/${subject_id}.${call_target}.sage.with_AF_MAF_filter.vcf.gz ]]; then
		write_status "Removing previous with_AF_MAF_filter file"
		rm ${output_dir}/${subject_id}.${call_target}.sage.with_AF_MAF_filter.vcf.gz ${output_dir}/${subject_id}.${call_target}.sage.with_AF_MAF_filter.vcf.gz.tbi
	fi

	if [[ -f ${output_dir}/${subject_id}.${call_target}.sage.AF_MAF_filtered.vcf.gz ]]; then
		write_status "Removing previous AF_MAF_filtered file"
		rm ${output_dir}/${subject_id}.${call_target}.sage.AF_MAF_filtered.vcf.gz ${output_dir}/${subject_id}.${call_target}.sage.AF_MAF_filtered.vcf.gz.tbi
		rm ${output_dir}/snpEff_summary_${snpeff_anno_db}.AF_MAF_filtered.html
		rm ${output_dir}/snpEff_summary_${snpeff_anno_db}.AF_MAF_filtered.genes.txt
	fi
	
	if [[ -f ${output_dir}/${subject_id}.${call_target}.sage.${tumor_names}_filtered.vcf.gz ]]; then
		write_status "Removing previous ${tumor_names}_filtered file"
		rm ${output_dir}/${subject_id}.${call_target}.sage.${tumor_names}_filtered.vcf.gz ${output_dir}/${subject_id}.${call_target}.sage.${tumor_names}_filtered.vcf.gz.tbi
		rm ${output_dir}/snpEff_summary_${snpeff_anno_db}.${tumor_names}_filtered.html
		rm ${output_dir}/snpEff_summary_${snpeff_anno_db}.${tumor_names}_filtered.genes.txt
	fi
	
	if [[ -f ${output_dir}/${subject_id}.${call_target}.sage.HMF_filtered.vcf.gz ]]; then
		write_status "Removing previous PON filtered file"
		rm ${output_dir}/${subject_id}.${call_target}.sage.HMF_filtered.vcf.gz ${output_dir}/${subject_id}.${call_target}.sage.HMF_filtered.vcf.gz.tbi
		rm ${output_dir}/snpEff_summary_${snpeff_anno_db}.HMF_filtered.html
		rm ${output_dir}/snpEff_summary_${snpeff_anno_db}.HMF_filtered.genes.txt
	fi
	
	if [[ -f ${output_dir}/${subject_id}.${call_target}.sage.relaxed_filtered.vcf.gz ]]; then
		write_status "Removing previous relaxed filtered filtered file"
		rm ${output_dir}/${subject_id}.${call_target}.sage.relaxed_filtered.vcf.gz ${output_dir}/${subject_id}.${call_target}.sage.relaxed_filtered.vcf.gz.tbi
	fi
}

annotate_mappability(){
	write_status "Annotating with SAGE germline mappabability"
	
	local input_bcf="$1"
	local output_bcf="$2"
	
	write_status "Running mappability annotation"
	bcftools annotate -a ${mappability_bed} -h ${mappability_hdr} -c CHROM,FROM,TO,-,MAPPABILITY --threads ${threads} -Ob -o ${output_bcf} ${input_bcf}
	
	write_status "Indexing mappability annotated BCF ${output_bcf} and removing input files"
	bcftools index ${output_bcf}
	rm ${input_bcf} ${input_bcf}.csi
}

annotate_PON_mappability(){
	write_status "Annotating with SAGE germline PON and mappabability"
	
	local input_bcf="$1"
	local output_bcf="$2"
	
	write_status "Running PON and mappability annotation"
	bcftools annotate -a ${germlinePon} -c PON_COUNT,PON_MAX -Ou ${input_bcf} \
	| bcftools annotate -a ${mappability_bed} -h ${mappability_hdr} -c CHROM,FROM,TO,-,MAPPABILITY -Ou - \
	| bcftools filter -e 'PON_COUNT!="." && INFO/TIER="HOTSPOT" && PON_MAX>=5 && PON_COUNT >= 5' -s PON -m+ -Ou - \
	| bcftools filter -e 'PON_COUNT!="." && INFO/TIER="PANEL" && PON_MAX>=5 && PON_COUNT >= 2' -s PON -m+ -Ou - \
	| bcftools filter -e 'PON_COUNT!="." && INFO/TIER!="HOTSPOT" && INFO/TIER!="PANEL" && PON_COUNT >= 2' -s PON -m+ --threads ${threads} -Ob -o ${output_bcf} -
	
	write_status "Indexing PON/mappability annotated BCF ${output_bcf} and removing input files"
	bcftools index ${output_bcf}
	rm ${input_bcf} ${input_bcf}.csi
}

annotate_clinvar(){
	write_status "Annotating with ClinVar CLNSIG and CLNSIGCONF"
	
	local input_bcf="$1"
	local output_bcf="$2"
	
	write_status "Running clinvar annotation of ${input_bcf}"
	bcftools annotate -a ${clinvar_vcf} -c INFO/CLNSIG,INFO/CLNSIGCONF --threads ${threads} -Ob -o ${output_bcf} ${input_bcf} 
	
	write_status "Indexing clinvar annotated BCF ${output_bcf} and removing input files"
	bcftools index ${output_bcf}
	
	rm ${input_bcf} ${input_bcf}.csi
}

annotate_blacklisted(){
	write_status "Annotating with blacklisted regionsF"
	
	local input_bcf="$1"
	local output_vcf="$2"
	
	local blacklisted_region_bcf="${tmp_output_dir}/blacklist_region_annotated.bcf"

	write_status "Annotating blacklisted regions"
	bcftools annotate -a ${known_blacklist_bed} -m BLACKLIST_BED -c CHROM,FROM,TO --threads ${threads} -Ob -o ${blacklisted_region_bcf} ${input_bcf}
	bcftools index ${blacklisted_region_bcf}
	rm ${input_bcf} ${input_bcf}.csi

	write_status "Annotating blacklisted variants and running snpEff with ${snpeff_anno_db} database"
	bcftools annotate -a ${known_blacklist_vcf} -m BLACKLIST_VCF --threads ${threads} -Ob -o ${output_vcf} ${blacklisted_region_bcf}
	tabix -p vcf ${output_vcf}
	rm ${blacklisted_region_bcf} ${blacklisted_region_bcf}.csi
}

remove_PON_filtered(){
	write_status "Filtering out PON somatic driver candidate variant tiers"	
	
	local input_bcf="$1"
	local output_vcf="$2"
	
	write_status "Running PON filter exclusion of ${input_bcf}"
	bcftools view -i 'FILTER!~"PON"' --threads ${threads} -Oz -o ${output_vcf} ${input_bcf}
	
	write_status "Indexing PON filtered VCF ${output_vcf} and removing input files"
	tabix -p vcf ${output_vcf}
	rm ${input_bcf} ${input_bcf}.csi
}


# Most MNVs do not get labelled in the next two annotationes
## Probably, would need to split MNVs into SNPs in input germline file,
## but it is only about 1.2% of variants.

annotate_dbSNP(){
	write_status "Annotating with ALFA variant frequencies"

	local input_vcf="$1"
	local output_vcf="$2"

	write_status "Annotating for dbSNP variants that are in subject vcf"
	bcftools annotate -a "${dbSNP_decomposed_bcf}" -c INFO/RS,INFO/VC,INFO/dbGaP_PopFreq,INFO/TOPMED,INFO/GnomAD,INFO/GnomAD_exomes,INFO/1000Genomes,INFO/TOMMO,INFO/KOREAN,INFO/Korea1K --threads ${threads} -Oz -o ${output_vcf} ${input_vcf}
	tabix -p vcf ${output_vcf}
}

annotate_ALFA_AF(){
	write_status "Annotating with ALFA variant frequencies"
	
	local input_vcf="$1"
	local output_vcf="$2"

	write_status "Annotating for ALFA vcf that are in subject vcf"
	bcftools annotate -a "${alfa_AF_MAF_txt}" -c CHROM,POS,ID,REF,ALT,AF_EAS,AF_TOT,MAF_EAS,MAF_TOT -h "${alfa_AF_MAF_annotation_header_vcf}" --threads ${threads} -Oz -o ${output_vcf} ${input_vcf}
	tabix -p vcf ${output_vcf}
}

annotate_ALFA_dbSNP(){
	write_status "Annotating with ALFA variant frequencies"
	
	local input_vcf="$1"
	local output_vcf="$2"

	local ALFA_annotated_bcf="${tmp_output_dir}/ALFA_annotated.bcf"
	
	write_status "Annotating for ALFA vcf that are in subject vcf"
	bcftools annotate -a "${alfa_AF_MAF_txt}" -c CHROM,POS,ID,REF,ALT,AF_EAS,AF_TOT,MAF_EAS,MAF_TOT -h "${alfa_AF_MAF_annotation_header_vcf}" --threads ${threads} -Ob -o ${ALFA_annotated_bcf} ${input_vcf}
	bcftools index ${ALFA_annotated_bcf}
	
	write_status "Annotating for dbSNP variants that are in subject vcf"
	bcftools annotate -a "${dbSNP_decomposed_bcf}" -c INFO/RS,INFO/VC,INFO/dbGaP_PopFreq,INFO/TOPMED,INFO/GnomAD,INFO/GnomAD_exomes,INFO/1000Genomes,INFO/TOMMO,INFO/KOREAN,INFO/Korea1K --threads ${threads} -Oz -o ${output_vcf} ${ALFA_annotated_bcf}
	tabix -p vcf ${output_vcf}
	
	rm ${ALFA_annotated_bcf} ${ALFA_annotated_bcf}.csi
}

annotate_with_annovar(){
	write_status "Annotating with annovar for allele frequencies"
	
	local input_vcf="$1"
	local output_prefix="$2"
	
	table_annovar.pl -thread ${threads} ${input_vcf} ${humandb_path} -buildver hg38 -out ${output_prefix} -remove -protocol gnomad30_genome -operation f -nastring . -vcfinput -polish
	bgzip ${output_prefix}.hg38_multianno.vcf
	tabix -p vcf ${output_prefix}.hg38_multianno.vcf.gz
}

filter_PASS_variants(){
	write_status "Filtering for variants with PASS in filter"
	
	local input_vcf="$1"
	local output_vcf="$2"
	local output_filter_string="$3"
	
	bcftools view -f PASS -Ov ${input_vcf} \
	| java -Xmx${JVM_MEM} -jar ${snpeff_tool_dir}/snpEff.jar ann -s snpEff_summary_${snpeff_anno_db}.${output_filter_string}.html ${snpeff_anno_db} \
	| bgzip -c > ${output_vcf}
	tabix -p vcf ${output_vcf}
}

filter_on_AF_MAF(){
	write_status "Filtering PASS VCF for AF and MAF; remove high MAF and low-MAF/high-AF (fixed variants)"

	local input_vcf="$1"
	local output_with_filter_vcf="$2"
	local output_vcf="$3"
	local output_filter_string="$4"

	bcftools filter -e "INFO/AF_EAS>=${filter_max_maf} & INFO/AF_EAS<${filter_max_af}" -s EAS_MAF -m+ -Ou ${input_vcf} \
	| bcftools filter -e "INFO/AF_EAS>=${filter_max_af}" -s fixed_EAS -m+ -Ou - \
	| bcftools filter -e "INFO/AF_TOT>=${filter_max_maf} & INFO/AF_TOT<${filter_max_af}" -s TOT_MAF -m+ -Ou - \
	| bcftools filter -e "INFO/AF_TOT>=${filter_max_af}" -s fixed_TOT -m+ -Ou - \
	| bcftools filter -e "INFO/TOMMO[1]>=${filter_max_maf} & INFO/TOMMO[1]<${filter_max_af}" -s TOMMO_MAF -m+ -Ou - \
	| bcftools filter -e "INFO/TOMMO[1]>=${filter_max_af}" -s fixed_TOMMO -m+ -Ou - \
	| bcftools filter -e "INFO/KOREAN[1]>=${filter_max_maf} & INFO/KOREAN[1]<${filter_max_af}" -s KOREAN_MAF -m+ -Ou - \
	| bcftools filter -e "INFO/KOREAN[1]>=${filter_max_af}" -s fixed_KOREAN -m+ -Ou - \
	| bcftools filter -e "INFO/Korea1K[1]>=${filter_max_maf} & INFO/Korea1K[1]<${filter_max_af}" -s Korea1K_MAF -m+ -Ou - \
	| bcftools filter -e "INFO/Korea1K[1]>=${filter_max_af}" -s fixed_Korea1K -m+ -Ou - \
	| bcftools filter -e "INFO/TOPMED[1]>=${filter_max_maf} & INFO/TOPMED[1]<${filter_max_af}" -s TOPMED_MAF -m+ -Ou - \
	| bcftools filter -e "INFO/TOPMED[1]>=${filter_max_af}" -s fixed_TOPMED -m+ -Ou - \
	| bcftools filter -e "INFO/GnomAD[1]>=${filter_max_maf} & INFO/GnomAD[1]<${filter_max_af}" -s GnomAD_MAF -m+ -Ou - \
	| bcftools filter -e "INFO/GnomAD[1]>=${filter_max_af}" -s fixed_GnomAD -m+ -Ou - \
	| bcftools filter -e "INFO/GnomAD_exomes[1]>=${filter_max_maf} & INFO/GnomAD_exomes[1]<${filter_max_af}" -s GnomAD_exomes_MAF -m+ -Ou - \
	| bcftools filter -e "INFO/GnomAD_exomes[1]>=${filter_max_af}" -s fixed_GnomAD_exomes -m+ -Ou - \
	| bcftools filter -e "INFO/1000Genomes[1]>=${filter_max_maf} & INFO/1000Genomes[1]<${filter_max_af}" -s 1000Genomes_MAF -m+ -Ou - \
	| bcftools filter -e "INFO/1000Genomes[1]>=${filter_max_af}" -s fixed_1000Genomes -m+ -Oz -o ${output_with_filter_vcf} -
	tabix -p vcf ${output_with_filter_vcf}
		
	bcftools view -f PASS -Ov ${output_with_filter_vcf} \
	| java -Xmx${JVM_MEM} -jar ${snpeff_tool_dir}/snpEff.jar ann -s snpEff_summary_${snpeff_anno_db}.${output_filter_string}.html ${snpeff_anno_db} \
	| bgzip -c > ${output_vcf} 
	tabix -p vcf ${output_vcf}
}

filter_on_chrM(){
	write_status "Filtering PASS VCF for chrM variants"

	local input_vcf="$1"
	local output_vcf="$2"
	local output_filter_string="$3"
	
	bcftools view -f PASS -r chrM -Ov ${input_vcf} \
	| java -Xmx${JVM_MEM} -jar ${snpeff_tool_dir}/snpEff.jar ann -s snpEff_summary_${snpeff_anno_db}.${output_filter_string}.html ${snpeff_anno_db} \
	| bgzip -c > ${output_vcf} 
	tabix -p vcf ${output_vcf}
}

filter_on_min_tumor_qual(){
	write_status "Filtering VCF for TIER min. tumor QUAL and saving PASS VCF"

	local input_vcf="$1"
	local output_vcf="$2"
	local filter_string="$3"
	local output_filter_string="$4"

	bcftools filter -e "INFO/TIER==\"HOTSPOT\" && QUAL<${filter_hotspot_min_tumor_qual}" -s ${filter_string} -m+ -Ou ${input_vcf} \
	| bcftools filter -e "INFO/TIER==\"PANEL\" && QUAL<${filter_panel_min_tumor_qual}" -s ${filter_string} -m+ -Ou - \
	| bcftools filter -e "INFO/TIER==\"HIGH_CONFIDENCE\" && QUAL<${filter_high_confidence_min_tumor_qual}" -s ${filter_string} -m+ -Ou - \
	| bcftools filter -e "INFO/TIER==\"LOW_CONFIDENCE\" && QUAL<${filter_low_confidence_min_tumor_qual}" -s ${filter_string} -m+ -Ou - \
	| bcftools view -f PASS -Ov - \
	| java -Xmx${JVM_MEM} -jar ${snpeff_tool_dir}/snpEff.jar ann -s snpEff_summary_${snpeff_anno_db}.${output_filter_string}.html ${snpeff_anno_db} \
	| bgzip -c > ${output_vcf} 
	tabix -p vcf ${output_vcf}
}

filter_on_tumor_vaf(){
	write_status "Filtering VCF for single tumor sample min. tumor QUAL and saving PASS VCF"

	local input_bcf="$1"
	local output_vcf="$2"
	local filter_string="$3"
	local output_filter_string="$4"

	bcftools filter -e "INFO/TIER==\"HOTSPOT\" && FMT/AF[0]<${hotspot_min_tumor_vaf}" -s ${filter_string} -m+ -Ou ${input_bcf} \
	| bcftools filter -e "INFO/TIER==\"PANEL\" && FMT/AF[0]<${panel_min_tumor_vaf}" -s ${filter_string} -m+ -Ou - \
	| bcftools filter -e "INFO/TIER==\"HIGH_CONFIDENCE\" && FMT/AF[0]<${high_confidence_min_tumor_vaf}" -s ${filter_string} -m+ -Ou - \
	| bcftools filter -e "INFO/TIER==\"LOW_CONFIDENCE\" && FMT/AF[0]<${low_confidence_min_tumor_vaf}" -s ${filter_string} -m+ -Ou - \
	| bcftools view -f PASS -Ov - \
	| java -Xmx${JVM_MEM} -jar ${snpeff_tool_dir}/snpEff.jar ann -s snpEff_summary_${snpeff_anno_db}.${output_filter_string}.html ${snpeff_anno_db} \
	| bgzip -c > ${output_vcf} 
	tabix -p vcf ${output_vcf}
}

filter_based_on_inclusion_string(){
	write_status "Annotating filtered PON file and saving genome-wide cutoff based PASS VCF"

	local input_vcf="$1"
	local output_vcf="$2"
	local inclusion_filter_string="$3"
	local filter_string="$4"
	local output_filter_string="$5"

	bcftools filter -i ${inclusion_filter_string} -s ${filter_string} -m+ -Ou ${input_vcf} \
	| bcftools view -f PASS -Ov - \
	| java -Xmx${JVM_MEM} -jar ${snpeff_tool_dir}/snpEff.jar ann -s snpEff_summary_${snpeff_anno_db}.${output_filter_string}.html ${snpeff_anno_db} \
	| bgzip -c > ${output_vcf} 
	tabix -p vcf ${output_vcf}
}
