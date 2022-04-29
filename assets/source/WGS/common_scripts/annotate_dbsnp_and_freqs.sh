#!/bin/bash

mkdir -p tmp
ann_tmp_bcf1=$(mktemp tmp/tmpXXXXX.bcf)

bcftools annotate --rename-chrs ${ucsc_to_ncbi_chr_conv_file} -Ob -o ${ann_tmp_bcf1} ${ann_input_vcf}
bcftools index ${tmp_bcf1}
		
		bcftools annotate -a ${dbsnp_common_vcf} -c COMMON -Ob -o ${tmp_bcf2} ${tmp_bcf1}
		bcftools index ${tmp_bcf2}
		
		bcftools annotate --rename-chrs ${ncbi_to_ucsc_chr_conv_file} -Ou ${tmp_bcf2} \
		| bcftools filter -e 'COMMON!=1' -s COMMON -m+ -Oz -o ${tmp_vcf1} - 
		tabix -p vcf ${tmp_vcf1}
		
		table_annovar.pl -thread ${threads} ${tmp_vcf1} ${humandb_path} -buildver hg38 -out ${base_output_dir}/${subject_id}.${call_target}.sage.PON_filtered -remove -protocol avsnp150,gnomad30_genome -operation f,f -nastring . -vcfinput -polish
		bgzip ${base_output_dir}/${subject_id}.${call_target}.sage.PON_filtered.hg38_multianno.vcf
		tabix -p vcf ${base_output_dir}/${subject_id}.${call_target}.sage.PON_filtered.hg38_multianno.vcf.gz
