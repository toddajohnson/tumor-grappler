#!/bin/bash

run_dir="${RUN_DIR}"
threads=${threads}
pon_output_file=${pon_output_file}

. ${run_dir}/config_files/common_config.sh

sage_vcf_dir="${run_dir}/result/vcf/SAGE/PON_mode"

target_pon_dir="${run_dir}/result/PON/SAGE-${sage_version}"
mkdir -p ${target_pon_dir}

cd ${sage_vcf_dir}

java -Xmx32g \
	-cp ${SAGE_JAR} \
	com.hartwig.hmftools.sage.pon.PonApplication \
	-in ${sage_vcf_dir} \
	-out ${target_pon_dir}/${pon_output_file} \
	-threads ${threads}
