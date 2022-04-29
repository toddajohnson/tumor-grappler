#!/bin/bash

run_dir="${RUN_DIR}"
threads=${threads}
breakpoint_filename=${breakpoint_filename}
breakend_filename=${breakend_filename}

. ${run_dir}/config_files/common_config.sh

gridss_vcf_dir="${run_dir}/result/vcf/GRIDSS"

target_pon_dir="${run_dir}/result/PON/GRIDSS-${gridss_version}"
mkdir -p ${target_pon_dir}

cd ${gridss_vcf_dir}

java -Xmx8g \
	-cp ${GRIDSS_JAR} \
	gridss.GeneratePonBedpe \
	$(ls -1 *.vcf.gz | awk ' { print "I=" $0 }') \
	O=${target_pon_dir}/${breakpoint_filename} \
	SBO=${target_pon_dir}/${breakend_filename} \
	NO=0 \
	THREADS=${threads} \
	INCLUDE_IMPRECISE_CALLS=false \
	REFERENCE_SEQUENCE=$ref_genome