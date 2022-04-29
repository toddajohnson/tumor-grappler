#!/bin/bash

RUN_DIR=${RUN_DIR}
threads=${threads}
jvm_mem=${jvm_mem}

. ${RUN_DIR}/config_files/common_config.sh

module use /usr/local/package/modulefiles/
module load java/8

export _JAVA_OPTIONS="-Xmx${jvm_mem}"

echo "Creating expected transcript counts file for Isofox"
java -jar ${ISOFOX_JAR} \
	-functions EXPECTED_TRANS_COUNTS \
	-output_dir ${isofox_dir} \
	-gene_transcripts_dir ${ensembl_cache_dir} \
	-read_length ${read_length} \
	-long_frag_limit 550 \
	-exp_rate_frag_lengths "50-0;75-0;100-0;125-0;150-0;200-0;250-0;300-0;400-0;550-0" \
	-ref_genome_version ${REFERENCE_GENOME} \
	-threads ${threads} 

echo "Finished creating expected transcript counts file for Isofox"