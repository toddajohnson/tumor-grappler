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
	-functions EXPECTED_GC_COUNTS \
	-output_dir ${isofox_dir} \
	-ref_genome ${REFERENCE_GENOME_FA} \
	-gene_transcripts_dir ${ensembl_cache_dir} \
	-read_length ${read_length} \
	-ref_genome_version ${REFERENCE_GENOME} \
	-threads ${threads} 

echo "Finished creating expected transcript counts file for Isofox"