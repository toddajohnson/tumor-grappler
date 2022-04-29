#!/bin/bash

jvm_mem=${jvm_mem}

module use /usr/local/package/modulefiles/
module load java/8

GRIDSS_VERSION=2.11.1
LINX_VERSION=1.15

install_dir="/home/tjohnson/tools/gridss-${GRIDSS_VERSION}-purple-linx"
export LINX_JAR=$install_dir/hmftools/sv-linx_v${LINX_VERSION}.jar

java -Xmx${jvm_mem} -cp ${LINX_JAR} com.hartwig.hmftools.linx.gene.GenerateEnsemblDataCache \
    -ensembl_db "mysql://ensembldb.ensembl.org:3306/homo_sapiens_core_104_38" -ensembl_user "anonymous" -ensembl_pass "" \
    -output_dir /home/tjohnson/reference/HMF_38/dbs/linx -ref_genome_version 38
   