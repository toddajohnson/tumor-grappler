#!/bin/bash

RUN_DIR=${RUN_DIR}
OVERHANG_BP=${OVERHANG_BP}

. ${RUN_DIR}/config_files/common_config.sh

threads=${threads}

## STAR 2.7.9a, samtools and bcftools 1.12, cutadapt are in /home/tjohnson/.local
export PATH="/home/tjohnson/.local/bin:${PATH}"


REFERENCE_GENOME_FA="/home/tjohnson/reference/HMF/38/refgenomes/Homo_sapiens.GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
REFERENCE_GENOME_GTF="/home/tjohnson/reference/GRCh38/gencode.v38.annotation.gtf"
STAR_GENOME_INDICES_DIR="/home/tjohnson/reference/GRCh38/STAR-2.7.9a_genome_indices_${OVERHANG_BP}bp/GCA_000001405.15_GRCh38_no_alt_analysis_set_gencode_v38"

mkdir -p ${STAR_GENOME_INDICES_DIR}

echo "Generating genomic indices for ${OVERHANG_BP}bp overhang using STAR-2.7.9a"
echo "Genome references: ${REFERENCE_GENOME_FA}"
echo "Annotation gtf: ${REFERENCE_GENOME_GTF}"

STAR --runThreadN ${threads} \
	--runMode genomeGenerate \
	--genomeDir ${STAR_GENOME_INDICES_DIR} \
	--genomeFastaFiles ${REFERENCE_GENOME_FA} \
	--sjdbGTFfile ${REFERENCE_GENOME_GTF} \
	--sjdbOverhang ${OVERHANG_BP}

echo "Created output in ${STAR_GENOME_INDICES_DIR}"
	