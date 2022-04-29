#!/bin/bash

# Variables passed from calling script
threads="${threads}"
gridss_version="${gridss_version}"
reference_genome_version="${reference_genome_version}"
run_dir="${RUN_DIR}"
jvmheap=${jvmheap}

. ${run_dir}/config_files/common_config.sh

BAM_INFO_FILE="${run_dir}/config_files/bam_info.csv"

## First line is header
let "LINENUM = ${SGE_TASK_ID} + 1"

BAM_INFO=$(sed -n "${LINENUM}p" ${BAM_INFO_FILE})
subject_id=$(echo "${BAM_INFO}" | awk -F ',' '{print $1}')
sample_name=$(echo ${BAM_INFO} | awk -F "," '{print $2}')
bam_file=$(echo ${BAM_INFO} | awk -F "," '{print $5}')

module use /usr/local/package/modulefiles/
module load java/8
module load bwa/0.7.17
module load R/4.0.2

export _JAVA_OPTIONS="-Xmx${jvmheap}"

install_dir="/home/tjohnson/tools/gridss-${gridss_version}-purple-linx"
gridss_jar="${install_dir}/gridss/gridss-${gridss_version}-gridss-jar-with-dependencies.jar"

## Need to have Perl packages in local directory available
eval $(perl -I${HOME}/perl5/lib/perl5 -Mlocal::lib)
## Needs kraken2, latest samtools, and circos in path
## Latest repeatmasker 4.1.1 needs python3 with h5py and latest RMBlast 2.10.0
## Installed and updated python 3.6.5 which is installed at:
## Installed RepeatMasker 4.1.1 using built-in curated Dfam 3.2
## samtools and bcftools 1.12 are installed now in /home/tjohnson/.local
export PATH="${install_dir}/gridss:/home/tjohnson/.local/bin:/home/tjohnson/.local/RepeatMasker:/home/tjohnson/tools/kraken2:/home/tjohnson/tools/circos-0.69-9/bin:${PATH}"
#export LD_LIBRARY_PATH="/home/tjohnson/.local/lib:${LD_LIBRARY_PATH}"
#export GRIDSS_JAR=${gridss_jar}


base_ref_dir="/home/tjohnson/reference"
	
if [[ "${reference_genome_version}" == "37" ]] ;  then
	ref_dir="${base_ref_dir}/HMF/37"
	ref_genome="${base_ref_dir}/HMF/37/refgenomes/Homo_sapiens.GRCh37/GRCh37.fa"
elif [[ "${reference_genome_version}" == "38" ]] ;  then
	ref_dir="${base_ref_dir}/HMF/38"
	ref_genome="${ref_dir}/refgenomes/Homo_sapiens.GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
fi

virusbreakenddb="${base_ref_dir}/VIRUSBreakend/virusbreakenddb"

sample_bam="${run_dir}/result/bam_final/${bam_file}"

main_output_dir="${run_dir}/result/VIRUSBreakend"
subject_output_dir="${main_output_dir}/${subject_id}"
sample_vbe_output_dir="${subject_output_dir}/${sample_name}"

virusbreakend_vcf=${sample_vbe_output_dir}/${sample_name}.vbe.vcf

mkdir -p $sample_vbe_output_dir

cd ${sample_vbe_output_dir}

write_status() { # Before logging initialised
	echo "$(date): $1" 1>&2
}

assert_file_exists() {
	if [[ ! -f "$1" ]] ; then
		write_status "File $1 not found. Specify using the command line argument --$2"
		exit $EX_NOINPUT
	fi
}
assert_directory_exists() {
	if [[ ! -d "$1" ]] ; then
		write_status "Directory $1 not found. Specify using the command line argument --$2"
		exit $EX_NOINPUT
	fi
}
assert_directory_exists $install_dir/gridss "install_dir"
assert_directory_exists $install_dir/hmftools "install_dir"
assert_file_exists $install_dir/gridss/gridss "install_dir"
assert_file_exists $install_dir/gridss/libgridss.R "install_dir"


write_status "sample_vbe_output_dir=$sample_vbe_output_dir"
write_status "virusbreakend_vcf=$virusbreakend_vcf"
write_status "install_dir=$install_dir"
write_status "bam=$sample_bam"
write_status "threads=$threads"
write_status "sample=$sample_name"
write_status "ref_genome=$ref_genome"
write_status "virusbreakenddb=$virusbreakenddb"
write_status "gridss_jar=$gridss_jar"


write_status "Running GRIDSS VIRUSBreakend"
virusbreakend \
	-j ${gridss_jar} \
	--kraken2db ${virusbreakenddb} \
	--output ${virusbreakend_vcf} \
	--reference ${ref_genome} \
	-w ${sample_vbe_output_dir} \
	-t $threads \
	${sample_bam}

if [[ ! -f $virusbreakend_vcf ]] ; then
	write_status "Error creating $virusbreakend_vcf. Aborting" 1>&2
	exit 1
fi
