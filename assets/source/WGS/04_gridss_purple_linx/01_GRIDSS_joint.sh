#!/bin/bash

GPL_PREFIX=$(basename "$PWD")

# Variables passed from calling script
run_dir="${RUN_DIR}"
threads="${PGL_THREAD_NUM}"
jvmheap="${PGL_JVM_MEM}"
job_nodes="${JOB_NODES}"
gridss_version="${GRIDSS_VERSION}"
use_gridss_suffix=${USE_GRIDSS_SUFFIX}
steps="${STEPS}"

subject_id="${SUBJECT_ID}"
joint_sample_labels=$(echo "${JOINT_SAMPLE_LABELS}" | tr " " ",")
joint_bams="${JOINT_BAMS}"

. ${run_dir}/config_files/common_config.sh


## Change of script to parallelize assembly across nodes
## Pass parameters job_index
let "job_index = ${SGE_TASK_ID} - 1"

if [[ -z ${threads} ]]; then
	echo "The number of threads was not passed."
fi
if [[ -z ${jvmheap} ]]; then
	echo "The maximum JVM heap size was not set."
fi
if [[ -z ${job_nodes} ]]; then
	echo "The number of nodes to use on the job was not passed."
fi
if [[ -z ${gridss_version} ]]; then
	echo "The GRIDSS version was not set."
fi
if [[ -z ${steps} ]]; then
	echo "The GRIDSS steps to perform were not set."
fi

if [[ -z ${run_dir} ]]; then
	echo "The run directory was not set."
fi

if [[ -z ${subject_id} ]]; then
	echo "There was no subject ID passed."
fi

if [[ -z ${joint_sample_labels} ]]; then
	echo "The joint sample labels were not passed."
fi

if [[ -z ${joint_bams} ]]; then
	echo "The joint bams list was not passed."
fi

module use /usr/local/package/modulefiles/
module load java/8
module load bwa/0.7.17
##module load repeatmasker/4.1.0
##module load blast+/2.9.0
##module load R/4.0.2

export _JAVA_OPTIONS="-Xmx${jvmheap}"

## Need to have Perl packages in local directory available
eval $(perl -I${HOME}/perl5/lib/perl5 -Mlocal::lib)
## Needs kraken2, latest samtools, and circos in path
## Latest repeatmasker 4.1.1 needs python3 with h5py and latest RMBlast 2.10.0
## Installed and updated python 3.6.5 which is installed at:
## Installed RepeatMasker 4.1.1 using built-in curated Dfam 3.2
## samtools and bcftools 1.12 are installed now in /home/tjohnson/.local
export PATH="/home/tjohnson/.local/bin:/home/tjohnson/.local/RepeatMasker:/home/tjohnson/tools/kraken2:/home/tjohnson/tools/circos-0.69-9/bin:${PATH}"

## load R from module after PATH to block R 4.1.0 that is installed in ~/.local/bin
module load R/4.0.2

install_dir="/home/tjohnson/tools/gridss-${gridss_version}-purple-linx"
export gridss_jar=$install_dir/gridss/gridss-${gridss_version}-gridss-jar-with-dependencies.jar

reference_genome_version="HG38"

base_ref_dir="/home/tjohnson/reference"
ref_dir="${base_ref_dir}/HMF/38"
ref_genome=${ref_dir}/refgenomes/Homo_sapiens.GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
blacklist=${ref_dir}/dbs/gridss/ENCFF001TDO.38.bed
gridss_properties=${ref_dir}/dbs/gridss/gridss.properties
rlib=${ref_dir}/rlib/
virusbreakenddb="${base_ref_dir}/VIRUSBreakend/virusbreakenddb"

viralreference="${ref_dir}/refgenomes/human_virus/human_virus.fa"

joint_sample_name=${subject_id}
gridss_main_output_dir=${run_dir}/result/${GPL_PREFIX}
subject_output_dir=${gridss_main_output_dir}/${subject_id}
subject_gridss_output_dir=${subject_output_dir}/gridss

assembly_bam=${subject_gridss_output_dir}/${joint_sample_name}.assembly.bam
gridss_driver_vcf=${subject_gridss_output_dir}/${joint_sample_name}.gridss.driver.vcf.gz
gridss_kraken2_vcf=${subject_gridss_output_dir}/${joint_sample_name}.gridss.kraken2.vcf.gz

## added rm_annotated and human viral reference annotations 7/7/2021
gridss_rm_annotated_vcf=${subject_gridss_output_dir}/${joint_sample_name}.gridss.rm.vcf.gz
##
gridss_unfiltered_vcf=${subject_gridss_output_dir}/${joint_sample_name}.gridss.unfiltered.vcf.gz

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

if [[ "${use_gridss_suffix}" == "true" ]] ; then
	assert_file_exists $install_dir/gridss/gridss.sh "install_dir"
	gridss_cmd=$install_dir/gridss/gridss.sh
	kraken2_annotation_cmd=$install_dir/gridss/gridss_annotate_vcf_kraken2.sh
	repeatmasker_annotation_cmd=$install_dir/gridss/gridss_annotate_vcf_repeatmasker.sh
else
	assert_file_exists $install_dir/gridss/gridss "install_dir"
	gridss_cmd=$install_dir/gridss/gridss
	kraken2_annotation_cmd=$install_dir/gridss/gridss_annotate_vcf_kraken2
	repeatmasker_annotation_cmd=$install_dir/gridss/gridss_annotate_vcf_repeatmasker
fi

assert_file_exists $install_dir/gridss/libgridss.R "install_dir"

cd ${subject_gridss_output_dir}

mkdir -p $subject_output_dir/logs $subject_output_dir/gridss

write_status "subject_gridss_output_dir=$subject_gridss_output_dir"
write_status "gridss_driver_vcf=$gridss_driver_vcf"
write_status "ref_dir=$ref_dir"
write_status "install_dir=$install_dir"
write_status "joint_sample_labels=$joint_sample_labels"
write_status "joint_bams=$joint_bams"
write_status "threads=$threads"
write_status "sample=$subject_id"
write_status "jvmheap=$jvmheap"
write_status "rlib=$rlib"
write_status "ref_genome=$ref_genome"
write_status "virusbreakenddb=$virusbreakenddb"
write_status "steps=$steps"
write_status "gridss_jar=$gridss_jar"

if [[ ! -f $gridss_driver_vcf ]]; then
	write_status "Running GRIDSS"
	${gridss_cmd} \
		-o ${gridss_driver_vcf} \
		-a $assembly_bam \
		-w ${subject_gridss_output_dir} \
		-r ${ref_genome} \
		-j ${gridss_jar} \
		-t $threads \
		-b ${blacklist} \
		-c ${gridss_properties} \
		--jobindex ${job_index} \
		--jobnodes ${job_nodes} \
		--jvmheap $jvmheap \
		--steps $steps \
		--labels ${joint_sample_labels} \
		${joint_bams}
		if [[ ! -f $gridss_driver_vcf ]] ; then
			write_status  "Error creating $gridss_driver_vcf. Aborting" 1>&2
			exit 1
		fi
elif [[ -f $gridss_driver_vcf ]] && [[ $steps=="recall" ]]; then
	write_status "Renaming previous VCF file"
	mv $gridss_driver_vcf ${subject_gridss_output_dir}/${joint_sample_name}.gridss.driver.previous.vcf.gz
	mv $gridss_kraken2_vcf ${subject_gridss_output_dir}/${joint_sample_name}.gridss.kraken2.previous.vcf.gz
	mv $gridss_rm_annotated_vcf ${subject_gridss_output_dir}/${joint_sample_name}.gridss.rm.previous.vcf.gz
	mv $gridss_unfiltered_vcf ${subject_gridss_output_dir}/${joint_sample_name}.gridss.unfiltered.previous.vcf.gz

	#rm ${gridss_driver_vcf}
	#rm ${$gridss_kraken2_vcf}
	#rm ${gridss_rm_annotated_vcf}
	#rm ${$gridss_unfiltered_vcf}
	
	write_status "Running GRIDSS to recall SVs"
	${gridss_cmd} \
		-o ${gridss_driver_vcf} \
		-a $assembly_bam \
		-w ${subject_gridss_output_dir} \
		-r ${ref_genome} \
		-j ${gridss_jar} \
		-t $threads \
		-b ${blacklist} \
		-c ${gridss_properties} \
		--jobindex ${job_index} \
		--jobnodes ${job_nodes} \
		--jvmheap $jvmheap \
		--steps call \
		--labels ${joint_sample_labels} \
		${joint_bams}
		if [[ ! -f $gridss_driver_vcf ]] ; then
			write_status  "Error creating $gridss_driver_vcf. Aborting" 1>&2
			exit 1
		fi
else
	write_status "Skipping GRIDSS - ${gridss_driver_vcf} exists"
fi

unset _JAVA_OPTIONS

if [[ "$steps" == "annotation" ]] || [[ "$steps" == "all" ]] || [[ "$steps" == "recall" ]]; then
	annotation_log_dir=${subject_gridss_output_dir}/logs_annotation
	mkdir -p ${annotation_log_dir}
	
	log_prefix=${annotation_log_dir}/${tumour_sample}_$(date +%Y%m%d_%H%M%S).$HOSTNAME.$$
	
	if [[ -s $gridss_kraken2_vcf ]] ; then
		write_status "Skipping GRIDSS Kraken2 Viral Annotation = ${gridss_kraken2_vcf} exists"
	else
		write_status "Running GRIDSS Kraken2 viral Annotation"
		${kraken2_annotation_cmd} \
			-o ${gridss_kraken2_vcf} \
			-j ${gridss_jar} \
			-t $threads \
			--kraken2db ${virusbreakenddb} \
			${gridss_driver_vcf} \
			2>&1 | tee $log_prefix.gridss.AnnotateInsertedSequence.log
		if [[ ! -s $gridss_kraken2_vcf ]] ; then
			write_status "Error creating $gridss_kraken2_vcf. Aborting" 1>&2
			exit 1
		fi	
	fi
	
	if [[ -s $gridss_rm_annotated_vcf ]] ; then
		write_status "Skipping GRIDSS RepeatMasker Annotations = ${gridss_rm_annotated_vcf} exists"
	else
		write_status "Running GRIDSS RepeatMasker Annotation"
		${repeatmasker_annotation_cmd} \
			-o ${gridss_rm_annotated_vcf} \
			-w ${subject_gridss_output_dir} \
			-j ${gridss_jar} \
			-t $threads \
			${gridss_kraken2_vcf} \
			2>&1 | tee $log_prefix.gridss.AnnotateInsertedSequence.log
		if [[ ! -s $gridss_rm_annotated_vcf ]] ; then
			write_status "Error creating $gridss_rm_annotated_vcf. Aborting" 1>&2
			exit 1
		fi		
	fi
	
	if [[ -s $gridss_unfiltered_vcf ]] ; then
		write_status "Skipping GRIDSS viral reference annotation = ${gridss_unfiltered_vcf} exists"
	else
		write_status "Running GRIDSS Annotations with viral reference"
		java -Xmx8G -Dsamjdk.create_index=true \
			-Dsamjdk.use_async_io_read_samtools=true \
			-Dsamjdk.use_async_io_write_samtools=true \
			-Dsamjdk.use_async_io_write_tribble=true \
			-Dsamjdk.buffer_size=4194304 \
			-cp ${gridss_jar} gridss.AnnotateInsertedSequence \
			REFERENCE_SEQUENCE=${viralreference} \
			INPUT=${gridss_rm_annotated_vcf} \
			OUTPUT=${gridss_unfiltered_vcf} \
			ALIGNMENT=APPEND \
			WORKER_THREADS=${threads} \
			 2>&1 | tee $log_prefix.gridss.AnnotateInsertedSequence.log
		if [[ ! -s $gridss_unfiltered_vcf ]] ; then
			write_status "Error creating $gridss_unfiltered_vcf. Aborting" 1>&2
			exit 1
		fi
	fi
		
fi
