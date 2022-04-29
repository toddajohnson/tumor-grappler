#!/bin/bash
#
# Stand-alone GRIPSS-PURPLE-LINX pipeline
#
## Modified to allow tumor-normal pair analysis in GRIPSS: 2/14/2021 (Todd A. Johnson)
## Modified to use pre-run GRIDSS joint SV output (no GRIDSS run in current script): 3/21/2021 (Todd A. Johnson)
## Modified to produce driver calls (provide germline & somatic vcfs, Actionable transcript bed files, etc.
#
# Example: ./gripss-purple-linx.sh -n /data/COLO829R_dedup.realigned.bam -t /data/COLO829T_dedup.realigned.bam -v /data/colo829snv.vcf.gz -s colo829 -v /data/COLO829v003T.somatic_caller_post_processed.vcf.gz
# docker run  gridss/gridss-purple-linx
EX_USAGE=64
EX_NOINPUT=66
EX_CANTCREAT=73
EX_CONFIG=78
set -o errexit -o pipefail -o noclobber -o nounset
! getopt --test > /dev/null 
if [[ ${PIPESTATUS[0]} -ne 4 ]]; then
	echo '`getopt --test` failed in this environment.'
	exit $EX_CONFIG
fi

run_dir=/data
ref_dir=/refdata
install_dir=/opt/
tumour_bam=""
normal_bam=""
snvvcf=""
threads=$(nproc)
sample=""
normal_sample=""
tumour_sample=""
jvmheap="25g"
ref_genome_version="38"
purple_args=""
amber_args=""
cobalt_args=""
linx_args=""
gripss_args=""
gridss_args=""
tumour_only="false"

validation_stringency="STRICT"
usage_msg="Usage: gripss-purple-linx.sh

Required command line arguments:
	--tumour_bam: tumour BAM file
	--normal_bam: matched normal BAM file
	--sample: sample name
Optional parameters:
	--output_dir: output directory. (/data)
	--ref_genome_version: reference genome. 37 or 38 ($ref_genome_version)
	--ref_dir: path to decompressed Hartwig reference data package. ($ref_dir)
	--snvvcf: A somatic SNV VCF with the AD genotype field populated.
	--nosnvvcf: Indicates a somatic SNV VCF will not be supplied. This will reduce the accuracy of PURPLE ASCN.
	--threads: number of threads to use. ($threads)
	--install_dir: root directory of gripss-purple-linx release ($install_dir)
	--normal_sample: sample name of matched normal ({sample}_N) 
	--tumour_sample: sample name of tumour. Must match the somatic \$snvvcf sample name. ({sample}_T) 
	--jvmheap: maximum java heap size for high-memory steps ($jvmheap)
	--gridss_args: additional arguments to GRIDSS
	--gripss_args: additional arguments to GRIPSS
	--amber_args: additional arguments to AMBER
	--cobalt_args: additional arguments to COBALT
	--purple_args: additional arguments to PURPLE
	--linx_args: additional arguments to LINX
	--validation_stringency: htsjdk validation_stringency ($validation_stringency)
	--help: print this message and exit
"
usage() {
	echo "$usage_msg" 1>&2
	exit $EX_USAGE
}

OPTIONS=v:o:t:n:s:r:b:h
LONGOPTS=snvvcf:,nosnvvcf,output_dir:,tumour_bam:,normal_bam:,sample:,threads:,jvmheap:,ref_dir:,ref_genome_version:,normal_sample:,tumour_sample:,rundir:,install_dir:,help,gridss_args:,gripss_args:,amber_args:,cobalt_args:,purple_args:,linx_args:
! PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
	# e.g. return value is 1
	#  then getopt has complained about wrong arguments to stdout
	exit 2
fi
eval set -- "$PARSED"
while true; do
	case "$1" in
		--rundir)
			run_dir="$2"
			shift 2
			;;
		-v|--snvvcf)
			snvvcf="$2"
			shift 2
			;;
		--nosnvvcf)
			snvvcf="nosnvvcf"
			shift 1
			;;
		--ref_genome_version)
			ref_genome_version="$2"
			shift 2
			;;
		-n|--normal_bam)
			normal_bam="$2"
			shift 2
			;;
		-o|--output_dir)
			run_dir="$2"
			shift 2
			;;
		-t|--tumour_bam)
			tumour_bam="$2"
			shift 2
			;;
		-s|--sample)
			sample="$2"
			shift 2
			;;
		--normal_sample)
			normal_sample="$2"
			shift 2
			;;
		--tumour_sample)
			tumour_sample="$2"
			shift 2
			;;
		--threads)
			printf -v threads '%d\n' "$2" 2>/dev/null
			printf -v threads '%d' "$2" 2>/dev/null
			shift 2
			;;
		--jvmheap)
			jvmheap="$2"
			shift 2
			;;
		--install_dir)
			install_dir="$2"
			shift 2
			;;
		--ref_dir)
			ref_dir="$2"
			shift 2
			;;
		-h|--help)
			usage
			exit 1
			;;
		--gridss_args)
			gridss_args="$2"
			shift 2
			;;
		--gripss_args)
			gripss_args="$2"
			shift 2
			;;
		--amber_args)
			amber_args="$2"
			shift 2
			;;
		--cobalt_args)
			cobalt_args="$2"
			shift 2
			;;
		--purple_args)
			purple_args="$2"
			shift 2
			;;
		--linx_args)
			linx_args="$2"
			shift 2
			;;
		--)
			shift
			break
			;;
		*)
			echo "Command line parsing error ($1)"
			echo "$@"
			exit 3
			;;
	esac
done
write_status() { # Before logging initialised
	echo "$(date): $1" 1>&2
}
# $1: variable containing filename
# $2: command line argument name
assert_file_exists() {
	if [[ ! -s "$1" ]] ; then
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

#assert_directory_exists $install_dir/gridss "install_dir"
assert_directory_exists $install_dir/hmftools "install_dir"
#assert_directory_exists $install_dir/gripss-purple-linx "install_dir"
#assert_file_exists $install_dir/gridss/gridss "install_dir"
#assert_file_exists $install_dir/gridss/libgridss.R "install_dir"

basedir=$(echo "$run_dir" | awk -F "/" '{print $6}')
studyname=$(echo "$basedir" | awk -F "_" '{print $1}')

viralreference=refgenomes/human_virus/human_virus.fa
blacklist=dbs/gridss/ENCFF001TDO.38.bed
repeatmasker=dbs/repeatmasker/38.fa.out.bed
#bafsnps=dbs/germline_het_pon/GermlineHetPon.vcf.gz
#bafsnps=dbs/germline_het_pon/BHD_het_sites_in_dbs.vcf.gz
bafsnps=dbs/germline_het_pon/het_sites_from_91_reference_samples.vcf.gz
gcprofile=dbs/gc/GC_profile.1000bp.38.cnp
gridss_properties=dbs/gridss/gridss.properties
## Replaced with latest HMF fusion files on 20211215
breakpoint_hotspot=dbs/knowledgebases/known_fusions.38.bedpe
#breakpoint_hotspot=dbs/knowledgebases/known_fusions.38_v3.bedpe
known_fusion_csv=dbs/knowledgebases/known_fusion_data.38.csv
#known_fusion_csv=dbs/knowledgebases/known_fusion_data.38_v3.csv
breakend_pon=dbs/gridss_pon/gridss_pon_single_breakend.JP.149x.bed
breakpoint_pon=dbs/gridss_pon/gridss_pon_breakpoint.JP.149x.bedpe
viral_hosts_csv=dbs/knowledgebases/viral_host_ref.csv
fragile_sites=dbs/knowledgebases/fragile_sites_hmf.38.csv
line_elements=dbs/knowledgebases/line_elements.38.csv
replication_origins=dbs/knowledgebases/heli_rep_origins.38.bed
tumor_only_diploid_bed=dbs/knowledgebases/DiploidRegions.38.bed
#ensembl_data_dir=dbs/ensembl_data_cache
## Updated to ensembl 104, Gencode 38
ensembl_data_dir=dbs/linx

## Updated base DriverGenePanel.38.tsv to current HMF version on 20211215
if [[ "$studyname" == "BHD" ]]; then
	driver_gene_panel=dbs/knowledgebases/DriverGenePanel.38.BHD.tsv
elif [[ "$studyname" == "BTC" ]]; then
	driver_gene_panel=dbs/knowledgebases/DriverGenePanel.38.BTC.tsv
elif [[ "$studyname" == "KeioOrganoid" ]]; then
	driver_gene_panel=dbs/knowledgebases/DriverGenePanel.38.BTC.tsv
elif [[ "$studyname" == "KakimiLab" ]]; then
	if [[ "$sample" == "GKSA003" ]]; then
		driver_gene_panel=dbs/knowledgebases/DriverGenePanel.38.GKSA003.tsv
	else
		driver_gene_panel=dbs/knowledgebases/DriverGenePanel.38.tsv
	fi
else
	driver_gene_panel=dbs/knowledgebases/DriverGenePanel.38.tsv
fi

known_hotspots_vcf=dbs/knowledgebases/KnownHotspots.vcf.gz
known_germline_hotspots_vcf=dbs/sage/38/KnownHotspots.germline.38.vcf.gz
known_somatic_hotspots_vcf=dbs/sage/38/KnownHotspots.somatic.38.vcf.gz
rlib=rlib/
case "$ref_genome_version" in
	"37" | "V37" | "HG37")
		ref_genome=refgenomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta
		"37"
		;;
	"38" | "hg38" | "HG38" | "V38")
		ref_genome=refgenomes/Homo_sapiens.GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
		;;
	*)
		write_status "Invalid reference genome version: $ref_genome_version"
		exit $EX_CONFIG
		;;
esac
write_status "Running reference genome version: $ref_genome_version"

rlib=$ref_dir/$rlib
ref_genome=$ref_dir/$ref_genome
viralreference=$ref_dir/$viralreference
blacklist=$ref_dir/$blacklist
repeatmasker=$ref_dir/$repeatmasker
gridss_properties=$ref_dir/$gridss_properties
bafsnps=$ref_dir/$bafsnps
gcprofile=$ref_dir/$gcprofile
breakpoint_hotspot=$ref_dir/$breakpoint_hotspot
breakend_pon=$ref_dir/$breakend_pon
breakpoint_pon=$ref_dir/$breakpoint_pon
viral_hosts_csv=$ref_dir/$viral_hosts_csv
known_fusion_csv=$ref_dir/$known_fusion_csv
fragile_sites=$ref_dir/$fragile_sites
line_elements=$ref_dir/$line_elements
replication_origins=$ref_dir/$replication_origins
ensembl_data_dir=$ref_dir/$ensembl_data_dir
driver_gene_panel=$ref_dir/$driver_gene_panel
known_hotspots_vcf=$ref_dir/${known_hotspots_vcf}
known_germline_hotspots_vcf=$ref_dir/${known_germline_hotspots_vcf}
known_somatic_hotspots_vcf=$ref_dir/${known_somatic_hotspots_vcf}
tumor_only_diploid_bed=$ref_dir/${tumor_only_diploid_bed}

if [[ "$normal_bam" == "" ]] && [[ "$normal_sample" == "" ]] ; then
	tumour_only="true"
else
	tumour_only="false"
fi

write_status "somatic vcf: $snvvcf"

if [[ "$snvvcf" == "nosnvvcf" ]] ; then
	write_status "No somatic SNV VCF supplied."
elif [[ ! -s "$snvvcf" ]] ; then
	write_status "Missing somatic SNV VCF. A SNV VCF with the AD genotype field populated is required."
	write_status "Use the script for generating this VCF with strelka if you have not already generated a compatible VCF."
	exit $EX_NOINPUT
fi
if [[ "${IGNORE_BAMS}" == "false" ]]; then
	if [[ ! -s "$tumour_bam" ]] ; then
		write_status "Missing tumour BAM: $tumour_bam"
		exit $EX_NOINPUT
	fi
fi
if [[ "${IGNORE_BAMS}" == "false" ]]; then
	if [[ "$tumour_only" != "true" ]] ; then
		if [[ ! -s "$normal_bam" ]] ; then
			write_status "Missing normal BAM"
			exit $EX_NOINPUT
		fi
	fi
fi
mkdir -p "$run_dir"
if [[ ! -d "$run_dir" ]] ; then
	write_status "Unable to create $run_dir"
	exit $EX_CANTCREAT
fi
if [[ ! -d "$ref_dir" ]] ; then
	write_status "Could not find reference data directory $ref_dir"
	exit $EX_NOINPUT
fi
if [[ ! -f "$ref_genome" ]] ; then
	write_status "Missing reference genome $ref_genome - specify with -r "
	exit $EX_NOINPUT
fi
if [[ -z "$sample" ]] ; then
	sample=$(basename $tumour_bam .bam)
fi
if [[ "$threads" -lt 1 ]] ; then
	write_status "Illegal thread count: $threads"
	exit $EX_CONFIG
fi
joint_sample_name=$sample
if [[ "$tumour_only" == "true" ]] ; then
	if [[ -z "$tumour_sample" ]] ; then
		tumour_sample=${sample}
	fi
else
	if [[ -z "$normal_sample" ]] ; then
		normal_sample=${sample}R
	fi
	if [[ -z "$tumour_sample" ]] ; then
		tumour_sample=${sample}T
	fi
fi
export R_LIBS="$rlib:${R_LIBS:-}"
base_path=$(dirname $(readlink $0 || echo $0))

### Find the jars
find_jar() {
	env_name=$1
	if [[ -f "${!env_name:-}" ]] ; then
		echo "${!env_name}"
	else
		write_status "Unable to find $2 jar. Specify using the environment variant $env_name"
		exit $EX_CONFIG
	fi
}

gridss_jar=$(find_jar GRIDSS_JAR gridss)
gripss_jar=$(find_jar GRIPSS_JAR gripss)
amber_jar=$(find_jar AMBER_JAR amber)
cobalt_jar=$(find_jar COBALT_JAR cobalt)
purple_jar=$(find_jar PURPLE_JAR purple)
linx_jar=$(find_jar LINX_JAR sv-linx)

for program in bwa samtools circos Rscript java ; do
	if ! which $program > /dev/null ; then
		write_status "Missing required dependency $program. $program must be on PATH"
		exit $EX_CONFIG
	fi
done
for rpackage in tidyverse devtools assertthat testthat NMF stringdist stringr argparser R.cache "copynumber" StructuralVariantAnnotation "VariantAnnotation" "rtracklayer" "BSgenome" "org.Hs.eg.db" ; do
	if ! Rscript -e "installed.packages()" | grep $rpackage > /dev/null ; then
		write_status "Missing R package $rpackage"
		exit $EX_CONFIG
	fi
done

#if ! java -Xms$jvmheap -cp $gridss_jar gridss.Echo ; then
#	write_status "Failure invoking java with --jvmheap parameter of \"$jvmheap\". Specify a JVM heap size (e.g. \"20g\") that is valid for this machine."
#	exit $EX_CONFIG
#fi

if [[ ! -s $ref_genome.bwt ]] ; then
	write_status "Missing bwa index for $ref_genome. Creating (this is a once-off initialisation step)"
	bwa index $ref_genome
fi

if [[ ! -s $ref_genome.bwt ]] ; then
	write_status "bwa index for $ref_genome not found."
	write_status "If you are running in a docker container, make sure refdata has been mounted read-write."
	exit $EX_NOINPUT
fi

mkdir -p $run_dir/logs $run_dir/gridss $run_dir/gripss $run_dir/amber $run_dir/purple
log_prefix=${run_dir}/logs/${tumour_sample}_$(date +%Y%m%d_%H%M%S).$HOSTNAME.$$

jvm_args=" \
	-Dsamjdk.reference_fasta=$ref_genome \
	-Dsamjdk.use_async_io_read_samtools=true \
	-Dsamjdk.use_async_io_write_samtools=true \
	-Dsamjdk.use_async_io_write_tribble=true \
	-Dsamjdk.buffer_size=4194304 \
	-Dsamjdk.async_io_read_threads=$threads"

write_status "studyname=$studyname"
write_status "run_dir=$run_dir"
write_status "ref_dir=$ref_dir"
write_status "install_dir=$install_dir"
write_status "driver_gene_panel=$driver_gene_panel"
write_status "tumour_bam=$tumour_bam"
write_status "normal_bam=$normal_bam"
write_status "snvvcf=$snvvcf"
write_status "threads=$threads"
write_status "sample=$sample"
write_status "normal_sample=$normal_sample"
write_status "tumour_sample=$tumour_sample"
write_status "jvmheap=$jvmheap"
write_status "ref_genome_version=$ref_genome_version"
write_status "rlib=$rlib"
write_status "ref_genome=$ref_genome"

gridss_dir=$run_dir/gridss
assembly_bam=$gridss_dir/$sample.assembly.bam
gridss_unfiltered_vcf=$gridss_dir/${sample}.gridss.unfiltered.vcf.gz


gripss_dir=$run_dir/gripss

if [[ "$gripss_args" == "gripss.somatic" ]] ; then
	gripss_somatic_vcf=$gripss_dir/${tumour_sample}.gripss.somatic.vcf.gz
	gripss_somatic_filtered_vcf=$gripss_dir/${tumour_sample}.gripss.somatic.filtered.vcf.gz
	
	if [[ "$tumour_only" == "false" ]] ; then
		write_status "Running GRIPSS in tumor-reference mode"
		gripss_args="-reference $normal_sample"
	elif [[ "$tumour_only" == "true" ]] ; then
		write_status "Running GRIPSS in tumor-only mode"
		gripss_args=""
	fi
	
	java -Xmx${jvmheap} -cp ${gripss_jar} com.hartwig.hmftools.gripss.GripssApplicationKt \
		-ref_genome ${ref_genome} \
		-breakpoint_hotspot ${breakpoint_hotspot} \
		-breakend_pon ${breakend_pon} \
		-breakpoint_pon ${breakpoint_pon} \
		-input_vcf ${gridss_unfiltered_vcf} \
		-output_vcf ${gripss_somatic_vcf} \
		-tumor ${tumour_sample} \
		$gripss_args 2>&1 | tee $log_prefix.gripss.log

	if [[ ! -s $gripss_somatic_vcf ]] ; then
		write_status "Error creating $gripss_somatic_vcf. Aborting"
		exit 1
	fi
	
	java -Xmx${jvmheap} -cp ${gripss_jar} com.hartwig.hmftools.gripss.GripssHardFilterApplicationKt \
		-input_vcf ${gripss_somatic_vcf} \
		-output_vcf ${gripss_somatic_filtered_vcf} \
		$gripss_args 2>&1 | tee $log_prefix.gripss.log

	if [[ ! -s ${gripss_somatic_filtered_vcf} ]] ; then
		write_status "Error creating ${gripss_somatic_filtered_vcf} - aborting"
		exit 1
	fi

elif [[ "$gripss_args" == "gripss" ]] ; then
	gripss_somatic_vcf=$gripss_dir/${tumour_sample}.gripss.vcf.gz
	gripss_somatic_filtered_vcf=$gripss_dir/${tumour_sample}.gripss.filtered.vcf.gz
	
	if [[ "$tumour_only" == "false" ]] ; then
		write_status "Running GRIPSS in tumor-reference mode"
		gripss_args="-reference $normal_sample"
	elif [[ "$tumour_only" == "true" ]] ; then
		write_status "Running GRIPSS in tumor-only mode"
		gripss_args=""
	fi
	
	# Replaced pon parameters to match GRIPSS 1.12 usage
	# Removed hardfiltering step as new GRIPSS does it in one call
	java -Xmx${jvmheap} \
		-jar ${gripss_jar} \
		-ref_genome ${ref_genome} \
		-pon_sgl_file ${breakend_pon} \
		-pon_sv_file ${breakpoint_pon} \
		-known_hotspot_file ${breakpoint_hotspot} \
		-vcf ${gridss_unfiltered_vcf} \
		-output_dir ${gripss_dir} \
		-sample ${tumour_sample} \
		$gripss_args 2>&1 | tee $log_prefix.gripss.log

	if [[ ! -s $gripss_somatic_vcf ]] ; then
		write_status "Error creating $gripss_somatic_vcf. Aborting"
		exit 1
	fi

	if [[ ! -s ${gripss_somatic_filtered_vcf} ]] ; then
		write_status "Error creating ${gripss_somatic_filtered_vcf} - aborting"
		exit 1
	fi

else
	write_status "Skipping GRIPSS - ${gripss_somatic_vcf} exists"
fi

## amber output for multiple tumours in one directory is probably OK
mkdir -p $run_dir/amber
amber_vcf=$run_dir/amber/$tumour_sample.amber.baf.vcf.gz
if [[ ! -s ${amber_vcf} ]] ; then
	write_status "Running AMBER"
	if [[ "$tumour_only" == "false" ]] ; then
		amber_args="-reference $normal_sample -reference_bam $normal_bam $amber_args"
	elif [[ "$tumour_only" == "true" ]] ; then
		amber_args="-tumor_only $amber_args"
		write_status "Running AMBER in tumour only mode."
	fi
	java -Xmx10G $jvm_args \
		-jar $amber_jar \
		-threads $threads \
		-tumor $tumour_sample \
		-tumor_bam $tumour_bam \
		-loci $bafsnps \
		-ref_genome $ref_genome \
		-validation_stringency $validation_stringency \
		-output_dir $run_dir/amber \
		$amber_args 2>&1 | tee $log_prefix.amber.log
	if [[ ! -s ${amber_vcf} ]] ; then
		write_status "Error running AMBER - aborting"
		exit 1
	fi
else
	write_status "Skipping AMBER - ${amber_vcf} exists"
fi

## cobalt output for multiple tumours in one directory is probably not OK
## Each run will output files about the same normal to the same directory -> conflict?
##mkdir -p $run_dir/cobalt
mkdir -p $run_dir/cobalt/$tumour_sample
cobalt_file=$run_dir/cobalt/$tumour_sample/$tumour_sample.cobalt.ratio.pcf
if [[ ! -s ${cobalt_file} ]] ; then
	write_status "Running COBALT"
	if [[ "$tumour_only" == "false" ]] ; then
		cobalt_args="-reference $normal_sample -reference_bam $normal_bam $cobalt_args"
	elif [[ "$tumour_only" == "true" ]] ; then
		cobalt_args="-tumor_only -tumor_only_diploid_bed ${tumor_only_diploid_bed} ${cobalt_args}"
		write_status "Running COBALT in tumour only mode."
	fi
	java -Xmx10G $jvm_args \
		-cp ${cobalt_jar} com.hartwig.hmftools.cobalt.CountBamLinesApplication \
		-threads ${threads} \
		-tumor ${tumour_sample} \
		-tumor_bam ${tumour_bam} \
		-ref_genome $ref_genome \
		-output_dir ${run_dir}/cobalt/$tumour_sample \
		-gc_profile ${gcprofile} \
		-threads ${threads} \
		$cobalt_args \
		2>&1 | tee $log_prefix.cobalt.log
	if [[ ! -s ${cobalt_file} ]] ; then
		write_status "Error running COBALT - aborting"
		exit 1
	fi
else
	write_status "Skipping COBALT - ${cobalt_file} exists"
fi

## Output for multiple tumours could probably co-exists in one directory
## But may be messy
##mkdir -p $run_dir/purple
mkdir -p $run_dir/purple/$tumour_sample

# circos requires /home/$LOGNAME to exist
if [[ -z "${LOGNAME:-}" ]] ; then
	export LOGNAME=$(whoami)
	mkdir -p /home/$LOGNAME
fi

## modified parameters on 20211215 for final purple 3.2
purple_vcf=$run_dir/purple/$tumour_sample/$tumour_sample.purple.sv.vcf.gz
if [[ ! -s ${purple_vcf} ]] ; then
	write_status "Running PURPLE"
	if [[ -s "$snvvcf" ]] ; then
		purple_args="-somatic_vcf $snvvcf $purple_args"
	fi

	if [[ "${tumour_only}" == "false" ]] ; then
		purple_args="-reference $normal_sample $purple_args"
	elif [[ "${tumour_only}" == "true" ]] ; then
		purple_args="-tumor_only $purple_args"
		write_status "Running PURPLE in tumour only mode."
	fi
#	java -Dorg.jooq.no-logo=true -Xmx10G $jvm_args \
	java -Dorg.jooq.no-logo=true -Xmx${jvmheap} $jvm_args \
		-jar ${purple_jar} \
		-output_dir $run_dir/purple/$tumour_sample \
		-tumor $tumour_sample \
		-amber $run_dir/amber \
		-cobalt $run_dir/cobalt/$tumour_sample \
		-gc_profile $gcprofile \
		-ref_genome $ref_genome \
		-ref_genome_version ${ref_genome_version} \
		-structural_vcf ${gripss_somatic_filtered_vcf} \
		-sv_recovery_vcf ${gripss_somatic_vcf} \
		-driver_catalog \
		-somatic_hotspots ${known_somatic_hotspots_vcf} \
		-germline_hotspots ${known_germline_hotspots_vcf} \
		-driver_gene_panel ${driver_gene_panel} \
		-circos circos \
		-threads ${threads} \
		$purple_args \
		2>&1 | tee $log_prefix.purple.log
	if [[ ! -s ${purple_vcf} ]] ; then
		write_status "Error running PURPLE - aborting"
		exit 1
	fi
else
	write_status "Skipping PURPLE - ${purple_vcf} exists"
fi

mkdir -p $run_dir/linx/$tumour_sample

# removed ref_genome arguments for 1.14
#	-ref_genome ${ref_genome} \

## modified parameters on 20211215 for final linx 1.17
if [[ "${tumour_only}" == "false" ]] ; then
	linx_vis_sv_data_file="$run_dir/linx/${tumour_sample}/${tumour_sample}.linx.vis_sv_data.tsv"
	if [[ ! -s ${linx_vis_sv_data_file} ]] ; then
		write_status "Running LINX"
#		java -Xmx8G -Xms4G -jar ${linx_jar} \
		java -Xmx${jvmheap} -Xms4G -jar ${linx_jar} \
			-sample ${tumour_sample} \
			-sv_vcf ${purple_vcf} \
			-purple_dir $run_dir/purple/$tumour_sample \
			-output_dir $run_dir/linx/$tumour_sample \
			-fragile_site_file ${fragile_sites} \
			-line_element_file ${line_elements} \
			-ensembl_data_dir ${ensembl_data_dir} \
			-check_fusions \
			-known_fusion_file ${known_fusion_csv} \
			-check_drivers \
			-driver_gene_panel ${driver_gene_panel} \
			-ref_genome_version ${ref_genome_version} \
			-write_vis_data \
			-log_verbose \
			$linx_args \
			2>&1 | tee $log_prefix.linx.log
	else
		write_status "Skipping LINX. ${tumour_sample}.linx.vis_sv_data.tsv exists"
	fi

	if [[ -s $run_dir/linx/$tumour_sample/plot ]] ; then
		write_status "Removing existing Linx plot output"
		rm $run_dir/linx/$tumour_sample/plot/*.png
		rm $run_dir/linx/$tumour_sample/circos/*.circos
		rm $run_dir/linx/$tumour_sample/circos/*.conf
	fi
	
	write_status "Generating LINX visualisations"
	java -cp ${linx_jar} com.hartwig.hmftools.linx.visualiser.SvVisualiser \
		-sample ${tumour_sample} \
		-ensembl_data_dir ${ensembl_data_dir} \
		-plot_out $run_dir/linx/$tumour_sample/plot/ \
		-data_out $run_dir/linx/$tumour_sample/circos/ \
		-vis_file_dir $run_dir/linx/$tumour_sample \
		-circos circos \
		-threads $threads \
		2>&1 | tee $log_prefix.linx_plot.log
else
	write_status "Skipping LINX. ${tumour_sample} is using tumor only mode"
fi		
trap - EXIT
exit 0 # success!
