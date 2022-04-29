#!/bin/bash

# Variables passed from calling script
run_dir="${RUN_DIR}"
threads="${threads}"
jvmheap=${jvmheap}
linx_plot_info_file=${linx_plot_info_file}
linx_config_file=${linx_config_file}
plot_sub_directory=${plot_sub_directory}

. ${run_dir}/config_files/common_config.sh
. ${linx_config_file}

module use /usr/local/package/modulefiles/
module load java/11
export PATH="/home/tjohnson/.local/bin:/home/tjohnson/tools/circos-0.69-9/bin:${PATH}"

## First line is header
let "LINENUM = ${SGE_TASK_ID} + 1"

linx_plot_info=$(sed -n "${LINENUM}p" ${linx_plot_info_file})
subject_id=$(echo "${linx_plot_info}" | awk -F "\t" '{print $1}')
sample_id=$(echo "${linx_plot_info}" | awk -F "\t" '{print $2}')
chromosome=$(echo "${linx_plot_info}" | awk -F "\t" '{print $3}')
gene_list=$(echo "${linx_plot_info}" | awk -F "\t" '{print $4}')
chrom_plot_str=$(echo "${linx_plot_info}" | awk -F "\t" '{print $5}')
cluster_plot_str=$(echo "${linx_plot_info}" | awk -F "\t" '{print $6}')

echo ${subject_id}
echo ${sample_id}
echo ${chromosome}
echo ${gene_list}
echo ${chrom_plot_str}
echo ${cluster_plot_str}

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

linx_jar=$(find_jar LINX_JAR sv-linx)

if [[ "${REFERENCE_GENOME}" == "HG38" ]] || [[ "${REFERENCE_GENOME}" == "38" ]] || [[ "${REFERENCE_GENOME}" == "GRCh38" ]] ; then
	reference_genome_version=38
elif [[ "${REFERENCE_GENOME}" == "HG37" ]] || [[ "${REFERENCE_GENOME}" == "37" ]] || [[ "${REFERENCE_GENOME}" == "hg19" ]]; then
	reference_genome_version=37
fi


ref_dir="/home/tjohnson/reference/HMF/${reference_genome_version}"
ensembl_cache_dir="${ref_dir}/dbs/linx"
	
data_out=$run_dir/result/SVs/Linx/circos/${plot_sub_directory}/$sample_id
mkdir -p ${data_out}

vis_data_dir=$run_dir/result/${gpl_prefix}/${subject_id}/linx/${sample_id}


if [[ "${cluster_plot_str}" != "none" ]] ; then
	echo "Constructing cluster Linx plot for ${sample_id}: cluster ${cluster_plot_str}"

	plot_out=$run_dir/result/SVs/Linx/plots/${plot_sub_directory}/clusters
	mkdir -p ${plot_out}
	
	java -cp ${linx_jar} com.hartwig.hmftools.linx.visualiser.SvVisualiser \
		-sample ${sample_id} \
		-ensembl_data_dir ${ensembl_cache_dir} \
		-ref_genome_version ${reference_genome_version} \
		-clusterId ${cluster_plot_str} \
		-plot_out ${plot_out} \
		-data_out ${data_out} \
		-vis_file_dir ${vis_data_dir} \
		-circos circos \
		-gene ${gene_list} \
		-fusion_legend_height_per_row ${fusion_legend_height_per_row} \
		-segment_relative_size ${segment_relative_size} \
		-outer_radius ${outer_radius} \
		-min_line_size ${min_line_size} -max_line_size ${max_line_size} \
		-min_label_size ${min_label_size} -max_label_size ${max_label_size} \
		-glyph_size ${glyph_size} \
		-exon_rank_radius ${exon_rank_radius}
fi

if [[ "${chrom_plot_str}" != "none" ]] ; then
	echo "Constructing chromosome Linx plot for ${sample_id}: cluster ${chrom_plot_str}"
	plot_out=$run_dir/result/SVs/Linx/plots/${plot_sub_directory}/chromosomes
	mkdir -p ${plot_out}
	
	java -cp ${linx_jar} com.hartwig.hmftools.linx.visualiser.SvVisualiser \
		-sample ${sample_id} \
		-ensembl_data_dir ${ensembl_cache_dir} \
		-ref_genome_version ${reference_genome_version} \
		-chromosome ${chrom_plot_str} \
		-plot_out ${plot_out} \
		-data_out ${data_out} \
		-vis_file_dir ${vis_data_dir} \
		-circos circos \
		-gene ${gene_list} \
		-fusion_legend_height_per_row ${fusion_legend_height_per_row} \
		-segment_relative_size ${segment_relative_size} \
		-outer_radius ${outer_radius} \
		-min_line_size ${min_line_size} -max_line_size ${max_line_size} \
		-min_label_size ${min_label_size} -max_label_size ${max_label_size} \
		-glyph_size ${glyph_size} \
		-exon_rank_radius ${exon_rank_radius}
fi

echo "Generating LINX visualisations for ${sample_id}: ${region_str}"

