#!/bin/bash

RUN_DIR="${RUN_DIR}"
. ${RUN_DIR}/config_files/common_config.sh

GPL_PREFIX=$(basename "$PWD")

## Note that file needs to end with a blank line - newline after last row
TUMOR_NORMAL_PAIR_INFO_FILE="${RUN_DIR}/config_files/tumor_normal_pair_info_GPL.csv"

PURPLE_GRIDSS_LINX_RESULT_DIR="${RUN_DIR}/result/${GPL_PREFIX}"

MOVED_RUN=false

FIG_OUTPUT_DIR="${RUN_DIR}/figures/PURPLE_GRIDSS_LINX/${GPL_PREFIX}"
FIG_TMP_COMBINED_DIR="${FIG_OUTPUT_DIR}/tmp_combined"
FIG_TMP_COPYNUMBER_DIR="${FIG_OUTPUT_DIR}/tmp_copynumber"
FIG_TMP_INPUT_CIRCOS_DIR="${FIG_OUTPUT_DIR}/tmp_input"
FIG_TMP_MAP_DIR="${FIG_OUTPUT_DIR}/tmp_map"
FIG_TMP_PURITY_RANGE_DIR="${FIG_OUTPUT_DIR}/tmp_purity_range"
FIG_TMP_SEGMENT_DIR="${FIG_OUTPUT_DIR}/tmp_segment"
FIG_TMP_SOMATIC_CLONALITY_DIR="${FIG_OUTPUT_DIR}/tmp_somatic_clonality"
FIG_TMP_SOMATIC_DIR="${FIG_OUTPUT_DIR}/tmp_somatic"
FIG_TMP_SOMATIC_RAINFALL_DIR="${FIG_OUTPUT_DIR}/tmp_somatic_rainfall"

mkdir -p ${FIG_OUTPUT_DIR}
mkdir -p ${FIG_TMP_COMBINED_DIR}
mkdir -p ${FIG_TMP_COPYNUMBER_DIR}
mkdir -p ${FIG_TMP_INPUT_CIRCOS_DIR}
mkdir -p ${FIG_TMP_MAP_DIR}
mkdir -p ${FIG_TMP_PURITY_RANGE_DIR}
mkdir -p ${FIG_TMP_SEGMENT_DIR}
mkdir -p ${FIG_TMP_SOMATIC_CLONALITY_DIR}
mkdir -p ${FIG_TMP_SOMATIC_DIR}
mkdir -p ${FIG_TMP_SOMATIC_RAINFALL_DIR}

COMBINED_CIRCOS_PNG_FILES=""
COPYNUMBER_PNG_FILES=""
INPUT_CIRCOS_PNG_FILES=""
MAP_PNG_FILES=""
PURITY_RANGE_PNG_FILES=""
SEGMENT_PNG_FILES=""
SOMATIC_CLONALITY_PNG_FILES=""
SOMATIC_PNG_FILES=""
SOMATIC_RAINFALL_PNG_FILES=""
		
while IFS="," read -r line; do
	SUBJECT_ID=$(echo ${line} | awk -F "," '{print $1}')
	TUMOR_SAMPLE_ID=$(echo ${line} | awk -F "," '{print $4}')

	if [[ "${SUBJECT_ID}" != "subject.id" ]] ; then
		FIG_FROM_DIR="${PURPLE_GRIDSS_LINX_RESULT_DIR}/${SUBJECT_ID}/purple/${TUMOR_SAMPLE_ID}/plot"


		COMBINED_CIRCOS_PNG_FROM="${FIG_FROM_DIR}/${TUMOR_SAMPLE_ID}.circos.png"
		COMBINED_CIRCOS_PNG_TO="${FIG_TMP_COMBINED_DIR}/${TUMOR_SAMPLE_ID}.png"
		cp ${COMBINED_CIRCOS_PNG_FROM} ${COMBINED_CIRCOS_PNG_TO}
		COMBINED_CIRCOS_PNG_FILES="${COMBINED_CIRCOS_PNG_FILES} ${COMBINED_CIRCOS_PNG_TO}"

		COPYNUMBER_PNG_FROM="${FIG_FROM_DIR}/${TUMOR_SAMPLE_ID}.copynumber.png"
		COPYNUMBER_PNG_TO="${FIG_TMP_COPYNUMBER_DIR}/${TUMOR_SAMPLE_ID}.png"
		cp ${COPYNUMBER_PNG_FROM} ${COPYNUMBER_PNG_TO}
		COPYNUMBER_PNG_FILES="${COPYNUMBER_PNG_FILES} ${COPYNUMBER_PNG_TO}"
		
		INPUT_CIRCOS_PNG_FROM="${FIG_FROM_DIR}/${TUMOR_SAMPLE_ID}.input.png"
		INPUT_CIRCOS_PNG_TO="${FIG_TMP_INPUT_CIRCOS_DIR}/${TUMOR_SAMPLE_ID}.png"
		cp ${INPUT_CIRCOS_PNG_FROM} ${INPUT_CIRCOS_PNG_TO}
		INPUT_CIRCOS_PNG_FILES="${INPUT_CIRCOS_PNG_FILES} ${INPUT_CIRCOS_PNG_TO}"
		
		MAP_PNG_FROM="${FIG_FROM_DIR}/${TUMOR_SAMPLE_ID}.map.png"
		MAP_PNG_TO="${FIG_TMP_MAP_DIR}/${TUMOR_SAMPLE_ID}.png"
		cp ${MAP_PNG_FROM} ${MAP_PNG_TO}
		MAP_PNG_FILES="${MAP_PNG_FILES} ${MAP_PNG_TO}"
			
		PURITY_RANGE_PNG_FROM="${FIG_FROM_DIR}/${TUMOR_SAMPLE_ID}.purity.range.png"
		PURITY_RANGE_PNG_TO="${FIG_TMP_PURITY_RANGE_DIR}/${TUMOR_SAMPLE_ID}.png"
		cp ${PURITY_RANGE_PNG_FROM} ${PURITY_RANGE_PNG_TO}
		PURITY_RANGE_PNG_FILES="${PURITY_RANGE_PNG_FILES} ${PURITY_RANGE_PNG_TO}"

		SEGMENT_PNG_FROM="${FIG_FROM_DIR}/${TUMOR_SAMPLE_ID}.segment.png"
		SEGMENT_PNG_TO="${FIG_TMP_SEGMENT_DIR}/${TUMOR_SAMPLE_ID}.png"
		cp ${SEGMENT_PNG_FROM} ${SEGMENT_PNG_TO}
		SEGMENT_PNG_FILES="${SEGMENT_PNG_FILES} ${SEGMENT_PNG_TO}"

		SOMATIC_CLONALITY_PNG_FROM="${FIG_FROM_DIR}/${TUMOR_SAMPLE_ID}.somatic.clonality.png"
		SOMATIC_CLONALITY_PNG_TO="${FIG_TMP_SOMATIC_CLONALITY_DIR}/${TUMOR_SAMPLE_ID}.png"
		cp ${SOMATIC_CLONALITY_PNG_FROM} ${SOMATIC_CLONALITY_PNG_TO}
		SOMATIC_CLONALITY_PNG_FILES="${SOMATIC_CLONALITY_PNG_FILES} ${SOMATIC_CLONALITY_PNG_TO}"
		
		SOMATIC_PNG_FROM="${FIG_FROM_DIR}/${TUMOR_SAMPLE_ID}.somatic.png"
		SOMATIC_PNG_TO="${FIG_TMP_SOMATIC_DIR}/${TUMOR_SAMPLE_ID}.png"
		cp ${SOMATIC_PNG_FROM} ${SOMATIC_PNG_TO}
		SOMATIC_PNG_FILES="${SOMATIC_PNG_FILES} ${SOMATIC_PNG_TO}"
		
		SOMATIC_RAINFALL_PNG_FROM="${FIG_FROM_DIR}/${TUMOR_SAMPLE_ID}.somatic.rainfall.png"
		SOMATIC_RAINFALL_PNG_TO="${FIG_TMP_SOMATIC_RAINFALL_DIR}/${TUMOR_SAMPLE_ID}.png"
		cp ${SOMATIC_RAINFALL_PNG_FROM} ${SOMATIC_RAINFALL_PNG_TO}
		SOMATIC_RAINFALL_PNG_FILES="${SOMATIC_RAINFALL_PNG_FILES} ${SOMATIC_RAINFALL_PNG_TO}"

		echo "Copied figures for ${TUMOR_SAMPLE_ID}"
	fi
done < ${TUMOR_NORMAL_PAIR_INFO_FILE}

echo "Merging combined circos png files into pdf."
convert ${COMBINED_CIRCOS_PNG_FILES} -gravity North -font "Liberation-Sans" -pointsize 64 -annotate 0 '%t' ${FIG_OUTPUT_DIR}/${GPL_PREFIX}.combined_circos.pdf

echo "Merging copynumber png files into pdf."
convert ${COPYNUMBER_PNG_FILES} -gravity  NorthEast -font "Liberation-Sans" -pointsize 32 -annotate 0 '%t' ${FIG_OUTPUT_DIR}/${GPL_PREFIX}.copynumber.pdf

echo "Merging circos input png files into pdf."
convert ${INPUT_CIRCOS_PNG_FILES} -gravity North -font "Liberation-Sans" -pointsize 64 -annotate 0 '%t' ${FIG_OUTPUT_DIR}/${GPL_PREFIX}.circos_input.pdf
		
echo "Merging map png files into pdf."
convert ${MAP_PNG_FILES} -gravity  NorthEast -font "Liberation-Sans" -pointsize 32 -annotate 0 '%t' ${FIG_OUTPUT_DIR}/${GPL_PREFIX}.map.pdf

echo "Merging purity png files into pdf."
convert ${PURITY_RANGE_PNG_FILES} -gravity  NorthEast -font "Liberation-Sans" -pointsize 32 -annotate 0 '%t' ${FIG_OUTPUT_DIR}/${GPL_PREFIX}.purity.range.pdf

echo "Merging segment png files into pdf."
convert ${SEGMENT_PNG_FILES} -gravity  NorthEast -font "Liberation-Sans" -pointsize 32 -annotate 0 '%t' ${FIG_OUTPUT_DIR}/${GPL_PREFIX}.segment.pdf

echo "Merging somatic clonality png files into pdf."
convert ${SOMATIC_CLONALITY_PNG_FILES} -gravity  NorthEast -font "Liberation-Sans" -pointsize 32 -annotate 0 '%t' ${FIG_OUTPUT_DIR}/${GPL_PREFIX}.somatic.clonality.pdf

echo "Merging somatic png files into pdf."
convert ${SOMATIC_PNG_FILES} -gravity NorthEast -font "Liberation-Sans" -pointsize 32 -annotate 0 '%t' ${FIG_OUTPUT_DIR}/${GPL_PREFIX}.somatic.pdf

echo "Merging somatic rainfall png files into pdf."
convert ${SOMATIC_RAINFALL_PNG_FILES} -gravity NorthEast -font "Liberation-Sans" -pointsize 32 -annotate 0 '%t' ${FIG_OUTPUT_DIR}/${GPL_PREFIX}.somatic.rainfall.pdf

echo "Finished merging.  Now cleaning up temporary files."
rm ${FIG_TMP_COMBINED_DIR}/*
rm ${FIG_TMP_COPYNUMBER_DIR}/*
rm ${FIG_TMP_INPUT_CIRCOS_DIR}/*
rm ${FIG_TMP_MAP_DIR}/*
rm ${FIG_TMP_PURITY_RANGE_DIR}/*
rm ${FIG_TMP_SEGMENT_DIR}/*
rm ${FIG_TMP_SOMATIC_CLONALITY_DIR}/*
rm ${FIG_TMP_SOMATIC_DIR}/*
rm ${FIG_TMP_SOMATIC_RAINFALL_DIR}/*

echo "Finished"
