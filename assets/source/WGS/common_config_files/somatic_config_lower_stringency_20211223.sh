#!/bin/bash

filter_hotspot_min_tumor_qual=70
filter_panel_min_tumor_qual=70
filter_high_confidence_min_tumor_qual=70
filter_low_confidence_min_tumor_qual=70

# min_tumor_vaf
hotspot_min_tumor_vaf=0.025
panel_min_tumor_vaf=0.025
high_confidence_min_tumor_vaf=0.025
low_confidence_min_tumor_vaf=0.025

# min_germline_depth
hotspot_min_germline_depth=10
panel_min_germline_depth=10
high_confidence_min_germline_depth=10
low_confidence_min_germline_depth=10

# min_germline_depth_allosome
hotspot_min_germline_depth_allosome=6
panel_min_germline_depth_allosome=6
high_confidence_min_germline_depth_allosome=6
low_confidence_min_germline_depth_allosome=6

# max_germline_vaf
hotspot_max_germline_vaf=0.0
panel_max_germline_vaf=0.0
high_confidence_max_germline_vaf=0.0
low_confidence_max_germline_vaf=0.0

# max_germline_rel_raw_base_qual
hotspot_max_germline_rel_raw_base_qual=0.04
panel_max_germline_rel_raw_base_qual=0.04
high_confidence_max_germline_rel_raw_base_qual=0.04
low_confidence_max_germline_rel_raw_base_qual=0.04

# remove common MAF SNPs and those near fixation
filter_max_maf=0.01
filter_max_af=0.99