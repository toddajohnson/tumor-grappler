#!/bin/bash

## added SAGE arguments
# min_tumor_qual
hotspot_min_tumor_qual=70
panel_min_tumor_qual=70
high_confidence_min_tumor_qual=70
low_confidence_min_tumor_qual=70

filter_hotspot_min_tumor_qual=70
filter_panel_min_tumor_qual=100
filter_high_confidence_min_tumor_qual=100
filter_low_confidence_min_tumor_qual=100

# min_tumor_vaf
hotspot_min_tumor_vaf=0.025
panel_min_tumor_vaf=0.025
high_confidence_min_tumor_vaf=0.025
low_confidence_min_tumor_vaf=0.025

# remove common MAF SNPs and those near fixation
filter_max_maf=0.01
filter_max_af=0.99
