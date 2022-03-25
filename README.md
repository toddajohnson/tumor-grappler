# tumor-grappler
Tumor-GRaPPLer is a workflow for integrating whole genome sequencing data analysis results, oriented around processing input/output from bioinformatics tools developed by the Hartwig Medical Foundation (HMF: https://github.com/hartwigmedical/hmftools) and the Papenfuss Lab (https://github.com/PapenfussLab).  In particular, Tumor-GRaPPLer is centered around processing data for GRIDSS, PURPLE, and LINX.

# Background
- GRIDSS (Genome rearrangement IDentification Software Suite)
- PURPLE (purity ploidy estimator)
- LINX ()

# Tools used by this workflow
## Sequence alignment and processing (DNA and RNA-seq)
- BWA 0.7.17
- samtools
- bedtools
- STAR 
- Isofox 1.3
- AMBER (B-allele frequency and contamination)
- COBALT (read depth coverage calculation)
## Annotation
- Annovar
- snpEff
## Mutation and structural variant (SV) calling
- SAGE 3.2
- GRIDSS 2.12.2 (Breakend/breakpoint detection from WGS data)
- GRIPSS 2.0 (GRIDSS SV filtering)
## Tumor genome variant integration
- PURPLE 3.2.2
- LINX 1.17
## Mutational signature analysis
- SigProfiler

## Copy number segment analysis

# Installation

# Running the workflow
