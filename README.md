# tumor-grappler
Tumor-GRaPPLer is a workflow for integrating whole genome sequencing data analysis results, oriented around processing input/output from bioinformatics tools developed by the Hartwig Medical Foundation (HMF: https://github.com/hartwigmedical/hmftools) and the Papenfuss Lab (https://github.com/PapenfussLab).  In particular, Tumor-GRaPPLer is centered around processing data for GRIDSS, PURPLE, and LINX.

# Background
- GRIDSS (Genome rearrangement IDentification Software Suite)
- PURPLE (PURity PLoidy Estimator)
- LINX (Tool for annotation, interpretation, and visualizing of SVs)

# WGS and RNA-seq analysis workflows


![WGS analysks workflow](/assets/images/cancer_genome_analysis.png)


![RNA-seq analysis workflow](/assets/images/tumor_transcript_analysis.png | width=50)

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
- Somatic Alterations in Genome (SAGE) 2.8
- GRIDSS 2.12.2 (Breakend/breakpoint detection from WGS data)
- GRIDSS Post Somatic Software (GRIPSS) 1.11/2.0
## Tumor genome variant integration
- PURPLE 3.2.2
- LINX 1.17
## Mutational signature analysis
- SigProfiler

## Copy number segment analysis

# Installation

# Running the workflow
