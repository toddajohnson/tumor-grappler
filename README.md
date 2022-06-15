# tumor-grappler
Tumor-GRaPPLer is a workflow for integrating whole genome sequencing data analysis results, oriented around processing input/output from bioinformatics tools developed by the Hartwig Medical Foundation (HMF: https://github.com/hartwigmedical/hmftools) and the Papenfuss Lab (https://github.com/PapenfussLab).  In particular, Tumor-GRaPPLer is centered around processing output from GRIDSS, PURPLE, and LINX, which analyzes small variants (germline and somatic SNVs and INDELs), structural variants, and copy number (CN) data to determine sample quality parameters and candidate driver variants.  Downstream WGS analyses include GISTIC2, CNApp, CHORD, SigProfiler, and maftools, as well as in-house R code to process and summarize the data.

Tumor-grappler is a repository for code used in a number of the studies in which I am involved, so at the moment it is a work in progress as towards documentation.

# Background
- GRIDSS (Genome rearrangement IDentification Software Suite)
- PURPLE (PURity PLoidy Estimator)
- LINX (Tool for annotation, interpretation, and visualizing of SVs)

# WGS and RNA-seq analysis workflows

<img src="https://github.com/toddajohnson/tumor-grappler/blob/main/assets/images/cancer_genome_analysis.png" width="800">
<img src="https://github.com/toddajohnson/tumor-grappler/blob/main/assets/images/tumor_transcript_analysis.png" width="500">

# Tools used by this workflow
## Sequence alignment and processing (DNA and RNA-seq)
- BWA 0.7.17
- samtools
- bedtools
- STAR 2.7.9a
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
- tumor-grappler/CNA_analysis (R code for CNA frequency summarization, Circos plot creation, chromosome-arm level CNA plot)
- GISTIC2
- CNApp
## Other analyses
- CHORD (Homologous recombination deficiency estimation)
- maftools
- 
# Installation

# Running the workflow
