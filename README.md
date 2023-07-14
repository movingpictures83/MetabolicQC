# MetabolicQC
# Language: R
# Input: TXT
# Output: DIR 
# Tested with: PluMA 1.1, R 4.0.0
# Dependencies: [1] ggplot2_3.4.0               MultiAssayExperiment_1.16.0 SummarizedExperiment_1.20.0 Biobase_2.50.0 GenomicRanges_1.42.0        GenomeInfoDb_1.26.7 IRanges_2.24.1              S4Vectors_0.28.1 BiocGenerics_0.36.1         MatrixGenerics_1.2.1 matrixStats_0.63.0

PluMA plugin to produce quality control graphs for metabolimics data.

Input is a tab-delimited text file of keyword-value pairs:

assay: File of observables and counts
metadata: Sample metadata
annotations: Chemical annotations

Output will be produced in the user-specified output directory.
