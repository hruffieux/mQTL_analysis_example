# mQTL_analysis_example - Replication of a real metabolite quantitative trait loci (mQTL) analysis using simulated data. 

## Overview 

The analysis uses our variational inference procedure for combined
predictor and outcome selection (Helene Ruffieux, Anthony C. Davison,
Jorg Hager, Irina Irincheeva, 2016, arXiv:1609.03400), which is
implemented in the R package `locus`.

## Data

The data can be downloaded 
[here](https://dx.doi.org/10.6084/m9.figshare.4509755.v1); they cannot be 
generated from the scripts as the sample minor allele frequencies and 
correlation structure of confidential SNP data are used to emulate 
real conditions.

## Algorithm

The package `locus` used for the analysis may be installed with the 
devtools command 
`devtools::install_github("hruffieux/locus")`

## Workflow

The scripts should be executed in the following order:

1. `replicate_real_data_pb.R` (analysis using `locus`) and 
   `replicate_real_data_pb_varbvs.R` (analysis using the single-trait 
   variational method "varbvs" by Carbonetto and Stephens, 2012, 
   Bayesian Analysis 7);

2. `replicate_perm.R` (analysis on permuted data with locus) and 
   `replicate_perm_varbvs.R` (idem but with `varbvs`);

3. `replicate_FDR_estimation.R` (FDR-based comparison of locus and varbvs)
   and `manhattan.R` (manhattan plots for inferences with `locus` and `varbvs`).

     