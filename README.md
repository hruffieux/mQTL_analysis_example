# Replication of an mQTL analysis using the **locus** method on simulated data 

## Overview 

A molecular quantitative trait locus analysis is performed using our 
variational inference procedure for combined predictor and outcome 
selection implemented in the R package **locus**. This example reproduces 
the application given in H. Ruffieux, A. C. Davison, J. Hager, I. Irincheeva, 
Efficient inference for genetic association studies with multiple outcomes, 
*Biostatistics*, 2017.

## Data

The SNP and metabolic expressions level datasets cannot be provided for privacy 
reason. However, these are used to simulate data that best emulate real 
conditions. The simulated dataset can be downloaded 
[here](https://dx.doi.org/10.6084/m9.figshare.4509755.v1); it cannot be 
directly generated from the scripts (though the code is provided) as it is based 
on the sample minor allele frequencies and correlation structure of the 
confidential SNP data.

## Algorithm

The package **locus** used for the analysis may be installed with the 
`devtools` command 
```R
devtools::install_github("hruffieux/locus", ref = "v0.5.0")
```
where `ref = v0.5.0` indicates the git tag corresponding to the package 
version we used.

## Workflow

The scripts should be executed in the following order:

1. `replicate_real_data_pb.R` (analysis using `locus`) and 
   `replicate_real_data_pb_varbvs.R` (analysis using the single-trait 
   variational method "varbvs" by Carbonetto and Stephens, 2012, 
   Bayesian Analysis 7);

2. `replicate_perm.R` (analysis on permuted data with `locus`) and 
   `replicate_perm_varbvs.R` (idem but with `varbvs`); and

3. `replicate_FDR_estimation.R` (FDR-based comparison of `locus` and 
   `varbvs`) and `manhattan.R` (manhattan plots for inferences with 
   `locus` and `varbvs`).

## Issues

To report an issue, please use the 
[issue tracker](https://github.com/hruffieux/mQTL_analysis_example/issues) 
at github.com.