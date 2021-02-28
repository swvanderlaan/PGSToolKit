#!/bin/bash

# This is a wrapper script for PRS-CS
#
#   T. Ge, C. Y. Chen, Y. Ni, Y. C. A. Feng, and J. W. Smoller, “Polygenic prediction 
#   via Bayesian regression and continuous shrinkage priors,” Nature Communications, 
#   vol. 10, no. 1, 2019, doi: 10.1038/s41467-019-09718-5.
#
# For additional documentation on PRS-CS see https://github.com/getian107/PRScs

# Pre-processing as GWAS statistics should be in the following format: SNP - A1 - A2 - BETA/OR - P


# Write the posterior SNP effect size estimates for each chromosome to the specified directory
${PYTHONPATH} ${PRSCS} \
--ref_dir = path_to_ref/ldblk_1kg_eur \
--bim_prefix = path_to_bim/test \
--sst_file = path_to_sumstats/sumstats.txt \
--n_gwas = 200000 \
--chrom = 22 \
--phi = 1e-2 \
--out_dir = path_to_output/eur

# Compute the PRS score for each individual in the target sample using PLINK ? 
















