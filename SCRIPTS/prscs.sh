#!/bin/bash

# This is a wrapper script for PRS-CS
#
#   T. Ge, C. Y. Chen, Y. Ni, Y. C. A. Feng, and J. W. Smoller, “Polygenic prediction 
#   via Bayesian regression and continuous shrinkage priors,” Nature Communications, 
#   vol. 10, no. 1, 2019, doi: 10.1038/s41467-019-09718-5.
#
# For additional documentation on PRS-CS see https://github.com/getian107/PRScs

# Write the posterior SNP effect size estimates for each chromosome to the specified directory
${PYTHONPATH} ${PRSCS} \
--ref_dir = path_to_ref/ldblk_1kg_eur \
--bim_prefix = path_to_bim/test \
--sst_file = ${BASEDATA} \
--n_gwas = 200000 \
--chrom = 22 \
--phi = 1e-2 \
--out_dir = path_to_output/eur

# Compute the PRS score using PLINK --score
# Documentation: https://www.cog-genomics.org/plink/1.9/score
${PLINK} --score <filename> \
${PLINK_VARIANTID_COL} \
${PLINK_ALLELE_COL} \
${PLINK_SCORE_COL} \
${PLINK_HEADER} \
${PLINK_SUM} \
${PLINK_IMPUTATION} \
${PLINK_INCLUDE_CNT} \
${PLINK_DOSAGE}














