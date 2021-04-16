#!/bin/bash
#
# Wrapper script for the --score function of PLINK
#
#   S. Purcell et al., “PLINK: A tool set for whole-genome association and population-based 
#   linkage analyses,” American Journal of Human Genetics, vol. 81, no. 3, 2007, doi: 10.1086/519795.
#
# For additional documentation on PLINK --score see 
# https://www.cog-genomics.org/plink/2.0/score

# Compute the PRS score using PLINK --score

for file in ${VALIDATIONDATA}/${VALIDATIONPREFIX}*.bgen; do
    ${PLINK} \
    --bgen $file ${PLINK_REF_POS} \
    --sample ${SAMPLE_FILE} \
    --out ${PRSDIR}/plink2_$(basename $file .bgen) \
    --score ${WEIGHTS_FILE} \
    ${PLINK_VARIANT_ID_COL} \
    ${PLINK_ALLELE_COL} \
    ${PLINK_SCORE_COL} \
    ${PLINK_HEADER} \
    ${PLINK_SETTINGS}
done
