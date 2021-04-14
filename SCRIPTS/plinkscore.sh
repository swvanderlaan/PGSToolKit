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
    echo "${file}"
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

    # ${PLINK_SUM} \
    # ${PLINK_IMPUTATION} \
    # ${PLINK_INCLUDE_CNT} \
    # ${PLINK_DOSAGE}

# # Command for PLINK 2.0
# ${PLINK} \
# --bgen /hpc/dhl_ec/data/_ae_originals/AEGS_COMBINED_EAGLE2_1000Gp3v5HRCr11/aegs.qc.1kgp3hrcr11.idfix.chr10.bgen ref-first \
# --sample /hpc/dhl_ec/aligterink/ProjectFiles/sample_data/20210309.LOOKUP.AEGS123.MFfix.lim.sample \
# --score ${PRSDIR}/weights.txt \
# 1 \
# 2 \
# 3 \
# header \
# --out ${PRSDIR}/plink2

