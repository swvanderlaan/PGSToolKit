#!/bin/bash
CHR=$1
$PLINK=$2
$STUDYDIR=$3
echo "> calculating frequencies for chromosome ${CHR}..."
$PLINK \
--bgen ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.8bit.chr${CHR}.bgen 'ref-first' \
--sample ${STUDYDATADIR}/AEGS_QC_imputation_2023/OUTPUT/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.forPGSToolKit.sample \
--freq --out ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.8bit.chr${CHR}.FREQ; 