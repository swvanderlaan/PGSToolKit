#!/bin/bash
CHR=$1
$PLINK=$2
$BGENIX=$3
$STUDYDIR=$4
echo "> converting chromosome ${CHR}..."
echo "...8-bits version and indexing"
$PLINK --vcf ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.chr${CHR}.vcf.gz --export bgen-1.2 bits=8 --out ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.8bit.chr${CHR}
$BGENIX -index -g ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.8bit.chr${CHR}.bgen 
echo ""
echo "...16-bits version (default) and indexing"
$PLINK --vcf ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.chr${CHR}.vcf.gz --export bgen-1.2 --out ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.16bit.chr${CHR}
$BGENIX -index -g ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.16bit.chr${CHR}.bgen
