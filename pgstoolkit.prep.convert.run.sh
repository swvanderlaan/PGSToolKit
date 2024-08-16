#!/bin/bash
# Array job: convert VCF to bgen format
CHR=${SLURM_ARRAY_TASK_ID}  # The chromosome number automatically assigned by SLURM
PLINK=$1
BGENIX=$2
STUDYDIR=$3
# Set variables - for manual per chromosome individual job version
# CHR=$1  # The chromosome number
# PLINK=$2
# BGENIX=$3
# STUDYDIR=$4
# Converting
echo "> converting chromosome ${CHR}..."
echo "...8-bits version and indexing"
if [ "$CHR" -eq 23 ]; then
    $PLINK --vcf ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.chrX.vcf.gz --export bgen-1.2 bits=8 --update-sex ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.psam --out ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.8bit.chr${CHR}
    $BGENIX -index -g ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.8bit.chr${CHR}.bgen
else
    $PLINK --vcf ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.chr${CHR}.vcf.gz --export bgen-1.2 bits=8 --update-sex ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.psam --out ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.8bit.chr${CHR}
    $BGENIX -index -g ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.8bit.chr${CHR}.bgen
fi
 
echo ""
echo "...16-bits version (default) and indexing"
if [ "$CHR" -eq 23 ]; then
    $PLINK --vcf ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.chrX.vcf.gz --export bgen-1.2 --update-sex ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.psam --out ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.16bit.chr${CHR}
    $BGENIX -index -g ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.16bit.chr${CHR}.bgen
else
    $PLINK --vcf ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.chr${CHR}.vcf.gz --export bgen-1.2 --update-sex ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.psam --out ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.16bit.chr${CHR}
    $BGENIX -index -g ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.16bit.chr${CHR}.bgen
fi

