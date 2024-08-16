#!/bin/bash
# Load the required conda environment and check if conda activate was successful
echo "Loading required mamba environment containing the pgstoolkit installation..."
eval "$(conda shell.bash hook)"
conda activate pgstoolkit

if [ $? -ne 0 ]; then
    echo "Error: Failed to activate pgstoolkit environment." >&2
    exit 1
fi
echo "> Checking existence of relevant apps..."
bgenix -help
bcftools --version
vcftools --version
samtools --version

# Array job: convert VCF to bgen format
CHR=${SLURM_ARRAY_TASK_ID}  # The chromosome number automatically assigned by SLURM
PLINK=$1
STUDYDIR=$2
# Set variables - for manual per chromosome individual job version
# CHR=$1  # The chromosome number
# PLINK=$2
# STUDYDIR=$3
# Converting
echo "> converting chromosome ${CHR}..."
echo "...8-bits version and indexing"
if [ "$CHR" -eq 23 ]; then
    $PLINK --vcf ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.chrX.vcf.gz --export bgen-1.2 bits=8 --update-sex ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.psam --out ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.8bit.chr${CHR}
    bgenix -index -g ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.8bit.chr${CHR}.bgen
else
    $PLINK --vcf ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.chr${CHR}.vcf.gz --export bgen-1.2 bits=8 --update-sex ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.psam --out ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.8bit.chr${CHR}
    bgenix -index -g ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.8bit.chr${CHR}.bgen
fi
 
echo ""
echo "...16-bits version (default) and indexing"
if [ "$CHR" -eq 23 ]; then
    $PLINK --vcf ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.chrX.vcf.gz --export bgen-1.2 --update-sex ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.psam --out ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.16bit.chr${CHR}
    bgenix -index -g ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.16bit.chr${CHR}.bgen
else
    $PLINK --vcf ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.chr${CHR}.vcf.gz --export bgen-1.2 --update-sex ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.psam --out ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.16bit.chr${CHR}
    bgenix -index -g ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.16bit.chr${CHR}.bgen
fi

# Deactivate the conda environment
conda deactivate