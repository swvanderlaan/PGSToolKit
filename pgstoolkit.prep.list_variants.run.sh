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

# Array job: list variants
CHR=${SLURM_ARRAY_TASK_ID}  # The chromosome number automatically assigned by SLURM
STUDYDIR=$1
# Set variables - for manual per chromosome individual job version
# CHR=$1  # The chromosome number 
# STUDYDIR=$2
# List the variants
bgenix -g ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.8bit.chr${CHR}.bgen \
-i ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.8bit.chr${CHR}.bgen.bgi \
-list | grep -v "#" > ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.8bit.chr${CHR}.variantlist.txt; 

# Deactivate the conda environment
conda deactivate