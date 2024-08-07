#!/bin/bash
# Array job: list variants
# CHR=${SLURM_ARRAY_TASK_ID}  # The chromosome number automatically assigned by SLURM
# BGENIX=$1
# STUDYDIR=$2
# Set variables
CHR=$1  # The chromosome number 
BGENIX=$2
STUDYDIR=$3
# List the variants
$BGENIX - g ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.8bit.chr${CHR}.bgen \
-i ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.8bit.chr${CHR}.bgen.bgi \
-list | grep -v "#" > ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.8bit.chr${CHR}.variantlist.txt; 
