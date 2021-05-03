#!/bin/bash

# Wrapper script for PRSice2: https://www.prsice.info

${RSCRIPT} ${PRSICE2_R} --prsice ${PRSICE2_SH} \
--dir ${PRSDIR} \
--bar-levels ${PRSICE_BARLEVELS} \
--base ${BASEDATA} \
--target ${VALIDATIONDATA}/${VALIDATIONPREFIX}\# \
--thread ${PRSICE_THREADS} \
--${BF_STAT,,} \
--binary-target ${BF_TARGET_TYPE} \
--type bgen \
--snp ${BF_ID_COL} \
--chr ${BF_CHR_COL} \
--bp ${BF_POS_COL} \
--A1 ${BF_EFFECT_COL} \
--A2 ${BF_NON_EFFECT_COL} \
--stat ${BF_STAT_COL} \
--pvalue ${BF_PVALUE_COL} \
--pheno-file ${PHENOTYPEFILE} \
--pheno-col ${DUMMY_PHENOTYPE} \
${PRSICE_PLOTTING} \
${PRSICE_SETTINGS} \
--out ${PRSICE_OUTPUTNAME}

# --keep /hpc/dhl_ec/aligterink/ProjectFiles/sample_data/inclusion.txt \

# PRSice_v233.R --prsice $(command -v prsice_v233) \
# --dir ${PRSICEDIR} \
# --seed ${PRSICESEED} \
# --bar-levels ${PRSICEBARLEVELS} \
# --base ${BASEDATA_disc} \
# --target ${TARGETDATA} \
# --thread ${PRSICETHREADS} \
# ${STATTYPE} \
# --binary-target ${TARGETTYPE} \
# --type ${DATATYPE} \
# --snp ${SNPID} \
# --chr ${CHRID} \
# --bp ${BPID} \
# --A1 ${A1ID} \
# --A2 ${A2ID} \
# --stat ${STATID} \
# --pvalue ${PVALUEID} \
# --cov-file ${COVARIATESFILE} \
# --cov-col ${COVARIATES} \
# --cov-factor ${COVFACTOR} \
# --pheno-file ${PHENOTYPEFILE} \
# --pheno-col ${PHENOTYPE} \
# --keep ${EXCLUSION} \
# ${PRSICEPLOTTING} \
# ${PRSICESETTINGS} \
# ${PRSICEOUTPUTNAME}