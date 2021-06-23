#!/bin/bash

# Wrapper script for PRSice2
#
#   S. W. Choi and P. F. O’Reilly, “PRSice-2: Polygenic Risk Score software for biobank-scale data,” 
#   GigaScience, vol. 8, no. 7, 2019, doi: 10.1093/gigascience/giz082.
#
# For additional info on PRSice2 see https://www.prsice.info

# Set the optional LD reference dataset
if [[ ${LDDATA} != "" ]]; then
    LDDATA="--ld "${LDDATA}\#
fi

# Set the target type par
if [[ ${PRSICE_PHENOTYPE_BINARY} == "TRUE" ]]; then
    PRSICE_PHENOTYPE_BINARY="T"
elif [[ ${PRSICE_PHENOTYPE_BINARY} == "FALSE" ]]; then
	PRSICE_PHENOTYPE_BINARY="F"
fi

# Set the extract file par
if [[ ${PRSICE_EXTRACT} != "" ]]; then
    PRSICE_EXTRACT="--extract "${PRSICE_EXTRACT}
fi

# Set the exclude file par
if [[ ${PRSICE_EXCLUDE} != "" ]]; then
    PRSICE_EXCLUDE="--exclude "${PRSICE_EXCLUDE}
fi

# Run PRSice
${RSCRIPT} ${PRSICE2_R} --prsice ${PRSICE2_SH} \
--dir ${PRSDIR} \
--base ${BASEDATA} \
--target ${TARGETDATA}\# \
--thread ${PRSICE_THREADS} \
--${BF_STAT,,} \
--binary-target ${PRSICE_PHENOTYPE_BINARY} \
--type bgen \
--perm ${PRSICE_PERM} \
--clump-kb ${PRSICE_CLUMP_KB} \
--clump-r2 ${PRSICE_CLUMP_R2} \
--clump-p ${PRSICE_CLUMP_P} \
${LDDATA} \
--snp ${BF_ID_COL} \
--chr ${BF_CHR_COL} \
--bp ${BF_POS_COL} \
--A1 ${BF_EFFECT_COL} \
--A2 ${BF_NON_EFFECT_COL} \
--stat ${BF_STAT_COL} \
--pvalue ${BF_PVALUE_COL} \
--pheno-file ${PHENOTYPEFILE} \
--pheno-col ${PRSICE_PHENOTYPE} \
${PRSICE_EXTRACT} \
${PRSICE_EXCLUDE} \
${PRSICE_SETTINGS} \
--out ${PRSICE_OUTPUTNAME}
