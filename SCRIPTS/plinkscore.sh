#!/bin/bash
#
# Wrapper script for the --score function of PLINK
#
#   S. Purcell et al., “PLINK: A tool set for whole-genome association and population-based 
#   linkage analyses,” American Journal of Human Genetics, vol. 81, no. 3, 2007, doi: 10.1086/519795.
#
# For additional documentation on PLINK --score see https://www.cog-genomics.org/plink/2.0/score

# Function for finding the index of column names in the header of the weights file
find_index() {
    i=1
    for col in "${HEADER_ARRAY[@]}"; do
        [[ $col == "$COLUMN" ]] && { echo "$i"; break; }
        (( ++i ))
    done
}

# Plink can't read zipped files, if the file ends with .gz we will unzip the file
if [[ ${WEIGHTS_FILE} == *.gz ]]; then
    echo "Compressed weight file will be decompressed as PLINK can't read compressed files."
    gzip -d -c ${WEIGHTS_FILE} > ${PRSDIR}/plink2_unzipped_basefile.txt
    WEIGHTS_FILE=${PRSDIR}/plink2_unzipped_basefile.txt
fi

# PLINK_HEADER indicates whether the weights file has a header, if not then the PLINK column parameters
# are assumed to be indices instead of column names. This is required as some PRS methods (e.g. PRS-CS)
# generate weights files without headers.
if [[ ${PLINK_HEADER} == "FALSE" ]]; then
    PLINK_VARIANT_ID_COL_NR=${SNP_COL}
    PLINK_ALLELE_COL_NR=${EFFECT_COL}
    PLINK_SCORE_COL_NR=${SCORE_COL}
    HEADER_PAR=""

elif [[ ${PLINK_HEADER} == "TRUE" ]]; then

    # Get header of weights file
    HEADER_ARRAY=($(cat ${WEIGHTS_FILE} | head -1))

    # Find indexes of required columns
    COLUMN=${SNP_COL}
    PLINK_VARIANT_ID_COL_NR=$(find_index)
    COLUMN=${EFFECT_COL}
    PLINK_ALLELE_COL_NR=$(find_index)
    COLUMN=${SCORE_COL}
    PLINK_SCORE_COL_NR=$(find_index)

    HEADER_PAR="header"

else
    echo "BF_PLINK_COLS_ARE_INDEX parameter not recognized, exiting..."
    exit 1

fi

# Compute the PRS score using PLINK --score
for file in ${VALIDATIONDATA}/${VALIDATIONPREFIX}*.bgen; do
    ${PLINK} \
    --bgen $file ${REF_POS} \
    --sample ${SAMPLE_FILE} \
    --out ${PRSDIR}/plink2_$(basename $file .bgen) \
    --score ${WEIGHTS_FILE} cols=+scoresums \
    ${PLINK_VARIANT_ID_COL_NR} \
    ${PLINK_ALLELE_COL_NR} \
    ${PLINK_SCORE_COL_NR} \
    ${HEADER_PAR} \
    ${PLINK_SETTINGS}
done
