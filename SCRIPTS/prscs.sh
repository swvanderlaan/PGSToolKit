#!/bin/bash
#
# This is a wrapper script for PRS-CS
#
#   T. Ge, C. Y. Chen, Y. Ni, Y. C. A. Feng, and J. W. Smoller, “Polygenic prediction 
#   via Bayesian regression and continuous shrinkage priors,” Nature Communications, 
#   vol. 10, no. 1, 2019, doi: 10.1038/s41467-019-09718-5.
#
# For additional documentation on PRS-CS see https://github.com/getian107/PRScs

# If a .bim file is not available, we will generate one for each file in the target dataset 
# and then merge these files into a single .bim file
if [ "$BIM_FILE_AVAILABLE" == "NO" ]; then
    for file in $VALIDATIONDATA/$VALIDATIONPREFIX*.bgen; do
        ${PLINK} \
        --bgen $file ${REF_POS} \
        --sample ${SAMPLE_FILE} \
        --make-just-bim \
        --out ${PRSDIR}/$(basename $file .bgen)
    done

    # Create combined .bim file
    BIM_FILE_PATH=${PRSDIR}/COMBINED_VALIDATION.bim
    cat ${PRSDIR}/*.bim > $BIM_FILE_PATH    

elif [ "$BIM_FILE_AVAILABLE" == "YES" ]; then
    echo ".bim file: ${BIM_FILE_PATH}"
else
    echo "BIM_FILE_AVAILABLE argument \"$BIM_FILE_AVAILABLE\" not recognized"
    exit 1
fi

# Perform PRS-CS
${PYTHONPATH} ${PRSCS} \
--ref_dir=${LDDATA} \
--bim_prefix=${BIM_FILE_PATH%.*} \
--sst_file=${PARSED_BASEDATA} \
--n_gwas=${BF_SAMPLE_SIZE} \
--out_dir=${PRSDIR}/single_chr_file

# Combine the ouput files into a single file
cat ${PRSDIR}/single_chr_file*.txt > ${WEIGHTS_FILE}
