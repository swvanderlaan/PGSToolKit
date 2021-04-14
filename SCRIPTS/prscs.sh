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

# tmp arguments for chr21
BIM_FILE_AVAILABLE="n"
VALIDATIONPREFIX=

if [ "$BIM_FILE_AVAILABLE" == "NO" ]; then
    for file in $VALIDATIONDATA/$VALIDATIONPREFIX*.bgen; do
        ${PLINK} \
        --bgen $file ref-first \
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

# Parse the base file to PRS-CS format (and also unzip it)
PARSED_BASEDATA=${PRSDIR}/basefile_PRScs_format.txt

# ${RSCRIPTPATH} ${PRSTOOLKITSCRIPTS}/basefile_PRScs_formatter.R \
# -i ${BASEDATA} -o ${PARSED_BASEDATA} -d ${SNP_COL} \
# -r ${A1_COL} -a ${A2_COL} -z ${MEASURE} -m ${MSR_COL} -v ${PVALUE_COL}



${PYTHONPATH} ${PRSCS} \
--ref_dir=${LDDATA} \
--bim_prefix=${BIM_FILE_PATH%.*} \
--sst_file=${PARSED_BASEDATA} \
--n_gwas=${GWAS_SAMPLE_SIZE} \
--chrom=${PRSCS_CHROM} \
--out_dir=${PRSDIR} \
--a=${PRSCS_PARAM_A} \
--b=${PRSCS_PARAM_B} \
--phi=${PRSCS_PARAM_PHI} \
--n_iter=${PRSCS_MCMC_ITERATIONS} \
--n_burnin=${PRSCS_MCMC_BURNIN} \
--seed=${PRSCS_SEED}




# 
# ${PYTHONPATH} ${PRSCS} \
# --ref_dir=${LDDATA} \
# --bim_prefix=${BIM_FILE_PATH%.*} \
# --sst_file=${PARSED_BASEDATA} \
# --n_gwas=${GWAS_SAMPLE_SIZE} \
# --chrom=${PRSCS_CHROM} \
# --out_dir=${PRSDIR} \
# --a=${PRSCS_PARAM_A} \
# --b=${PRSCS_PARAM_B} \
# --phi=${PRSCS_PARAM_PHI} \
# --n_iter=${PRSCS_MCMC_ITERATIONS} \
# --n_burnin=${PRSCS_MCMC_BURNIN} \
# --seed=${PRSCS_SEED}

#--thin=${PRSCS_MCM_THINNING} \


# ${PYTHONPATH} ${PRSCS} \
# --ref_dir=path_to_ref/ldblk_1kg_eur \
# --bim_prefix=/hpc/local/CentOS7/dhl_ec/software/PRScs_04Jan2021/test_data/test \
# --sst_file=/hpc/local/CentOS7/dhl_ec/software/PRScs_04Jan2021/test_data/sumstats.txt \
# --n_gwas=200000 \
# --chrom=22 \
# --phi=1e-2 \
# --out_dir=path_to_output/eur


# Compute the PRS score using PLINK --score


# ${PYTHONPATH} ${PRSCS} --help