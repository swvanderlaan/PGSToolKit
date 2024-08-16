#!/bin/bash
#################################################################################################
### PARAMETERS SLURM
#SBATCH --job-name=prepPGSTK                                  		# the name of the job
#SBATCH --output=/hpc/dhl_ec/svanderlaan/projects/prepPGSTK.log 	    # the log file of this job
#SBATCH --error=/hpc/dhl_ec/svanderlaan/projects/prepPGSTK.errors	# the error file of this job
#SBATCH --time=00:30:00                                             # the amount of time the job will take: -t [min] OR -t [days-hh:mm:ss]
#SBATCH --mem=8G                                                   # the amount of memory you think the script will consume, found on: https://wiki.bioinformatics.umcutrecht.nl/bin/view/HPC/SlurmScheduler
#SBATCH --gres=tmpspace:32G                                        # the amount of temporary diskspace per node
#SBATCH --mail-user=s.w.vanderlaan-2@umcutrecht.nl                  # where should be mailed to?
#SBATCH --mail-type=FAIL                                            # when do you want to receive a mail from your job?  Valid type values are NONE, BEGIN, END, FAIL, REQUEUE
                                                                    # or ALL (equivalent to BEGIN, END, FAIL, INVALID_DEPEND, REQUEUE, and STAGE_OUT), 
                                                                    # Multiple type values may be specified in a comma separated list. 
####    Note:   You do not have to specify workdir: 
####            'Current working directory is the calling process working directory unless the --chdir argument is passed, which will override the current working directory.'
####            TODO: select the type of interpreter you'd like to use
####            TODO: Find out whether this job should dependant on other scripts (##SBATCH --depend=[state:job_id])
####
#################################################################################################
###
### Excellent resource: https://2cjenn.github.io/PRS_Pipeline/
###

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

### Set the software
SOFTWARE="/hpc/local/Rocky8/dhl_ec/software"
PGSTK="${SOFTWARE}/PGSToolKit"

PLINK="${SOFTWARE}/plink2_linux_x86_64_20240105_alpha_5_10/plink2"

### Set the project directory
PROJECTDIR="/hpc/dhl_ec/svanderlaan/projects/"
PGSDIR="/hpc/dhl_ec/svanderlaan/projects/polygenicscores"

### Set the study directory
STUDYDATADIR="/hpc/dhl_ec/data/_ae_originals/"
### b38 -- TOPMed imputed
STUDYDIR="${STUDYDATADIR}/AEGS_QC_imputation_2023/aegscombo/_topmed_r3_f10_b38"
### b37 version -- deprecated
### STUDYDIR="${STUDYDATADIR}/AEGS_COMBINED_EAGLE2_1000Gp3v5HRCr11"

### Job settings
PREP_TIME="02:00:00"
PREP_MEM="10G"
PREP_MAIL="s.w.vanderlaan-2@umcutrecht.nl"
PREP_MAILTYPE="FAIL"

echo "========================================================================"
echo "STEP 1 PGSToolKit pipeline"
echo "========================================================================"
echo ""
echo "Note that this part is only run once to convert the VCF to bgen format,
create the variant list, calculate frequencies, and create the input for PGSToolKit."
echo ""
echo "1. -- Converting VCF to bgen including indexing. This is done for 8-bit and 16-bit versions."
### b38
### This data includes the correct ID-type (chr#:BP, e.g. chr1:10711)
# convert_job_id=$(sbatch --parsable --array=1-23 --time ${PREP_TIME} --mem ${PREP_MEM} --mail-user ${PREP_MAIL} --mail-type ${PREP_MAILTYPE} ${PGSTK}/pgstoolkit.prep.convert.run.sh $PLINK $STUDYDIR)
# For-loop individual job version
# for CHR in {1..23}; do
#     echo "> submitting conversion job for chromosome ${CHR}..."
#     convert_job_id=$(sbatch --parsable --time ${PREP_TIME} --mem ${PREP_MEM} --mail-user ${PREP_MAIL} --mail-type ${PREP_MAILTYPE} ${PGSTK}/pgstoolkit.prep.convert.run.sh $CHR $PLINK $STUDYDIR)
# done

echo ""
echo "2. -- Creating variant lists."
# list_variants_job_id=$(sbatch --parsable --array=1-23 --dependency=afterok:${convert_job_id} --time ${PREP_TIME} --mem ${PREP_MEM} --mail-user ${PREP_MAIL} --mail-type ${PREP_MAILTYPE} ${PGSTK}/pgstoolkit.prep.list_variants.run.sh $STUDYDIR)
list_variants_job_id=$(sbatch --parsable --array=1-23 --time ${PREP_TIME} --mem ${PREP_MEM} --mail-user ${PREP_MAIL} --mail-type ${PREP_MAILTYPE} ${PGSTK}/pgstoolkit.prep.list_variants.run.sh $STUDYDIR)
# For-loop individual job version
# for CHR in {1..23}; do
#     echo "> submitting job to list variants for chromosome ${CHR}..."
#     list_variants_job_id=$(sbatch --parsable --time ${PREP_TIME} --mem ${PREP_MEM} --mail-user ${PREP_MAIL} --mail-type ${PREP_MAILTYPE} ${PGSTK}/pgstoolkit.prep.list_variants.run.sh $CHR $STUDYDIR)
# done

# Variantlist looks like this
# alternate_ids	rsid	chromosome	position	number_of_alleles	first_allele	alternative_alleles
# rs367896724	rs367896724	1	10177	2	A	AC
# rs540431307	rs540431307	1	10235	2	T	TA
# rs555500075	rs555500075	1	10352	2	T	TA
# rs537182016	rs537182016	1	10539	2	C	A

echo ""
echo "3. -- Concatenating variant list."
concat_variants_job_id=$(sbatch --parsable --dependency=afterok:${list_variants_job_id} --time ${PREP_TIME} --mem ${PREP_MEM} --mail-user ${PREP_MAIL} --mail-type ${PREP_MAILTYPE} ${PGSTK}/pgstoolkit.prep.concat_variants.run.sh $STUDYDIR)

# Concatenated variantlist looks like this
# variantid alt_ids rsid_aegs chromosome position number_of_alleles first_allele alt_alleles
# 1:10177 rs367896724 rs367896724 1 10177 2 A AC
# 1:10235 rs540431307 rs540431307 1 10235 2 T TA
# 1:10352 rs555500075 rs555500075 1 10352 2 T TA
# 1:10539 rs537182016 rs537182016 1 10539 2 C A
# 1:10616 rs376342519 rs376342519 1 10616 2 CCGCCGTTGCAAAGGCGCGCCG C

echo ""
echo "4. -- Calculating frequencies."
freq_variants_job_id=$(sbatch --parsable --dependency=afterok:${concat_variants_job_id} --time ${PREP_TIME} --mem ${PREP_MEM} --mail-user ${PREP_MAIL} --mail-type ${PREP_MAILTYPE} ${PGSTK}/pgstoolkit.prep.freq.run.sh $PLINK $STUDYDIR)

# Deactivate the conda environment
conda deactivate