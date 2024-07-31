#!/bin/bash
#################################################################################################
### PARAMETERS SLURM
#SBATCH --job-name=runPGSTK                                  		# the name of the job
#SBATCH --output=/hpc/dhl_ec/svanderlaan/projects/runPGSTK.log 	    # the log file of this job
#SBATCH --error=/hpc/dhl_ec/svanderlaan/projects/runPGSTK.errors	# the error file of this job
#SBATCH --time=02:15:00                                             # the amount of time the job will take: -t [min] OR -t [days-hh:mm:ss]
#SBATCH --mem=64G                                                   # the amount of memory you think the script will consume, found on: https://wiki.bioinformatics.umcutrecht.nl/bin/view/HPC/SlurmScheduler
#SBATCH --gres=tmpspace:128G                                        # the amount of temporary diskspace per node
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

### Set the project directory
PROJECTDIR="/hpc/dhl_ec/svanderlaan/projects/"
PGSDIR="/hpc/dhl_ec/svanderlaan/projects/polygenicscores"
PGS_CAD="${PGSDIR}/Inouye_bioRxiv_2018"
PGS_CAD_UKBB="${PGSDIR}/UKBB_GWAS1KG_2017"

SOFTWARE="/hpc/local/Rocky8/dhl_ec/software"

BGENIX="/hpc/local/Rocky8/dhl_ec/bin/bgenix"
PLINK="${SOFTWARE}/plink2_linux_x86_64_20240105_alpha_5_10/plink2"

STUDYDATADIR="/hpc/dhl_ec/data/_ae_originals/"
### b38 -- TOPMed imputed
STUDYDIR="${STUDYDATADIR}/AEGS_QC_imputation_2023/aegscombo/_topmed_r3_f10_b38"
### b37 version -- deprecated
### STUDYDIR="${STUDYDATADIR}/AEGS_COMBINED_EAGLE2_1000Gp3v5HRCr11"

echo "========================================================================"
echo "STEP 2 PGSToolKit pipeline"
echo "========================================================================"
echo ""
echo "Note that this part is specific to a given polygenic score to be created."
echo ""

### First, we need to create the input for PGSToolKit

# # Head INOUYE
# # chr	position	rsid	allele1	allele2	effect_allele	beta
# #  1	  2245570	rs2843152	C	G	G	-2.76009e-02
# #  1	 22132518	rs35465346	A	G	G	 2.39340e-02
# #  1	 38386727	rs28470722	A	G	G	-1.74935e-02
# 
echo "variantid rsid rsid_aegs chromosome position effect_allele other_allele beta P_ukbb" > ${PGS_CAD}/metaGRS_hg19_20180205.foo
zcat ${PGS_CAD}/metaGRS_hg19_20180205.txt.gz | \
parseTable --col chr,position,rsid,allele1,allele2,effect_allele,beta | \
awk '{ if($6 == $5) { print $1":"$2, $3, $3, $1, $2, $6, $4, $7, "NA" } else { print $1":"$2, $3, $3, $1, $2, $6, $5, $7, "NA" } }' | tail -n +2 >> ${PGS_CAD}/metaGRS_hg19_20180205.foo

mergeTablesv2 \
--file1 ${PGS_CAD_UKBB}/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.4pvalupdate.txt \
--file2 ${PGS_CAD}/metaGRS_hg19_20180205.foo \
--index variantid --format NORM --replace > ${PGS_CAD}/metaGRS_hg19_20180205.4PGSTK.foo

mergeTablesv2 \
--file1 ${STUDYDIR}/aegs.qc.1kgp3hrcr11.idfix.rsid.8bit.variantlist.txt.gz \
--file2 ${PGS_CAD}/metaGRS_hg19_20180205.4PGSTK.foo \
--index variantid --format GZIP1 --replace > ${PGS_CAD}/metaGRS_hg19_20180205.4PGSTK.txt

### Second, we need to create the configuration file for PGSToolKit
### There is a template available in the PGSToolKit directory, and we will copy this to create one 
### for our specific polygenic score example for CAD.
### In this example, the file is called "${PROJECTDIR}/PGS/pgstoolkit.cad.config".

### Third, we can run PGSToolKit
echo ""
echo "Running PGSToolKit."

#################################################################################################
### PARAMETERS SLURM YOU SHOULD PROVIDE
### --job-name=pgsTK												# the name of the job
### --output=/hpc/dhl_ec/svanderlaan/projects/pgsTK.log 			# the log file of this job
### --error=/hpc/dhl_ec/svanderlaan/projects/pgsTK.errors			# the error file of this job
### --time=12:15:00													# the amount of time the job will take: -t [min] OR -t [days-hh:mm:ss]
### --mem=48G														# the amount of memory you think the script will consume, found on: https://wiki.bioinformatics.umcutrecht.nl/bin/view/HPC/SlurmScheduler
### --gres=tmpspace:128G											# the amount of temporary diskspace per node
### --mail-user=s.w.vanderlaan-2@umcutrecht.nl									# where should be mailed to?
### --mail-type=FAIL												# when do you want to receive a mail from your job?  Valid type values are NONE, BEGIN, END, FAIL, REQUEUE
																	# or ALL (equivalent to BEGIN, END, FAIL, INVALID_DEPEND, REQUEUE, and STAGE_OUT), 
																	# Multiple type values may be specified in a comma separated list. 
####	Note:	You do not have to specify workdir: 
####			'Current working directory is the calling process working directory unless the --chdir argument is passed, which will override the current working directory.'
####			TODO: select the type of interpreter you'd like to use
####			TODO: Find out whether this job should dependant on other scripts (##SBATCH --depend=[state:job_id])
####
#################################################################################################

### Submit the job for PGSToolKit -- uncomment after you have prepared the input and created the configuration file (see above)
# sbatch --job-name=pgsTK_CAD --output=${PROJECTDIR}/PGS/pgsTK_CAD.log --error=${PROJECTDIR}/PGS/pgsTK_CAD.errors  \
# --time=12:15:00 --mem=48G --gres=tmpspace:128G \
# --mail-user=s.w.vanderlaan-2@umcutrecht.nl --mail-type=FAIL \
# ${SOFTWARE}/PGSToolKit/pgstoolkit.sh ${PROJECTDIR}/PGS/pgstoolkit.cad.config

