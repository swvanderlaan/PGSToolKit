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

PROJECTDIR="/hpc/dhl_ec/svanderlaan/projects/"
PGSDIR="/hpc/dhl_ec/svanderlaan/projects/polygenicscores"
PGS_CAD="${PGSDIR}/Inouye_bioRxiv_2018"
PGS_CAD_UKBB="${PGSDIR}/UKBB_GWAS1KG_2017"

SOFTWARE="/hpc/local/CentOS7/dhl_ec/software"

BGENIX="/hpc/local/CentOS7/dhl_ec/bin/bgenix"
PLINK="${SOFTWARE}/plink2_linux_x86_64_20200124_alpha_v2.3_final/plink2"

AEGSDATADIR="/hpc/dhl_ec/data/_ae_originals/"
AEGSDIR="${AEGSDATADIR}/AEGS_COMBINED_EAGLE2_1000Gp3v5HRCr11"

echo ""
echo "Converting VCF to bgen."
### This data includes the correct ID-type (chr:BP, e.g. 1:10711)
# for CHR in $(seq 1 22 ); do
#   
#   echo "> converting chromosome ${CHR}..."
#   echo "...8-bits version and indexing"
#   $PLINK --vcf ${AEGSDIR}/aegs.qc.1kgp3hrcr11.idfix.chr${CHR}.vcf.gz --export bgen-1.2 bits=8 --out ${AEGSDIR}/aegs.qc.1kgp3hrcr11.idfix.8bit.chr${CHR}
#   $BGENIX -index -g ${AEGSDIR}/aegs.qc.1kgp3hrcr11.idfix.8bit.chr${CHR}.bgen 
#   
#   echo ""
#   echo "...16-bits version (default) and indexing"
#   $PLINK --vcf ${AEGSDIR}/aegs.qc.1kgp3hrcr11.idfix.chr${CHR}.vcf.gz --export bgen-1.2 --out ${AEGSDIR}/aegs.qc.1kgp3hrcr11.idfix.16bit.chr${CHR}
#   $BGENIX -index -g ${AEGSDIR}/aegs.qc.1kgp3hrcr11.idfix.16bit.chr${CHR}.bgen
#   
# done

echo ""
echo "Creating variant lists."
# for CHR in $(seq 1 22); do 
# 	bgenix_20230202 \
# 	-g ${AEGSDIR}/aegs.qc.1kgp3hrcr11.idfix.rsid.8bit.chr${CHR}.bgen \
# 	-i ${AEGSDIR}/aegs.qc.1kgp3hrcr11.idfix.rsid.8bit.chr${CHR}.bgen.bgi \
# 	-list | grep -v "#" > ${AEGSDIR}/aegs.qc.1kgp3hrcr11.idfix.rsid.8bit.chr${CHR}.variantlist.txt; 
# done

# Variantlist looks like this
# alternate_ids	rsid	chromosome	position	number_of_alleles	first_allele	alternative_alleles
# rs367896724	rs367896724	1	10177	2	A	AC
# rs540431307	rs540431307	1	10235	2	T	TA
# rs555500075	rs555500075	1	10352	2	T	TA
# rs537182016	rs537182016	1	10539	2	C	A

echo ""
echo "Concatenating variant list."

# echo "variantid alt_ids rsid_aegs chromosome position number_of_alleles first_allele alt_alleles" > ${AEGSDIR}/aegs.qc.1kgp3hrcr11.idfix.rsid.8bit.variantlist.txt
# 
# for CHR in $(seq 1 22); do 
# 	cat ${AEGSDIR}/aegs.qc.1kgp3hrcr11.idfix.rsid.8bit.chr${CHR}.variantlist.txt | \
# 	awk '{ print $3":"$4, $1, $2, $3, $4, $5, $6, $7}' | tail -n +2 >> ${AEGSDIR}/aegs.qc.1kgp3hrcr11.idfix.rsid.8bit.variantlist.txt 
# done
# gzip -v ${AEGSDIR}/aegs.qc.1kgp3hrcr11.idfix.rsid.8bit.variantlist.txt

# Concatenated variantlist looks like this
# variantid alt_ids rsid_aegs chromosome position number_of_alleles first_allele alt_alleles
# 1:10177 rs367896724 rs367896724 1 10177 2 A AC
# 1:10235 rs540431307 rs540431307 1 10235 2 T TA
# 1:10352 rs555500075 rs555500075 1 10352 2 T TA
# 1:10539 rs537182016 rs537182016 1 10539 2 C A
# 1:10616 rs376342519 rs376342519 1 10616 2 CCGCCGTTGCAAAGGCGCGCCG C

echo ""
echo "Calculating frequencies."

# for CHR in $(seq 1 22); do 
# 	echo "> calculating frequencies for chromosome ${CHR}..."
# 	$PLINK \
# 	--bgen ${AEGSDIR}/aegs.qc.1kgp3hrcr11.idfix.rsid.8bit.chr${CHR}.bgen 'ref-first' \
# 	--sample ${AEGSDATADIR}//AEGS_COMBINED_EAGLE2_1000Gp3v5HRCr11_pheno/20200419.QC.AEGS123.forPGSToolKit.sample \
# 	--freq --out ${AEGSDIR}/aegs.qc.1kgp3hrcr11.idfix.rsid.8bit.chr${CHR}.FREQ; 
# done

echo ""
echo "Creating input for PGSToolKit."

# # Head INOUYE
# # chr	position	rsid	allele1	allele2	effect_allele	beta
# #  1	  2245570	rs2843152	C	G	G	-2.76009e-02
# #  1	 22132518	rs35465346	A	G	G	 2.39340e-02
# #  1	 38386727	rs28470722	A	G	G	-1.74935e-02
# 
# echo "variantid rsid rsid_aegs chromosome position effect_allele other_allele beta P_ukbb" > ${PGS_CAD}/metaGRS_hg19_20180205.foo
# zcat ${PGS_CAD}/metaGRS_hg19_20180205.txt.gz | \
# parseTable --col chr,position,rsid,allele1,allele2,effect_allele,beta | \
# awk '{ if($6 == $5) { print $1":"$2, $3, $3, $1, $2, $6, $4, $7, "NA" } else { print $1":"$2, $3, $3, $1, $2, $6, $5, $7, "NA" } }' | tail -n +2 >> ${PGS_CAD}/metaGRS_hg19_20180205.foo
# 
# mergeTablesv2 \
# --file1 ${PGS_CAD_UKBB}/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.4pvalupdate.txt \
# --file2 ${PGS_CAD}/metaGRS_hg19_20180205.foo \
# --index variantid --format NORM --replace > ${PGS_CAD}/metaGRS_hg19_20180205.4PGSTK.foo
# 
# mergeTablesv2 \
# --file1 ${AEGSDIR}/aegs.qc.1kgp3hrcr11.idfix.rsid.8bit.variantlist.txt.gz \
# --file2 ${PGS_CAD}/metaGRS_hg19_20180205.4PGSTK.foo \
# --index variantid --format GZIP1 --replace > ${PGS_CAD}/metaGRS_hg19_20180205.4PGSTK.txt

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
### --mail-user=yourname@domain.com									# where should be mailed to?
### --mail-type=FAIL												# when do you want to receive a mail from your job?  Valid type values are NONE, BEGIN, END, FAIL, REQUEUE
																	# or ALL (equivalent to BEGIN, END, FAIL, INVALID_DEPEND, REQUEUE, and STAGE_OUT), 
																	# Multiple type values may be specified in a comma separated list. 
####	Note:	You do not have to specify workdir: 
####			'Current working directory is the calling process working directory unless the --chdir argument is passed, which will override the current working directory.'
####			TODO: select the type of interpreter you'd like to use
####			TODO: Find out whether this job should dependant on other scripts (##SBATCH --depend=[state:job_id])
####
#################################################################################################

### INOUYE CAD
sbatch --job-name=pgsTK_CAD --output=${PROJECTDIR}/PGS/pgsTK_CAD.log --error=${PROJECTDIR}/PGS/pgsTK_CAD.errors  \
--time=12:15:00 --mem=48G --gres=tmpspace:128G \
--mail-user=s.w.vanderlaan-2@umcutrecht.nl --mail-type=FAIL \
${SOFTWARE}/PGSToolKit/pgstoolkit.sh ${PROJECTDIR}/PGS/pgstoolkit.cad.config

