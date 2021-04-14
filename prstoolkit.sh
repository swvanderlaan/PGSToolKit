#!/bin/bash
#SBATCH --job-name="PRS workflow"                                 #Job name
#SBATCH --output=/home/dhl_ec/aligterink/logs/batchjob_%j.log     # Standard output and error log
#SBATCH --mail-user=a.j.ligterink@umcutrecht.nl                   # Mail
#SBATCH --mail-type=NONE		                                  # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --time=20:00:00                                           # Time limit hrs:min:sec
#SBATCH --mem=4G	                                              # RAM required per node

set -- "${@:1:2}" "/home/dhl_ec/aligterink/PRSToolKit/prstoolkit.config" "/hpc/dhl_ec/aligterink/ProjectFiles/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.SNPIDFIX2.txt.gz"

# Setting colouring
NONE='\033[00m'
OPAQUE='\033[2m'
FLASHING='\033[5m'
BOLD='\033[1m'
ITALIC='\033[3m'
UNDERLINE='\033[4m'
STRIKETHROUGH='\033[9m'

RED='\033[01;31m'
GREEN='\033[01;32m'
YELLOW='\033[01;33m'
PURPLE='\033[01;35m'
CYAN='\033[01;36m'
WHITE='\033[01;37m'

### Regarding changing the 'type' of the things printed with 'echo'
### Refer to: 
### - http://askubuntu.com/questions/528928/how-to-do-underline-bold-italic-strikethrough-color-background-and-size-i
### - http://misc.flogisoft.com/bash/tip_colors_and_formatting
### - http://unix.stackexchange.com/questions/37260/change-font-in-echo-command

# Creating color-functions
function echobold { #'echobold' is the function name
    echo -e "${BOLD}${1}${NONE}" # this is whatever the function needs to execute, note ${1} is the text for echo
}
function echoitalic { 
    echo -e "${ITALIC}${1}${NONE}" 
}
function echonooption { 
    echo -e "${OPAQUE}${RED}${1}${NONE}"
}
function echoerrorflash { 
    echo -e "${RED}${BOLD}${FLASHING}${1}${NONE}" 
}
function echoerror { 
    echo -e "${RED}${1}${NONE}"
}
function echocyan { #'echobold' is the function name
    echo -e "${CYAN}${1}${NONE}" # this is whatever the function needs to execute.
}
function echoerrornooption { 
    echo -e "${YELLOW}${1}${NONE}"
}
function echoerrorflashnooption { 
    echo -e "${YELLOW}${BOLD}${FLASHING}${1}${NONE}"
}

script_copyright_message() {
	echo ""
	THISYEAR=$(date +'%Y')
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	echo "+ The MIT License (MIT)                                                                                 +"
	echo "+ Copyright (c) 2016-${THISYEAR} Sander W. van der Laan                                                        +"
	echo "+                                                                                                       +"
	echo "+ Permission is hereby granted, free of charge, to any person obtaining a copy of this software and     +"
	echo "+ associated documentation files (the \"Software\"), to deal in the Software without restriction,         +"
	echo "+ including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, +"
	echo "+ and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, +"
	echo "+ subject to the following conditions:                                                                  +"
	echo "+                                                                                                       +"
	echo "+ The above copyright notice and this permission notice shall be included in all copies or substantial  +"
	echo "+ portions of the Software.                                                                             +"
	echo "+                                                                                                       +"
	echo "+ THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT     +"
	echo "+ NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND                +"
	echo "+ NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES  +"
	echo "+ OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN   +"
	echo "+ CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                            +"
	echo "+                                                                                                       +"
	echo "+ Reference: http://opensource.org.                                                                     +"
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
}

script_arguments_error() {
	echoerror "$1" # Additional message
	echoerror "- Argument #1 is path_to/filename of the configuration file."
	echoerror "- Argument #2 is path_to/filename of the list of GWAS summary statistics-files with arbitrarily chosen names, path_to, and file-names."
	echoerror ""
	echoerror "An example command would be: prstoolkit.sh [arg1] [arg2]"
	echoerror "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
 	echo ""
	script_copyright_message
	exit 1
}

echobold "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echobold "                          PRSToolKit: A TOOLKIT TO CALCULATE POLYGENIC SCORES"
echoitalic "                        --- prepare data and calculate polygenic scores ---"
echobold ""
echobold "* Version:      v1.0.1"
echobold ""
echobold "* Last update:  2021-04-07"
echobold "* Written by:   Sander W. van der Laan | s.w.vanderlaan@gmail.com"
echobold "				  Anton Ligterink | anton.ligterink@gmail.com"
echobold "* Description:  Prepare files and calculate polygenic scores. This tool will do the following:"
echobold "                - Perform quality control on GWAS summary data"
echobold "                - Calculate polygenic risk scores from GWAS summary statistcs (LDpred/PRSICE/RapidoPGS/PRScs)"
echobold "				  - Calculate polygenic scores using premade scoring systems (PLINK)"
echobold "* References:   RapidoPGS: https://cran.r-project.org/web/packages/RapidoPGS/vignettes/Computing_RapidoPGS.html"
echobold "                LDpred: https://github.com/bvilhjal/ldpred      PRSICE: https://www.prsice.info/"
echobold "                PRScs: https://github.com/getian107/PRScs       PLINK: http://zzz.bwh.harvard.edu/plink/"
echobold "* REQUIRED: "
echobold "  - A high-performance computer cluster with a SLURM system"
echobold "  - R v3.6+, Python 3.7+"
echobold "  - Required R 3.6+ modules: [RapidoPGS], [data.table], [optparse]"
echobold "  - Required Python 3.7+ modules: [pandas], [scipy], [numpy], [plinkio], [h5py]"
### ADD-IN: function to check requirements...
### This might be a viable option! https://gist.github.com/JamieMason/4761049
echobold ""
echobold "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Today's: "$(date)
TODAY=$(date +"%Y%m%d")

##########################################################################################
### SET THE SCENE FOR THE SCRIPT
##########################################################################################

##########################################################################################
### Command-line argument check
if [[ $# -lt 2 ]]; then 
	echo ""
	echoerror "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	echoerrorflash "               *** Oh, oh, computer says no! Number of arguments found "$#". ***"
	echoerror "You must supply [2] arguments when running *** PRSToolKit ***!"
	script_arguments_error
else
	echo "These are the "$#" arguments that passed:"
	echo "The configuration file..................: "$(basename ${1}) # argument 1
	echo "The GWAS-files file.....................: "$(basename ${2}) # argument 2
fi 

##########################################################################################
### Loading command-line arguments and configuration file
source "$1" # Depends on arg1.
CONFIGURATIONFILE="$1" # Depends on arg1 -- but also on where it resides!!!
BASEDATA="$2" # Depends on arg2 -- all the GWAS dataset information

##########################################################################################
### Parameter validation
if [[ ${PRSMETHOD} == "PLINK" || ${PRSMETHOD} == "PRSCS" || ${PRSMETHOD} == "LDPRED" || ${PRSMETHOD} == "PRSICE2" || ${PRSMETHOD} == "RAPIDOPGS" ]]; then
	echo ""
	echo "Calculating Polygenic Risk Scores using ${PRSMETHOD}."

elif [[ ${PRSMETHOD} == "MANUAL" ]]; then
	echo ""
	echoerrornooption "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	echoerrornooption ""
	echoerrorflashnooption "               *** Oh, computer says no! This option is not available yet. ***"
	echoerrornooption "Unfortunately using [${PRSMETHOD}] as a PRS calculation method is not possible yet."
	echoerrornooption "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	### The wrong arguments are passed, so we'll exit the script now!
	echo ""
	script_copyright_message
	exit 1

else
	echo ""
	echoerror "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	echoerror ""
	echoerrorflash "                  *** Oh, computer says no! Argument not recognised. ***"
	echoerror "You have the following options to calculate polygenic scores:"
	echoerror "LDPRED   -- use LDpred to LD-correct GWAS summary statistics and calculate polygenic scores [default]."
	echonooption "MANUAL   -- use an in-house developed Rscript to calculate a polygenic score using a limited set of variants."
	echonooption "PLINK    -- use PLINK to calculate scores, GWAS summary statistics will be LD-pruned using p-value thresholds (traditional approach, slow)."
	echonooption "PRSICE   -- use PRSice to calculate scores, GWAS summary statistics will be LD-pruned using p-value thresholds (new approach, fast)."
	echonooption "(Opaque: *not implemented yet*)"
	echoerror "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	### The wrong arguments are passed, so we'll exit the script now!
	echo ""
	script_copyright_message
	exit 1
fi


##########################################################################################
### SETTING UP NECESSARY DIRECTORIES
echo ""
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Checking for the existence of the output directory [ ${OUTPUTDIRNAME} ]."
if [ ! -d ${PROJECTDIR}/${OUTPUTDIRNAME} ]; then
	echo "> Output directory doesn't exist - Mr. Bourne will create it for you."
	mkdir -v ${PROJECTDIR}/${OUTPUTDIRNAME}
else
	echo "> Output directory already exists."
fi
OUTPUTDIR=${PROJECTDIR}/${OUTPUTDIRNAME}

echo ""
echo "Checking for the existence of the subproject directory [ ${SUBPROJECTDIRNAME} ]."
if [ ! -d ${OUTPUTDIR}/${SUBPROJECTDIRNAME} ]; then
	echo "> Subproject directory doesn't exist - Mr. Bourne will create it for you."
	mkdir -v ${OUTPUTDIR}/${SUBPROJECTDIRNAME}
else
	echo "> Subproject directory already exists."
fi
SUBPROJECTDIR=${OUTPUTDIR}/${SUBPROJECTDIRNAME}

echo ""
echo "Checking for the existence of the parsed data directory [ ${SUBPROJECTDIRNAME}/PARSED ]."
if [ ! -d ${SUBPROJECTDIR}/PARSED ]; then
	echo "> Parsed data directory doesn't exist - Mr. Bourne will create it for you."
	mkdir -v ${SUBPROJECTDIR}/PARSED
else
	echo "> Parsed data directory already exists."
fi
PARSEDDIR=${SUBPROJECTDIR}/PARSED

# echo ""
# echo "Checking for the existence of the ${PRSMETHOD} directory [ ${PROJECTDIR}/${PRSMETHOD} ]."
# if [ ! -d ${PROJECTDIR}/${PRSMETHOD} ]; then
# 	echo "> ${PRSMETHOD} directory doesn't exist - Mr. Bourne will create it for you."
# 	mkdir -v ${PROJECTDIR}/${PRSMETHOD}
# else
# 	echo "> ${PRSMETHOD} directory already exists."
# fi
# PRSDIR=${PROJECTDIR}/${PRSMETHOD}

echo ""
echo "Checking for the existence of the ${TEMPDIRNAME} directory [ ${PROJECTDIR}/${TEMPDIRNAME} ]."
if [ ! -d ${PROJECTDIR}/${TEMPDIRNAME} ]; then
	echo "> ${TEMPDIRNAME} directory doesn't exist - Mr. Bourne will create it for you."
	mkdir -v ${PROJECTDIR}/${TEMPDIRNAME}
else
	echo "> ${TEMPDIRNAME} directory already exists."
fi
PRSDIR=${PROJECTDIR}/${TEMPDIRNAME}

echo ""
echo "Checking for the existence of the ${LOGDIRNAME} directory [ ${PROJECTDIR}/${LOGDIRNAME} ]."
if [ ! -d ${PROJECTDIR}/${LOGDIRNAME} ]; then
	echo "> ${LOGDIRNAME} directory doesn't exist - Mr. Bourne will create it for you."
	mkdir -v ${PROJECTDIR}/${LOGDIRNAME}
else
	echo "> ${LOGDIRNAME} directory already exists."
fi
LOGDIR=${PROJECTDIR}/${LOGDIRNAME}

OUTPUTNAME=${PROJECTNAME}_${PRSMETHOD}_$(date +%Y-%b-%d--%H-%M)

echo ""
echo "The scene is properly set, and directories are created! 🖖"
echo "PRSToolKit directory............................................: "${PRSTOOLKITDIR}
echo "LD reference data directory.....................................: "${LDDATA}
echo "Validation data directory.......................................: "${VALIDATIONDATA}
echo "Main directory..................................................: "${PROJECTDIR}
echo "Main analysis output directory..................................: "${OUTPUTDIR}
echo "Subproject's analysis output directory..........................: "${SUBPROJECTDIR}
echo "Parsed data directory...........................................: "${PARSEDDIR}
echo "Working directory...............................................: "${PRSDIR}
echo "We are processing these GWAS-summary statistics(s)..............: "${BASEDATA}
# while IFS='' read -r GWASCOHORT || [[ -n "$GWASCOHORT" ]]; do
# 	LINE=${GWASCOHORT}
# 	COHORT=$(echo "${LINE}" | awk '{ print $1 }')
# 	echo " * ${COHORT}"
# done < ${BASEDATA}
echo ""
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

# check_slurm_codes() {
# 	echo "${CURRENT_JOBID}"
# 	echo "$(sacct -j ${CURRENT_JOBID})"
# }
QC_DEPENDENCY=""
if [[ ${QC} == "YES" ]]; then

	echo ""
	echobold "#========================================================================================================"
	echobold "#== OPTIONAL QUALITY CONTROL IS IN EFFECT [DEFAULT]"
	echobold "#========================================================================================================"
	echobold "#"
	echo ""
	echo "We will perform quality control on the base dataset, using the following thresholds:"
	echo "      Minor allele frequency: ${MAF}"
	echo "      Imputation score:       ${INFO}"

	# Start the quality control job
	QC_RESULTS_FILE=${PRSDIR}/QC_results.txt
	QC_JOBID=$(sbatch --parsable --wait --job-name=PRS_QC --time ${RUNTIME_QC} --mem ${MEMORY_QC} -o ${LOGDIR}/${OUTPUTNAME}_QC.log ${PRSTOOLKITSCRIPTS}/QC.R -b ${BASEDATA} -s ${STATS_FILE} -o ${PRSDIR}/QCd_basefile.txt.gz -m ${MAF} -i ${INFO} -a ${BF_SNP_COL} -d ${STATS_ID_COL} -c ${STATS_MAF_COL} -t ${STATS_INFO_COL} -r ${QC_RESULTS_FILE})
	QC_DEPENDENCY="--dependency=afterany:${QC_JOBID}"
	
	echo ""
	cat ${QC_RESULTS_FILE}
	echo ""

	# We will now use the quality controlled file as our new base file
	BASEDATA=${PRSDIR}/QCd_basefile.txt.gz

elif [[ ${QC} == "NO" ]]; then
	echo ""
	echo "Quality control will not be performed"
	echo ""
else 
	echo ""
	echo "QC parameter should be [YES/NO], not \"${QC}\""
	echo ""
	exit 1
fi

echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo ""
echo "Risk score computation will now commence."
echo ""
##########################################################################################
### Run the individual PRS method scripts

RESULTS_FILE=${OUTPUTDIR}/${OUTPUTNAME}_RESULTS.txt

if [[ ${PRSMETHOD} == "PLINK" ]]; then

	# Calculate the effect sizes for each chromosome
	### TODO: make this run in parallel and make it prettier
	PLINKSCORE_JOBID=$(sbatch --parsable --wait --job-name=PRS_PLINKSCORE ${QC_DEPENDENCY} --time ${RUNTIME_PLINKSCORE} --mem ${MEMORY_PLINKSCORE} -o ${LOGDIR}/${OUTPUTNAME}_PLINK_score.log --export=ALL,VALIDATIONDATA=${VALIDATIONDATA},VALIDATIONPREFIX=${VALIDATIONPREFIX},PLINK=${PLINK},PLINK_REF_POS=${PLINK_REF_POS},SAMPLE_FILE=${SAMPLE_FILE},PRSDIR=${PRSDIR},WEIGHTS_FILE=${BASEDATA},PLINK_VARIANT_ID_COL=${PLINK_VARIANT_ID_COL},PLINK_ALLELE_COL=${PLINK_ALLELE_COL},PLINK_SCORE_COL=${PLINK_SCORE_COL},PLINK_HEADER=${PLINK_HEADER},PLINK_SETTINGS=${PLINK_SETTINGS} ${PRSTOOLKITSCRIPTS}/plinkscore.sh)
	PLINKSCORE_DEPENDENCY="--dependency=afterany:${PLINKSCORE_JOBID}"

	# Sum the effect sizes to calculate the final score
	PLINKSUM_JOBID=$(sbatch --parsable --wait --job-name=PRS_PLINKSCORE ${PLINKSCORE_DEPENDENCY} --time ${RUNTIME_PLINKSUM} --mem ${MEMORY_PLINKSUM} -o ${LOGDIR}/${OUTPUTNAME}_PLINK_sum.log ${PRSTOOLKITSCRIPTS}/sum_plink_scores.R -s ${SAMPLE_FILE} -f ${SAMPLE_FID_COL} -i ${SAMPLE_IID_COL} -d ${PRSDIR} -p ${PLINKSCORE_PREFIX} -f 1 -i 2 -r 6 -o ${RESULTS_FILE})


elif [[ ${PRSMETHOD} == "PRSCS" ]]; then
	# source "${PRSTOOLKITSCRIPTS}/prscs.sh"
	echo "Not implemented yet!"
	 
elif [[ ${PRSMETHOD} == "LDPRED" ]]; then
	echo "Not implemented yet!"

elif [[ ${PRSMETHOD} == "PRSICE2" ]]; then
	# source "${PRSTOOLKITSCRIPTS}/prsice.sh"
	echo "Not implemented yet!"


elif [[ ${PRSMETHOD} == "RAPIDOPGS" ]]; then

	#Set PLINK parameters
	PLINK_VARIANT_ID_COL=1
	PLINK_ALLELE_COL=2
	PLINK_SCORE_COL=3
	PLINK_HEADER="header"
	WEIGHTS_FILE="${PRSDIR}/${OUTPUTNAME}_weights.txt"

	# Calculate weights using RapidoPGS
	RAPIDO_JOBID=$(sbatch --parsable --wait --job-name=PRS_RAPIDO ${QC_DEPENDENCY} --time ${RUNTIME_RAPIDO} --mem ${MEMORY_RAPIDO} -o ${LOGDIR}/${OUTPUTNAME}_RAPIDO.log ${PRSTOOLKITSCRIPTS}/rapidopgs.R -k ${PRSDIR} -o ${WEIGHTS_FILE} -b ${BASEDATA} -d ${BF_BUILD} -i ${BF_SNP_COL} -c ${BF_CHR_COL} -p ${BF_POS_COL} -r ${BF_A1_COL} -a ${BF_A2_COL} -f ${BF_FRQ_COL} -w ${BF_WFRQ} -m ${BF_MSR_COL} -e ${BF_SE_COL} -s ${BF_SBJ_COL})
	RAPIDO_DEPENDENCY="--dependency=afterany:${RAPIDO_JOBID}"

	# ${RSCRIPTPATH} ${PRSTOOLKITSCRIPTS}/rapidopgs.R -k ${PRSDIR} -o ${SCORE_FILE} -b ${BASEDATA} -d ${BF_BUILD} -i ${BF_SNP_COL} \
	# -c ${BF_CHR_COL} -p ${BF_POS_COL} -r ${BF_A1_COL} -a ${BF_A2_COL} -f ${BF_FRQ_COL} -w ${BF_WFRQ} -m ${BF_MSR_COL} -e ${BF_SE_COL} -s ${BF_SBJ_COL}

	
	# Calculate individual scores using PLINK
	PLINKSCORE_JOBID=$(sbatch --parsable --wait --job-name=PRS_PLINKSCORE ${RAPIDO_DEPENDENCY} --time ${RUNTIME_PLINKSCORE} --mem ${MEMORY_PLINKSCORE} -o ${LOGDIR}/${OUTPUTNAME}_PLINK_score.log --export=ALL,VALIDATIONDATA=${VALIDATIONDATA},VALIDATIONPREFIX=${VALIDATIONPREFIX},PLINK=${PLINK},PLINK_REF_POS=${PLINK_REF_POS},SAMPLE_FILE=${SAMPLE_FILE},PRSDIR=${PRSDIR},WEIGHTS_FILE=${BASEDATA},PLINK_VARIANT_ID_COL=${PLINK_VARIANT_ID_COL},PLINK_ALLELE_COL=${PLINK_ALLELE_COL},PLINK_SCORE_COL=${PLINK_SCORE_COL},PLINK_HEADER=${PLINK_HEADER},PLINK_SETTINGS=${PLINK_SETTINGS} ${PRSTOOLKITSCRIPTS}/plinkscore.sh)
	PLINKSCORE_DEPENDENCY="--dependency=afterany:${PLINKSCORE_JOBID}"

	# Sum the effect sizes to calculate the final score
	PLINKSUM_JOBID=$(sbatch --parsable --wait --job-name=PRS_PLINKSCORE ${PLINKSCORE_DEPENDENCY} --time ${RUNTIME_PLINKSUM} --mem ${MEMORY_PLINKSUM} -o ${LOGDIR}/${OUTPUTNAME}_PLINK_sum.log ${PRSTOOLKITSCRIPTS}/sum_plink_scores.R -s ${SAMPLE_FILE} -f ${SAMPLE_FID_COL} -i ${SAMPLE_IID_COL} -d ${PRSDIR} -p ${PLINKSCORE_PREFIX} -f 1 -i 2 -r 6 -o ${RESULTS_FILE})

	# source "${PRSTOOLKITSCRIPTS}/plinkscore.sh"
	# ${RSCRIPTPATH} ${PRSTOOLKITSCRIPTS}/sum_plink_scores.R
fi

# Make all created files accessible
# for file in ${PRSDIR}/*; do
#     # readlink -f $file
# 	chmod -Rv a+rwx $file
#     #chmod -Rv a+rwx ${file%.vcf.gz}.bgen
# done
echo ""
echo "${PRSMETHOD} risk score calculation has finished."
echo "The risk scores were stored in [ ${RESULTS_FILE} ]"
echo ""
echocyan "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echobold "Wow. I'm all done buddy. What a job 😱 ! let's have a 🍻🍻 ... 🖖 "







# ${RSCRIPTPATH} /home/dhl_ec/aligterink/scripts/QC_testing.R
# ${RSCRIPTPATH} ${PRSTOOLKITSCRIPTS}/QC.R


# ${RSCRIPTPATH} /home/dhl_ec/aligterink/PRSToolKit/SCRIPTS/Loose_scripts/vcf_chrbp_to_rs.R

# QCTOOL="/hpc/local/CentOS7/dhl_ec/software/qctool_v2.0.8/build/release/qctool_v2.0.8"
# VCFFILE="/hpc/dhl_ec/data/_ae_originals/AEGS_COMBINED_EAGLE2_1000Gp3v5HRCr11/aegs.qc.1kgp3hrcr11.idfix.chr21.vcf.gz"
# MAPFILE="/hpc/dhl_ec/aligterink/ProjectFiles/vcf_qctool_files/mapfiles/chr21_whole2.txt"

# OUTPUT="/hpc/dhl_ec/aligterink/ProjectFiles/vcf_qctool_files/outputfiles/chr21_whole2.vcf.gz"
# ${QCTOOL} -g ${VCFFILE} -og ${OUTPUT} -map-id-data ${MAPFILE} -compare-variants-by position

# file=/hpc/dhl_ec/aligterink/ProjectFiles/vcf_qctool_files/outputfiles/chr21_whole2.vcf.gz
# out=/hpc/dhl_ec/aligterink/ProjectFiles/vcf_qctool_files/outputfiles/chr21_whole2.bgen
# /hpc/local/CentOS7/dhl_ec/software/qctool_v204 -g $file -vcf-genotype-field GP -og $out

# ${RSCRIPTPATH} /home/dhl_ec/aligterink/PRSToolKit/SCRIPTS/Loose_scripts/sumstats_summary.R







# if [[ ${VALIDATIONFORMAT} == "VCF" ]]; then
# 	echo ""
# 	echo "The validation dataset is encoded in the [${VALIDATIONFORMAT}] file-format; PRSToolKit will procede "
# 	echo "immediately after optional QC."

# elif [[ ${VALIDATIONFORMAT} == "OXFORD" || ${VALIDATIONFORMAT} == "PLINK" ]]; then

# 	echoerrornooption "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
# 	echoerrornooption ""
# 	echoerrorflashnooption "               *** Oh, computer says no! This option is not available yet. ***"
# 	echoerrornooption "Unfortunately the [${VALIDATIONFORMAT}] file-format is not supported."
# 	echoerrornooption "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
# 	### The wrong arguments are passed, so we'll exit the script now!
# 	echo ""
# 	script_copyright_message
# 	exit 1
	
# else
# 	echoerror "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
# 	echoerror ""
# 	echoerrorflash "                  *** Oh, computer says no! Argument not recognised. ***"
# 	echoerror "You can indicate the following validation file-formats:"
# 	echoerror "OXFORD   -- file format used by IMPUTE2 (.gen/.bgen) [default]."
# 	echonooption "VCF   -- VCF file format, version 4.2 is expected."
# 	echoerror "PLINK    -- PLINK file format; PRSToolKit can immediately use this."
# 	echonooption "(Opaque: *not implemented yet*)"
# 	echoerror "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
# 	### The wrong arguments are passed, so we'll exit the script now!
# 	echo ""
# 	script_copyright_message
# 	exit 1
# fi
##########################################################################################
### SETTING DIRECTORIES (from configuration file).

# # Where PRSToolKit resides
# PRSTOOLKITDIR=${PRSTOOLKITDIR} # from configuration file

# # Data information
# LDDATA=${LDDATA} # from configuration file
# VALIDATIONDATA="${VALIDATIONDATA}/${VALIDATIONFILE}" # from configuration file






# 	if [[ ${VALIDATIONQC} == "YES" ]]; then

# 		echo ""
# 		echobold "#========================================================================================================"
# 		echobold "#== OPTIONAL VALIDATION QUALITY CONTROL IS IN EFFECT [DEFAULT]"
# 		echobold "#========================================================================================================"
# 		echobold "#"
# 		echo ""
# 		echo "We will perform quality control on the validation dataset [${VALIDATIONNAME}] which is in [${VALIDATIONFORMAT}]-format."
	
# 		### Example head of STATS-file 15 16 19
# 		### SNPID RSID Chr BP A_allele B_allele MinorAllele MajorAllele AA AB BB AA_calls AB_calls BB_calls MAF HWE missing missing_calls Info CAF
# 		### --- 1:10177:A:AC 01 10177 A AC AC A 554.34 731.78 239.87 46 56 8 0.39696 0.018899 6.3195e-06 0.92792 0.35024 0.396962
# 		### --- 1:10235:T:TA 01 10235 T TA TA T 1524.2 1.8055 0 1523 0 0 0.00059159 4.8216e-17 0 0.0019659 0.26078 0.000591577
# 		### --- rs145072688:10352:T:TA 01 10352 T TA TA T 490.88 755.85 279.26 46 55 15 0.43066 0.14578 1.0239e-05 0.92398 0.34431 0.430661
# 		### --- 1:10505:A:T 01 10505 A T T A 1525.7 0.34198 0 1525 0 0 0.00011205 -0 0 0.00065531 0.2532 0.000112048
# # 		
# # 		echo "SNPID RSID Chr BP alleleA alleleB HWE Info CAF" > ${SUBPROJECTDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.keep.txt
# # 		zcat ${VALIDATIONDATA}/aegs_combo_1kGp3GoNL5_RAW.stats.gz | tail -n +2 | awk ' $15 > '$MAF' && $16 > '$HWE' && $19 > '$INFO' ' >> ${SUBPROJECTDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.keep.txt
# # 		
# # 		cat ${SUBPROJECTDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.keep.txt | awk '{ print $2 }' > ${SUBPROJECTDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.keeptofilter.txt
# # 		
# 		for CHR in $(seq 1 22) X; do 
# 		### FOR DEBUGGING
# 		### for CHR in 22; do
		
# 			echo ""
# 			echo "* processing chromosome ${CHR} and extracting relevant variants."
# 			echo "${QCTOOL} -g ${VALIDATIONDATA}/${VALIDATIONFILE}${CHR}.gen.gz -s ${VALIDATIONDATA}/${VALIDATIONFILE}${CHR}.sample -excl-rsids ${SUBPROJECTDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.keeptofilter.txt -og ${PARSEDDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.chr${CHR}.gen -os ${PARSEDDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.chr${CHR}.sample" > ${PARSEDDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.chr${CHR}.filter.sh
# 			qsub -S /bin/bash -N PRS.FILTER.${VALIDATIONNAME}.chr${CHR} -o ${PARSEDDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.chr${CHR}.filter.log -e ${PARSEDDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.chr${CHR}.filter.errors -l h_vmem=${QMEMFILTER} -l h_rt=${QRUNTIMEFILTER} -wd ${PARSEDDIR} ${PARSEDDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.chr${CHR}.filter.sh
		
# 			echo ""
# 			echo "* converting to PLINK-binary format."
# 			echo "${QCTOOL} -g ${PARSEDDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.chr${CHR}.gen -s ${PARSEDDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.chr${CHR}.sample -threshhold ${THRESHOLD} -og ${PARSEDDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.chr${CHR} -ofiletype binary_ped" > ${PARSEDDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.chr${CHR}.convert.sh
# 			qsub -S /bin/bash -N PRS.CONVERT.${VALIDATIONNAME}.chr${CHR} -hold_jid PRS.FILTER.${VALIDATIONNAME}.chr${CHR} -o ${PARSEDDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.chr${CHR}.convert.log -e ${PARSEDDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.chr${CHR}.convert.errors -l h_vmem=${QMEMCONVERT} -l h_rt=${QRUNTIMECONVERT} -wd ${PARSEDDIR} ${PARSEDDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.chr${CHR}.convert.sh
		
# 			echo ""
# 			echo "* deleting old files."
# 			echo "rm -v ${PARSEDDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.chr${CHR}.gen ${PARSEDDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.chr${CHR}.sample" > ${PARSEDDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.chr${CHR}.clean.sh
# 			qsub -S /bin/bash -N PRS.CLEAN.${VALIDATIONNAME}.chr${CHR} -hold_jid PRS.CONVERT.${VALIDATIONNAME}.chr${CHR} -o ${PARSEDDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.chr${CHR}.clean.log -e ${PARSEDDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.chr${CHR}.clean.errors -l h_vmem=${QMEM} -l h_rt=${QRUNTIME} -wd ${PARSEDDIR} ${PARSEDDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.chr${CHR}.clean.sh
		
# 		done

# 	elif [[ ${VALIDATIONQC} == "NO" ]]; then

# 		echo ""
# 		echobold "#========================================================================================================"
# 		echobold "#== PARSING DATA TO PLINK-BINARY FORMAT"
# 		echobold "#========================================================================================================"
# 		echobold "#"
# 		echo ""
# 		echo "We will parse the validation dataset [${VALIDATIONNAME}], which is in [${VALIDATIONFORMAT}]-format, to PLINK-binary format."
	
# 		for CHR in $(seq 1 22) X; do 
# 		### FOR DEBUGGING
# 		### for CHR in 22; do
		
# 			echo ""
# 			echo "* processing chromosome ${CHR} and converting to PLINK-binary format."
# 			echo "${QCTOOL} -g ${VALIDATIONDATA}/${VALIDATIONFILE}${CHR}.gen.gz -s ${VALIDATIONDATA}/${VALIDATIONFILE}${CHR}.sample -threshhold ${THRESHOLD} -og ${PARSEDDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.chr${CHR} -ofiletype binary_ped" > ${PARSEDDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.chr${CHR}.convert.sh
# 			qsub -S /bin/bash -N PRS.CONVERT.${VALIDATIONNAME}.chr${CHR} -o ${PARSEDDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.chr${CHR}.convert.log -e ${PARSEDDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.chr${CHR}.convert.errors -l h_vmem=${QMEMCONVERT} -l h_rt=${QRUNTIMECONVERT} -wd ${PARSEDDIR} ${PARSEDDIR}/${VALIDATIONNAME}.${POPULATION}.${REFERENCE}.chr${CHR}.convert.sh
		
# 		done
	
# 	fi

# 








# 	echobold "#========================================================================================================"
# 	echobold "#== VALIDATION QUALITY CONTROL"
# 	echobold "#========================================================================================================"
# 	echobold "#"
# 	### REQUIRED: VEGAS/VEGAS2 settings.
# 	### Note: we do `cd ${VEGASDIR}` because VEGAS is making temp-files in a special way, 
# 	###       adding a date-based number in front of the input/output files.
# 	echo "Creating VEGAS input files..." 
# 	mkdir -v ${SUBPROJECTDIR}/vegas
# 	VEGASRESULTDIR=${SUBPROJECTDIR}/vegas
# 	chmod -Rv a+rwx ${VEGASRESULTDIR}
# 	echo "...per chromosome."

# 	while IFS='' read -r GWASCOHORT || [[ -n "$GWASCOHORT" ]]; do
# 			LINE=${GWASCOHORT}
# 			COHORT=$(echo "${LINE}" | awk '{ print $1 }')
# 			echo "     * ${COHORT}"
		
# 			if [ ! -d ${VEGASRESULTDIR}/${COHORT} ]; then
# 				echo "> VEGAS results directory doesn't exist - Mr. Bourne will create it for you."
# 				mkdir -v ${VEGASRESULTDIR}/${COHORT}
# 				chmod -Rv a+rwx ${VEGASRESULTDIR}/${COHORT}
# 			else
# 				echo "> VEGAS results directory already exists."
# 				chmod -Rv a+rwx ${VEGASRESULTDIR}/${COHORT}
# 			fi

# 		for CHR in $(seq 1 23); do
# 			if [[ $CHR -le 22 ]]; then 
# 				echo "Processing chromosome ${CHR}..."
# 				echo "zcat ${PARSEDDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.FINAL.txt.gz | ${SCRIPTS}/parseTable.pl --col ${VARIANTID},CHR,P | awk ' \$2==${CHR} ' | awk '{ print \$1, \$3 }' | tail -n +2 > ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.forVEGAS.txt " > ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.createVEGAS.sh
# # 				qsub -S /bin/bash -N VEGAS2.${PROJECTNAME}.chr${CHR}.create -hold_jid gwas.wrapper -o ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.createVEGAS.log -e ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.createVEGAS.errors -l h_vmem=${QMEMVEGAS} -l h_rt=${QRUNTIMEVEGAS} -wd ${VEGASRESULTDIR} ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.createVEGAS.sh
			
# 				echo "cd ${VEGASRESULTDIR}/${COHORT} " > ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.runVEGAS.sh
# 				echo "$VEGAS2 -G -snpandp ${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.forVEGAS.txt -custom ${VEGAS2POP}.chr${CHR} -glist ${VEGAS2GENELIST} -upper ${VEGAS2UPPER} -lower ${VEGAS2LOWER} -chr ${CHR} -out ${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.fromVEGAS " >> ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.runVEGAS.sh
# # 				qsub -S /bin/bash -N VEGAS2.${PROJECTNAME}.chr${CHR} -hold_jid VEGAS2.${PROJECTNAME}.chr${CHR}.create -o ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.runVEGAS.log -e ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.runVEGAS.errors -l h_vmem=${QMEMVEGAS} -l h_rt=${QRUNTIMEVEGAS} -wd ${VEGASRESULTDIR} ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.runVEGAS.sh
	
# 			elif [[ $CHR -eq 23 ]]; then  
# 				echo "Processing chromosome X..."
# 				echo "zcat ${PARSEDDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.FINAL.txt.gz | ${SCRIPTS}/parseTable.pl --col ${VARIANTID},CHR,P | awk ' \$2==\"X\" ' | awk '{ print \$1, \$3 }' | tail -n +2 > ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.forVEGAS.txt " > ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.createVEGAS.sh
# # 				qsub -S /bin/bash -N VEGAS2.${PROJECTNAME}.chr${CHR}.create -hold_jid gwas.wrapper -o ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.createVEGAS.log -e ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.createVEGAS.errors -l h_vmem=${QMEMVEGAS} -l h_rt=${QRUNTIMEVEGAS} -wd ${VEGASRESULTDIR} ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.createVEGAS.sh
			
# 				echo "cd ${VEGASRESULTDIR}/${COHORT} " > ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.runVEGAS.sh
# 				echo "$VEGAS2 -G -snpandp ${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.forVEGAS.txt -custom ${VEGAS2POP}.chr${CHR} -glist ${VEGAS2GENELIST} -upper ${VEGAS2UPPER} -lower ${VEGAS2LOWER} -chr ${CHR} -out ${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.fromVEGAS " >> ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.runVEGAS.sh
# # 				qsub -S /bin/bash -N VEGAS2.${PROJECTNAME}.chr${CHR} -hold_jid VEGAS2.${PROJECTNAME}.chr${CHR}.create -o ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.runVEGAS.log -e ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.runVEGAS.errors -l h_vmem=${QMEMVEGAS} -l h_rt=${QRUNTIMEVEGAS} -wd ${VEGASRESULTDIR} ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.runVEGAS.sh
	
# 			else
# 				echo "*** ERROR *** Something is rotten in the City of Gotham; most likely a typo. Double back, please."	
# 				exit 1
# 			fi

# 		done
# 	done < ${BASEDATA}
#
	
# script_copyright_message
