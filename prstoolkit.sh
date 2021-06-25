#!/bin/bash
#SBATCH --job-name="PRS workflow"                                 #Job name
#SBATCH --output=/home/dhl_ec/aligterink/logs/batchjob_%j.log     # Standard output and error log
#SBATCH --mail-user=a.j.ligterink@umcutrecht.nl                   # Mail
#SBATCH --mail-type=NONE		                                  # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --time=50:00:00                                           # Time limit hrs:min:sec
#SBATCH --mem=4G	                                              # RAM required per node

##########################################################################################
### Setting colouring
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

# Creating color-functions
function echobold { 
    echo -e "${BOLD}${1}${NONE}"
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
function echocyan {
    echo -e "${CYAN}${1}${NONE}"
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
	echoerror ""
	echoerror "An example command would be: prstoolkit.sh [arg1]"
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
echobold "* Last update:  2021-06-21"
echobold "* Written by:   Sander W. van der Laan | s.w.vanderlaan@gmail.com"
echobold "	        Anton Ligterink"
echobold "* Description:  Prepare files and calculate polygenic scores. This tool will do the following:"
echobold "                - Perform quality control on GWAS summary data"
echobold "                - Calculate polygenic risk scores from GWAS summary statistcs (PRSice/RapidoPGS/PRS-CS)"
echobold " 		- Calculate polygenic scores using premade scoring systems (PLINK)"
echobold "* References:   RapidoPGS: https://rdrr.io/cran/RapidoPGS/man/rapidopgs_single.html"
echobold "                LDpred: https://github.com/bvilhjal/ldpred      PRSICE: https://www.prsice.info/"
echobold "                PRScs: https://github.com/getian107/PRScs       PLINK: http://zzz.bwh.harvard.edu/plink/"
echobold "* REQUIRED: "
echobold "  - A high-performance computer cluster with a SLURM system"
echobold "  - R v3.6+, Python 3.7+"
echobold "  - Required R 3.6+ modules: [RapidoPGS], [data.table], [optparse], [remotes], [method], [tools], [ggplot2], [grDevices], [RColorBrewer]"
echobold "  - Required Python 3.7+ modules: [gzip], [argparse], [scipy], [h5py]"
echobold ""
echobold "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Today's: "$(date)
TODAY=$(date +"%Y%m%d")

##########################################################################################
### SET THE SCENE FOR THE SCRIPT
##########################################################################################

##########################################################################################
### Command-line argument check
if [[ $# -lt 1 ]]; then 
	echo ""
	echoerror "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	echoerrorflash "               *** Oh, oh, computer says no! Number of arguments found "$#". ***"
	echoerror "You must supply [1] argument when running *** PRSToolKit ***!"
	script_arguments_error
else
	echo "These are the "$#" arguments that passed:"
	echo "The configuration file..................: "$(basename ${1}) # argument 1
fi 

##########################################################################################
### Loading command-line arguments and configuration file
source "$1" # Depends on arg1.
CONFIGURATIONFILE="$1" # Depends on arg1 -- but also on where it resides!!!

##########################################################################################
### Set the output file name
OUTPUTNAME=${PROJECTNAME}_${PRSMETHOD}_$(date +%Y-%b-%d--%H-%M)_job_${SLURM_JOB_ID}

##########################################################################################
### Parameter validation
if [[ ${PRSMETHOD} == "PLINK" || ${PRSMETHOD} == "PRSCS" || ${PRSMETHOD} == "PRSICE" || ${PRSMETHOD} == "RAPIDOPGS" || ${PRSMETHOD} == "NONE" ]]; then
	echo ""
	echo "Calculating Polygenic Risk Scores using ${PRSMETHOD}."

elif [[ ${PRSMETHOD} == "MANUAL" || ${PRSMETHOD} == "LDPRED" ]]; then
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
	echoerror "PLINK	  -- use PLINK to calculate PGS from posterior SNP effect sizes"
	echoerror "PRSICE	  -- use PRSice to calculate PGS, GWAS summary statistics will be LD-pruned using P-value thresholds"
	echoerror "PRSCS	  -- use PRS-CS to calculate PGS by inferring posterior SNP effect sizes under continuous shrinkage priors"
	echoerror "RAPIDOPGS -- use RapidoPGS to calculate PGS from GWAS summary statistics without requiring an external LD panel"
	echonooption "LDPRED    -- use LDpred to LD-correct GWAS summary statistics and calculate polygenic scores *not implemented yet*"
	echonooption "MANUAL    -- use an in-house developed Rscript to calculate a polygenic score using a limited set of variants *not implemented yet*"
	echoerror "NONE	  -- don't perform PRS, instead only peform quality control"
	echoerror "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	### The wrong arguments are passed, so we'll exit the script now!
	echo ""
	script_copyright_message
	exit 1
fi

##########################################################################################
### Setting up necessary directories
echo ""
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Checking for the existence of the output directory [ ${OUTPUT_DIRNAME} ]."
if [ ! -d ${PROJECT_DIR}/${OUTPUT_DIRNAME} ]; then
	echo "> Output directory doesn't exist - Mr. Bourne will create it for you."
	mkdir -v ${PROJECT_DIR}/${OUTPUT_DIRNAME}
else
	echo "> Output directory already exists."
fi
OUTPUTDIR=${PROJECT_DIR}/${OUTPUT_DIRNAME}

echo ""
echo "Checking for the existence of the subproject directory [ ${SUBPROJECT_DIR_NAME} ]."
if [ ! -d ${OUTPUTDIR}/${SUBPROJECT_DIR_NAME} ]; then
	echo "> Subproject directory doesn't exist - Mr. Bourne will create it for you."
	mkdir -v ${OUTPUTDIR}/${SUBPROJECT_DIR_NAME}
else
	echo "> Subproject directory already exists."
fi
SUBPROJECT_DIR=${OUTPUTDIR}/${SUBPROJECT_DIR_NAME}

echo ""
echo "Checking for the existence of the ${LOG_DIRNAME} directory [ ${PROJECT_DIR}/${LOG_DIRNAME} ]."
if [ ! -d ${PROJECT_DIR}/${LOG_DIRNAME} ]; then
	echo "> ${LOG_DIRNAME} directory doesn't exist - Mr. Bourne will create it for you."
	mkdir -v ${PROJECT_DIR}/${LOG_DIRNAME}
else
	echo "> ${LOG_DIRNAME} directory already exists."
fi
LOGDIR=${PROJECT_DIR}/${LOG_DIRNAME}

echo ""
echo "Checking for the existence of the main working directory [ ${PROJECT_DIR}/${MAIN_WORKDIR_NAME} ]."
if [ ! -d ${PROJECT_DIR}/${MAIN_WORKDIR_NAME} ]; then
	echo "> ${MAIN_WORKDIR_NAME} directory doesn't exist - Mr. Bourne will create it for you."
	mkdir -v ${PROJECT_DIR}/${MAIN_WORKDIR_NAME}
else
	echo "> ${MAIN_WORKDIR_NAME} directory already exists."
fi
MAIN_WORKDIR=${PROJECT_DIR}/${MAIN_WORKDIR_NAME}

echo ""
echo "Checking for the existence of the working directory for this project [ ${MAIN_WORKDIR}/${OUTPUTNAME} ]."
if [ ! -d ${MAIN_WORKDIR}/${OUTPUTNAME} ]; then
	echo "> ${OUTPUTNAME} doesn't exist - Mr. Bourne will create it for you."
	mkdir -v ${MAIN_WORKDIR}/${OUTPUTNAME}
else
	echo "> ${OUTPUTNAME} already exists."
fi
PRSDIR=${MAIN_WORKDIR}/${OUTPUTNAME}

echo ""
echo "The scene is properly set, and directories are created! 🖖"
echo "PRSToolKit directory............................................: "${PRSTOOLKITDIR}
echo "LD reference data directory.....................................: "${LDDATA}
echo "Validation data directory.......................................: "${VALIDATIONDATA}
echo "Main directory..................................................: "${PROJECT_DIR}
echo "Main analysis output directory..................................: "${OUTPUTDIR}
echo "Subproject's analysis output directory..........................: "${SUBPROJECT_DIR}
echo "Log directory...................................................: "${LOGDIR}
echo "Main working directory..........................................: "${MAIN_WORKDIR}
echo "Working directory for this project..............................: "${PRSDIR}
echo ""
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo ""

##########################################################################################
### Save the configuration file to the results folder depending on the SAVE_CONFIG parameter
if [[ ${SAVE_CONFIG} == "TRUE" ]]; then	
	echo "SAVE_CONFIG parameter is active, the configuration file will be saved in [ ${SUBPROJECT_DIR}/${OUTPUTNAME}.config ]."
	cp ${CONFIGURATIONFILE} ${SUBPROJECT_DIR}/${OUTPUTNAME}.config

elif [[ ${SAVE_CONFIG} == "FALSE" || ${SAVE_CONFIG} == "" ]]; then
	true

else
	echo "SAVE_CONFIG parameter not recognized, the configuration file will be saved in [ ${SUBPROJECT_DIR}/${OUTPUTNAME}.config ]."
	cp ${CONFIGURATIONFILE} ${SUBPROJECT_DIR}/${OUTPUTNAME}.config

fi
echo ""
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

##########################################################################################
### Perform quality control if enabled
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
	QC_SUMMARY_FILE=${PRSDIR}/${OUTPUTNAME}_QC_results.txt
	QC_OUTPUT=${PRSDIR}/QCd_basefile.txt.gz
	QC_JOBID=$(sbatch --parsable --wait --job-name=PRS_QC --time ${RUNTIME_QC} --mem ${MEMORY_QC} -o ${LOGDIR}/${OUTPUTNAME}_QC.log ${PRSTOOLKITSCRIPTS}/QC.py -b ${BASEDATA} -i ${BF_ID_COL} -s ${STATS_FILE} -a ${STATS_ID_COL} -c ${STATS_MAF_COL} -d ${STATS_INFO_COL} -m ${MAF} -n ${INFO} -o ${QC_OUTPUT} -r ${QC_SUMMARY_FILE})
	QC_DEPENDENCY="--dependency=afterok:${QC_JOBID}"

	echo ""
	cat ${QC_SUMMARY_FILE}
	echo ""

	# We will now use the quality controlled file as our new base file
	BASEDATA=${QC_OUTPUT}

elif [[ ${QC} == "NO" ]]; then
	echo ""
	echo "Quality control will not be performed"
	echo ""
else 
	echo ""
	echo "QC parameter should be [YES/NO], not \"${QC}\". Exiting..."
	echo ""
	exit 1
fi

echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo ""
echo "Risk score computation will now commence."
echo ""

# Polygenic scores will be stored in this file
RESULTS_FILE=${SUBPROJECT_DIR}/${OUTPUTNAME}_results.txt

##########################################################################################
### Run the individual PRS method scripts
if [[ ${PRSMETHOD} == "PLINK" ]]; then

	# Calculate scores per chromosome for each individual in the target population
	### TODO: make this run in parallel
	PLINK_HEADER="TRUE"
	PLINKSCORE_JOBID=$(sbatch --parsable --wait --job-name=PRS_PLINKSCORE ${QC_DEPENDENCY} --time ${RUNTIME_PLINKSCORE} --mem ${MEMORY_PLINKSCORE} -o ${LOGDIR}/${OUTPUTNAME}_PLINK_score.log --export=ALL,VALIDATIONDATA=${VALIDATIONDATA},VALIDATIONPREFIX=${VALIDATIONPREFIX},PLINK=${PLINK},REF_POS=${VAL_REF_POS},SAMPLE_FILE=${SAMPLE_FILE},PRSDIR=${PRSDIR},WEIGHTS_FILE=${BASEDATA},SNP_COL=${BF_ID_COL},EFFECT_COL=${BF_EFFECT_COL},SCORE_COL=${BF_STAT_COL},PLINK_SETTINGS=${PLINK_SETTINGS},PLINK_HEADER=${PLINK_HEADER} ${PRSTOOLKITSCRIPTS}/plinkscore.sh)
	PLINKSCORE_DEPENDENCY="--dependency=afterok:${PLINKSCORE_JOBID}"

	# Sum the scores to calculate the final polygenic score
	PLINKSUM_JOBID=$(sbatch --parsable --wait --job-name=PRS_PLINKSUM ${PLINKSCORE_DEPENDENCY} --time ${RUNTIME_PLINKSUM} --mem ${MEMORY_PLINKSUM} -o ${LOGDIR}/${OUTPUTNAME}_PLINK_sum.log ${PRSTOOLKITSCRIPTS}/sum_plink_scores.R -s ${SAMPLE_FILE} -f 1 -i 2 -d ${PRSDIR} -p plink2_${VALIDATIONPREFIX} -r "SCORE1_SUM" -o ${RESULTS_FILE})

elif [[ ${PRSMETHOD} == "PRSCS" ]]; then

	# Parse the base file to PRS-CS format (and also unzip it)
	PARSED_BASEDATA=${PRSDIR}/basefile_PRScs_format.txt
	PRSCS_format_JOBID=$(sbatch --parsable --wait --job-name=PRS_PRScs_format ${QC_DEPENDENCY} --time ${RUNTIME_PRSCS_format} --mem ${MEMORY_PRSCS_format} -o ${LOGDIR}/${OUTPUTNAME}_PRScs_format.log ${PRSTOOLKITSCRIPTS}/basefile_PRScs_formatter.R -i ${BASEDATA} -o ${PARSED_BASEDATA} -d ${BF_ID_COL} -r ${BF_EFFECT_COL} -a ${BF_NON_EFFECT_COL} -z ${BF_STAT} -m ${BF_STAT_COL} -v ${BF_PVALUE_COL})
	PRSCS_format_DEPENDENCY="--dependency=afterok:${PRSCS_format_JOBID}"
	WEIGHTS_FILE=${PRSDIR}/PRScs_weights_combined.txt

	# Run PRS-CS
	PRSCS_JOBID=$(sbatch --parsable --wait --job-name=PRS_PRScs ${PRSCS_format_DEPENDENCY} --time ${RUNTIME_PRSCS} --mem ${MEMORY_PRSCS} -c ${PRSCS_CPUS} -o ${LOGDIR}/${OUTPUTNAME}_PRScs.log --export=ALL,VALIDATIONDATA=${VALIDATIONDATA},VALIDATIONPREFIX=${VALIDATIONPREFIX},SAMPLE_FILE=${SAMPLE_FILE},BIM_FILE_AVAILABLE=${BIM_FILE_AVAILABLE},BIM_FILE_PATH=${BIM_FILE_PATH},PRSDIR=${PRSDIR},REF_POS=${VAL_REF_POS},BF_SAMPLE_SIZE=${BF_SAMPLE_SIZE},PLINK=${PLINK},PYTHONPATH=${PYTHONPATH},LDDATA=${LDDATA},PRSCS=${PRSCS},PARSED_BASEDATA=${PARSED_BASEDATA},WEIGHTS_FILE=${WEIGHTS_FILE},PRSCS_THREADS=${PRSCS_THREADS} ${PRSTOOLKITSCRIPTS}/prscs.sh)
	PRSCS_DEPENDENCY="--dependency=afterok:${PRSCS_JOBID}"

	# Calculate individual scores using PLINK
	PLINK_HEADER="FALSE"
	PLINKSCORE_JOBID=$(sbatch --parsable --wait --job-name=PRS_PLINKSCORE ${PRSCS_DEPENDENCY} --time ${RUNTIME_PLINKSCORE} --mem ${MEMORY_PLINKSCORE} -o ${LOGDIR}/${OUTPUTNAME}_PLINK_score.log --export=ALL,VALIDATIONDATA=${VALIDATIONDATA},VALIDATIONPREFIX=${VALIDATIONPREFIX},PLINK=${PLINK},REF_POS=${VAL_REF_POS},SAMPLE_FILE=${SAMPLE_FILE},PRSDIR=${PRSDIR},WEIGHTS_FILE=${WEIGHTS_FILE},SNP_COL=2,EFFECT_COL=4,SCORE_COL=6,PLINK_SETTINGS=${PLINK_SETTINGS},PLINK_HEADER=${PLINK_HEADER} ${PRSTOOLKITSCRIPTS}/plinkscore.sh)
	PLINKSCORE_DEPENDENCY="--dependency=afterok:${PLINKSCORE_JOBID}"

	# # Sum the effect sizes to calculate the final score
	PLINKSUM_JOBID=$(sbatch --parsable --wait --job-name=PRS_PLINKSUM ${PLINKSCORE_DEPENDENCY} --time ${RUNTIME_PLINKSUM} --mem ${MEMORY_PLINKSUM} -o ${LOGDIR}/${OUTPUTNAME}_PLINK_sum.log ${PRSTOOLKITSCRIPTS}/sum_plink_scores.R -s ${SAMPLE_FILE} -f 1 -i 2 -d ${PRSDIR} -p plink2_${VALIDATIONPREFIX} -r "SCORE1_SUM" -o ${RESULTS_FILE})

elif [[ ${PRSMETHOD} == "LDPRED" ]]; then
	# cd ${PRSDIR}
	# ldpred_cpus=12
	# LDPRED_JOBID=$(sbatch --parsable --wait --job-name=PRS_LDPRED ${QC_DEPENDENCY} --time ${RUNTIME_LDPRED} --mem ${MEMORY_LDPRED} -c ${ldpred_cpus} -o ${LOGDIR}/${OUTPUTNAME}_LDPRED.log ${PRSTOOLKITSCRIPTS}/LDpred2.R)
	echo "not implemented"

elif [[ ${PRSMETHOD} == "PRSICE" ]]; then

	# Run PRSice
	PRSICE_OUTPUTNAME=${PRSDIR}/out
	PRSICE_JOBID=$(sbatch --parsable --wait --job-name=PRS_PRSICE ${QC_DEPENDENCY} --time ${RUNTIME_PRSICE} --mem ${MEMORY_PRSICE} -c ${PRSICE_CPUS} -o ${LOGDIR}/${OUTPUTNAME}_PRSICE.log --export=ALL,RSCRIPT=${RSCRIPT},PRSICE2_R=${PRSICE2_R},PRSICE2_SH=${PRSICE2_SH},PRSDIR=${PRSDIR},BASEDATA=${BASEDATA},TARGETDATA=${VALIDATIONDATA}/${VALIDATIONPREFIX},BF_STAT=${BF_STAT},PRSICE_PHENOTYPE_BINARY=${PRSICE_PHENOTYPE_BINARY},BF_ID_COL=${BF_ID_COL},BF_CHR_COL=${BF_CHR_COL},BF_POS_COL=${BF_POS_COL},BF_EFFECT_COL=${BF_EFFECT_COL},BF_NON_EFFECT_COL=${BF_NON_EFFECT_COL},BF_STAT_COL=${BF_STAT_COL},BF_PVALUE_COL=${BF_PVALUE_COL},PHENOTYPEFILE=${SAMPLE_FILE},PRSICE_PHENOTYPE=${PRSICE_PHENOTYPE},PRSICE_CLUMP_KB=${PRSICE_CLUMP_KB},PRSICE_CLUMP_P=${PRSICE_CLUMP_P},PRSICE_CLUMP_R2=${PRSICE_CLUMP_R2},PRSICE_PERM=${PRSICE_PERM},PRSICE_THREADS=${PRSICE_THREADS},PRSICE_SETTINGS="${PRSICE_SETTINGS}",PRSICE_OUTPUTNAME=${PRSICE_OUTPUTNAME},LDDATA="${LDDATA}",PRSICE_EXTRACT="${PRSICE_EXTRACT}",PRSICE_EXCLUDE="${PRSICE_EXCLUDE}" ${PRSTOOLKITSCRIPTS}/prsice.sh)
	
	# Copy the results from the output to the results folder
	awk '//{print $1,$2,$4 }' ${PRSICE_OUTPUTNAME}.best > ${RESULTS_FILE}
	echo "The scores best fitted to the \"${PRSICE_PHENOTYPE}\" phenotype were written to the results folder."
	echo "Note that PRSice has potentially stored additional scores at different P-value thresholds in the work directory depending on the provided parameters."
	echo "For more info visit https://www.prsice.info/step_by_step/"

elif [[ ${PRSMETHOD} == "RAPIDOPGS" ]]; then

	# File for storing the effect sizes computed by Rapido
	WEIGHTS_FILE="${PRSDIR}/Rapido_weights.txt"

	# Calculate weights
	RAPIDO_JOBID=$(sbatch --parsable --wait --job-name=PRS_RAPIDO ${QC_DEPENDENCY} --time ${RUNTIME_RAPIDO} --mem ${MEMORY_RAPIDO} -o ${LOGDIR}/${OUTPUTNAME}_RAPIDO.log ${PRSTOOLKITSCRIPTS}/rapidopgs.R -k ${PRSDIR} -o ${WEIGHTS_FILE} -b ${BASEDATA} -d ${BF_BUILD} -i ${BF_ID_COL} -c ${BF_CHR_COL} -p ${BF_POS_COL} -r ${BF_NON_EFFECT_COL} -a ${BF_EFFECT_COL} -f "${BF_FRQ_COL}" -m ${BF_STAT_COL} -e ${BF_SE_COL} -s ${BF_SAMPLE_SIZE} -n "${BF_SBJ_COL}" -g "${RP_filt_threshold}" -j "${RP_recalc}" -l ${BF_TARGET_TYPE} -v "${RP_ppi}" -x "${RP_prior}" -z "${RP_REF}")
	RAPIDO_DEPENDENCY="--dependency=afterok:${RAPIDO_JOBID}"
	
	# Calculate individual scores using PLINK
	PLINK_HEADER="TRUE"
	PLINKSCORE_JOBID=$(sbatch --parsable --wait --job-name=PRS_PLINKSCORE ${RAPIDO_DEPENDENCY} --time ${RUNTIME_PLINKSCORE} --mem ${MEMORY_PLINKSCORE} -o ${LOGDIR}/${OUTPUTNAME}_PLINK_score.log --export=ALL,VALIDATIONDATA=${VALIDATIONDATA},VALIDATIONPREFIX=${VALIDATIONPREFIX},PLINK=${PLINK},REF_POS=${VAL_REF_POS},SAMPLE_FILE=${SAMPLE_FILE},PRSDIR=${PRSDIR},WEIGHTS_FILE=${WEIGHTS_FILE},SNP_COL=SNPID,EFFECT_COL=REF,SCORE_COL=WEIGHT,PLINK_SETTINGS=${PLINK_SETTINGS},PLINK_HEADER=${PLINK_HEADER} ${PRSTOOLKITSCRIPTS}/plinkscore.sh)
	PLINKSCORE_DEPENDENCY="--dependency=afterok:${PLINKSCORE_JOBID}"

	# Sum the effect sizes to calculate the final score
	PLINKSUM_JOBID=$(sbatch --parsable --wait --job-name=PRS_PLINKSUM ${PLINKSCORE_DEPENDENCY} --time ${RUNTIME_PLINKSUM} --mem ${MEMORY_PLINKSUM} -o ${LOGDIR}/${OUTPUTNAME}_PLINK_sum.log ${PRSTOOLKITSCRIPTS}/sum_plink_scores.R -s ${SAMPLE_FILE} -f 1 -i 2 -d ${PRSDIR} -p plink2_${VALIDATIONPREFIX} -r "SCORE1_SUM" -o ${RESULTS_FILE})
fi

echo ""
echo "${PRSMETHOD} risk score calculation has finished."
echo ""

##########################################################################################
### Remove intermediate files depending on the KEEP_TEMP_FILES parameter
if [[ ${KEEP_TEMP_FILES} == "TRUE" ]]; then
	echo "KEEP_TEMP_FILES parameter is active, temporary files stored in [ ${PRSDIR} ] will not be removed."

elif [[ ${KEEP_TEMP_FILES} == "FALSE" || ${KEEP_TEMP_FILES} == "" ]]; then
	echo "KEEP_TEMP_FILES parameter is inactive, temporary files stored in [ ${PRSDIR} ] will now be removed."
	rm -r ${PRSDIR}

else
	echo "KEEP_TEMP_FILES parameter not recognized, temporary files stored in [ ${PRSDIR} ] will not be removed."

fi

echo ""
echocyan "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echobold "Wow. I'm all done buddy. What a job 😱 ! let's have a 🍻🍻 ... 🖖 "

script_copyright_message
