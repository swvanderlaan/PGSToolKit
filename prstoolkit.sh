#!/bin/bash
#
#$ -S /bin/bash
#$ -N calculate_prs_ldpred
# -hold_jid some_other_basic_bash_script
#$ -o /hpc/dhl_ec/svanderlaan/projects/polygenicscores/calculate_prs_ldpred.log
#$ -e /hpc/dhl_ec/svanderlaan/projects/polygenicscores/calculate_prs_ldpred.errors
#$ -l h_rt=02:00:00
#$ -l h_vmem=8G
# -l tmpspace=64G
#$ -M s.w.vanderlaan-2@umcutrecht.nl
#$ -m ea
#$ -cwd 

### BASH OPTIONS
# -S 				# the type of BASH you'd like to use
# -N 				# the name of this script
# -hold_jid 		# the current script (basic_bash_script) will hold until some_other_basic_bash_script has finished
# -o 				# the log file of this job
# -e 				# the error file of this job
# -l h_rt=02:00:00  # h_rt=[max time, e.g. 02:02:01] - this is the time you think the script will take
# -l h_vmem=8G  	# h_vmem=[max. mem, e.g. 45G] - this is the amount of memory you think your script will use
# -l tmpspace=64G  	# this is the amount of temporary space you think your script will use
# -M 				# you can send yourself emails when the job is done; "-M" and "-m" go hand in hand
# -m 	  			# you can choose: b=begin of job; e=end of job; a=abort of job; s=suspended job; n=no mail is send
# -cwd  			# set the job start to the current directory - so all the things in this script are relative to the current directory!!!

### Creating display functions
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
# errors no option
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
	echoerror "- Argument #2 is path_to/filename of the list of GWAS files with names."
	echoerror ""
	echoerror "An example command would be: prstoolkit.sh [arg1] [arg2]"
	echoerror "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
 	echo ""
	script_copyright_message
	exit 1
}

echobold "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echobold "                        PRSToolKit: A TOOLKIT TO CALCULATE POLYGENIC SCORES"
echoitalic "                    --- prepare and calculate polygenic scores using LDpred ---"
echobold ""
echobold "* Version:      v1.0.0"
echobold ""
echobold "* Last update:  2018-04-05"
echobold "* Written by:   Sander W. van der Laan | s.w.vanderlaan@gmail.com."
echobold "* Description:  Prepare files and calculate polygenic scores using LDpred. It will do the following:"
echobold "                - Perform QC prior of validation data to parsing to required format."
echobold "                - Automatically parse the imputed data (Oxford-format) to hard-coded PLINK-style format."
echobold "                - Perform LD-correct GWAS summary statistics (LDpred)."
echobold "                - Calculate polygenic scores (LDpred)."
echobold ""
echobold "* Reference   : https://github.com/bvilhjal/ldpred."
echobold "* REQUIRED: "
echobold "  - A high-performance computer cluster with a qsub system"
echobold "  - R v3.2+, Python 2.7+, Perl."
echobold "  - Required Python 2.7+ modules: [pandas], [scipy], [numpy], [plinkio], [h5py]."
echobold "  - Note: it will also work on a Mac OS X system with R and Python installed."
### ADD-IN: function to check requirements...
### This might be a viable option! https://gist.github.com/JamieMason/4761049
echobold ""
echobold "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Today's: "$(date)
TODAY=$(date +"%Y%m%d")

##########################################################################################
### SET THE SCENE FOR THE SCRIPT
##########################################################################################

### START of if-else statement for the number of command-line arguments passed ###
if [[ $# -lt 2 ]]; then 
	echo ""
	echoerror "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	echoerrorflash "               *** Oh, oh, computer says no! Number of arguments found "$#". ***"
	echoerror "You must supply [2] arguments when running *** PRSToolKit ***!"
	script_arguments_error
else
	echo "These are the "$#" arguments that passed:"
	echo "The configuration file.................: "$(basename ${1}) # argument 1
	
	### SETTING DIRECTORIES (from configuration file).
	# Loading the configuration file (please refer to the GBASToolKit-Manual for specifications of this file). 
	source "$1" # Depends on arg1.
	
	CONFIGURATIONFILE="$1" # Depends on arg1 -- but also on where it resides!!!
	GWASFILES="$2" # Depends on arg2 -- all the GWAS dataset information
	
	# Where GBASToolKit resides
	PRSTOOLKITDIR=${PRSTOOLKITDIR} # from configuration file
	
	# Data information
	LDDATA=${LDDATA} # from configuration file
	VALIDATIONDATA=${VALIDATIONDATA} # from configuration file

	##########################################################################################
	### CREATE THE OUTPUT DIRECTORIES
	echo ""
	echo "Checking for the existence of the output directory [ ${OUTPUTDIRNAME} ]."
	if [ ! -d ${PROJECTDIR}/${OUTPUTDIRNAME} ]; then
		echo "> Output directory doesn't exist - Mr. Bourne will create it for you."
		mkdir -v ${PROJECTDIR}/${OUTPUTDIRNAME}
	else
		echo "> Output directory already exists."
	fi
	OUTPUTDIR=${OUTPUTDIRNAME}
	
	echo ""
	echo "Checking for the existence of the subproject directory [ ${OUTPUTDIR}/${SUBPROJECTDIRNAME} ]."
	if [ ! -d ${PROJECTDIR}/${OUTPUTDIR}/${SUBPROJECTDIRNAME} ]; then
		echo "> Subproject directory doesn't exist - Mr. Bourne will create it for you."
		mkdir -v ${PROJECTDIR}/${OUTPUTDIR}/${SUBPROJECTDIRNAME}
	else
		echo "> Subproject directory already exists."
	fi
	SUBPROJECTDIR=${PROJECTDIR}/${OUTPUTDIR}/${SUBPROJECTDIRNAME}
	
	echo ""
	echo "Checking for the existence of the parsed data directory [ ${OUTPUTDIR}/${SUBPROJECTDIRNAME}/PARSED ]."
	if [ ! -d ${PROJECTDIR}/${OUTPUTDIR}/${SUBPROJECTDIRNAME}/PARSED ]; then
		echo "> Parsed data directory doesn't exist - Mr. Bourne will create it for you."
		mkdir -v ${PROJECTDIR}/${OUTPUTDIR}/${SUBPROJECTDIRNAME}/PARSED
	else
		echo "> Parsed data directory already exists."
	fi
	PARSEDDIR=${PROJECTDIR}/${OUTPUTDIR}/${SUBPROJECTDIRNAME}/PARSED

	##########################################################################################
	### SETTING UP NECESSARY DIRECTORIES
	echo ""
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	echo ""
	echo "The scene is properly set, and directories are created! 🖖"
	echo "PRSToolKit directory...........................: "${PRSTOOLKITDIR}
	echo "LD reference data directory....................: "${LDDATA}
	echo "Validation data directory......................: "${VALIDATIONDATA}
	echo "Main directory.................................: "${PROJECTDIR}
	echo "Main analysis output directory.................: "${OUTPUTDIR}
	echo "Subproject's analysis output directory.........: "${SUBPROJECTDIR}
	echo "Parsed data is stored in.......................: "${PARSEDDIR}
	echo "We are processing these cohort(s)..............:"
	while IFS='' read -r GWASCOHORT || [[ -n "$GWASCOHORT" ]]; do
		LINE=${GWASCOHORT}
		COHORT=$(echo "${LINE}" | awk '{ print $1 }')
		echo " * ${COHORT}"
	done < ${GWASFILES}
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	echo ""
	
	### Making raw data directories, unless they already exist. Depends on arg2.
	if [[ ${PRSMETHOD} = "LDPRED" ]]; then
	
		echo "Calculating Polygenic Risk Scores using LDpred (by Vilhjálmsson et al. AJHG 2016)."

	
	elif [[ ${PRSMETHOD} = "MANUAL" || ${PRSMETHOD} = "PLINK" || ${PRSMETHOD} = "PRSICE" ]]; then
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
	  	echoerror "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	  	echoerror ""
	  	echoerrorflash "                  *** Oh, computer says no! Argument not recognised. ***"
	  	echoerror "You have the following options as reference for the quality control"
	  	echoerror "and meta-analysis:"
	  	echoerror "LDPRED   -- use LDpred to LD-correct GWAS summary statistics and calculate polygenic scores [default]"
	  	echonooption "MANUAL   -- use an in-house developed Rscript to calculate a polygenic score using a limited set of variants"
	  	echonooption "PLINK    -- use PLINK to calculate scores, GWAS summary statistics will be LD-pruned using p-value thresholds (traditional approach, slow)"
	  	echonooption "PRSICE   -- use PRSice to calculate scores, GWAS summary statistics will be LD-pruned using p-value thresholds (new approach, fast)"
	  	echonooption "(Opaque: *not implemented yet*)"
	  	echoerror "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		### The wrong arguments are passed, so we'll exit the script now!
		echo ""
		script_copyright_message
		exit 1
	fi
# 	
# 	echobold "#========================================================================================================"
# 	echobold "#== GENE-BASED ANALYSIS OF META-ANALYSIS RESULTS USING VEGAS2"
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
# 	
# 	while IFS='' read -r GWASCOHORT || [[ -n "$GWASCOHORT" ]]; do
# 			LINE=${GWASCOHORT}
# 			COHORT=$(echo "${LINE}" | awk '{ print $1 }')
# 			echo "     * ${COHORT}"
# 			
# 			if [ ! -d ${VEGASRESULTDIR}/${COHORT} ]; then
# 				echo "> VEGAS results directory doesn't exist - Mr. Bourne will create it for you."
# 				mkdir -v ${VEGASRESULTDIR}/${COHORT}
# 				chmod -Rv a+rwx ${VEGASRESULTDIR}/${COHORT}
# 			else
# 				echo "> VEGAS results directory already exists."
# 				chmod -Rv a+rwx ${VEGASRESULTDIR}/${COHORT}
# 			fi
# 	
# 		for CHR in $(seq 1 23); do
# 			if [[ $CHR -le 22 ]]; then 
# 				echo "Processing chromosome ${CHR}..."
# 				echo "zcat ${PARSEDDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.FINAL.txt.gz | ${SCRIPTS}/parseTable.pl --col ${VARIANTID},CHR,P | awk ' \$2==${CHR} ' | awk '{ print \$1, \$3 }' | tail -n +2 > ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.forVEGAS.txt " > ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.createVEGAS.sh
# # 				qsub -S /bin/bash -N VEGAS2.${PROJECTNAME}.chr${CHR}.create -hold_jid gwas.wrapper -o ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.createVEGAS.log -e ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.createVEGAS.errors -l h_vmem=${QMEMVEGAS} -l h_rt=${QRUNTIMEVEGAS} -wd ${VEGASRESULTDIR} ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.createVEGAS.sh
# 				
# 				echo "cd ${VEGASRESULTDIR}/${COHORT} " > ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.runVEGAS.sh
# 				echo "$VEGAS2 -G -snpandp ${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.forVEGAS.txt -custom ${VEGAS2POP}.chr${CHR} -glist ${VEGAS2GENELIST} -upper ${VEGAS2UPPER} -lower ${VEGAS2LOWER} -chr ${CHR} -out ${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.fromVEGAS " >> ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.runVEGAS.sh
# # 				qsub -S /bin/bash -N VEGAS2.${PROJECTNAME}.chr${CHR} -hold_jid VEGAS2.${PROJECTNAME}.chr${CHR}.create -o ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.runVEGAS.log -e ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.runVEGAS.errors -l h_vmem=${QMEMVEGAS} -l h_rt=${QRUNTIMEVEGAS} -wd ${VEGASRESULTDIR} ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.runVEGAS.sh
# 		
# 			elif [[ $CHR -eq 23 ]]; then  
# 				echo "Processing chromosome X..."
# 				echo "zcat ${PARSEDDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.FINAL.txt.gz | ${SCRIPTS}/parseTable.pl --col ${VARIANTID},CHR,P | awk ' \$2==\"X\" ' | awk '{ print \$1, \$3 }' | tail -n +2 > ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.forVEGAS.txt " > ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.createVEGAS.sh
# # 				qsub -S /bin/bash -N VEGAS2.${PROJECTNAME}.chr${CHR}.create -hold_jid gwas.wrapper -o ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.createVEGAS.log -e ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.createVEGAS.errors -l h_vmem=${QMEMVEGAS} -l h_rt=${QRUNTIMEVEGAS} -wd ${VEGASRESULTDIR} ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.createVEGAS.sh
# 				
# 				echo "cd ${VEGASRESULTDIR}/${COHORT} " > ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.runVEGAS.sh
# 				echo "$VEGAS2 -G -snpandp ${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.forVEGAS.txt -custom ${VEGAS2POP}.chr${CHR} -glist ${VEGAS2GENELIST} -upper ${VEGAS2UPPER} -lower ${VEGAS2LOWER} -chr ${CHR} -out ${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.fromVEGAS " >> ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.runVEGAS.sh
# # 				qsub -S /bin/bash -N VEGAS2.${PROJECTNAME}.chr${CHR} -hold_jid VEGAS2.${PROJECTNAME}.chr${CHR}.create -o ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.runVEGAS.log -e ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.runVEGAS.errors -l h_vmem=${QMEMVEGAS} -l h_rt=${QRUNTIMEVEGAS} -wd ${VEGASRESULTDIR} ${VEGASRESULTDIR}/${COHORT}/${COHORT}.${PROJECTNAME}.${REFERENCE}.${POPULATION}.chr${CHR}.runVEGAS.sh
# 		
# 			else
# 				echo "*** ERROR *** Something is rotten in the City of Gotham; most likely a typo. Double back, please."	
# 				exit 1
# 			fi
# 
# 		done
# 	done < ${GWASFILES}
# 
# 	echobold "#========================================================================================================"
# 	echobold "#== GENE-BASED ANALYSIS OF META-ANALYSIS RESULTS USING MAGMA"
# 	echobold "#========================================================================================================"
# 	echobold "#"
	
	### END of if-else statement for the number of command-line arguments passed ###
fi 

script_copyright_message

###	UtrechtSciencePark Colours Scheme
###
### Website to convert HEX to RGB: http://hex.colorrrs.com.
### For some functions you should divide these numbers by 255.
###
###	No.	Color				HEX		RGB							CMYK					CHR		MAF/INFO
### --------------------------------------------------------------------------------------------------------------------
###	1	yellow				#FBB820 (251,184,32)				(0,26.69,87.25,1.57) 	=>	1 		or 1.0 > INFO
###	2	gold				#F59D10 (245,157,16)				(0,35.92,93.47,3.92) 	=>	2		
###	3	salmon				#E55738 (229,87,56) 				(0,62.01,75.55,10.2) 	=>	3 		or 0.05 < MAF < 0.2 or 0.4 < INFO < 0.6
###	4	darkpink			#DB003F ((219,0,63)					(0,100,71.23,14.12) 	=>	4		
###	5	lightpink			#E35493 (227,84,147)				(0,63,35.24,10.98) 	=>	5 		or 0.8 < INFO < 1.0
###	6	pink				#D5267B (213,38,123)				(0,82.16,42.25,16.47) 	=>	6		
###	7	hardpink			#CC0071 (204,0,113)					(0,0,0,0) 	=>	7		
###	8	lightpurple			#A8448A (168,68,138)				(0,0,0,0) 	=>	8		
###	9	purple				#9A3480 (154,52,128)				(0,0,0,0) 	=>	9		
###	10	lavendel			#8D5B9A (141,91,154)				(0,0,0,0) 	=>	10		
###	11	bluepurple			#705296 (112,82,150)				(0,0,0,0) 	=>	11		
###	12	purpleblue			#686AA9 (104,106,169)				(0,0,0,0) 	=>	12		
###	13	lightpurpleblue		#6173AD (97,115,173/101,120,180)	(0,0,0,0) 	=>	13		
###	14	seablue				#4C81BF (76,129,191)				(0,0,0,0) 	=>	14		
###	15	skyblue				#2F8BC9 (47,139,201)				(0,0,0,0) 	=>	15		
###	16	azurblue			#1290D9 (18,144,217)				(0,0,0,0) 	=>	16		 or 0.01 < MAF < 0.05 or 0.2 < INFO < 0.4
###	17	lightazurblue		#1396D8 (19,150,216)				(0,0,0,0) 	=>	17		
###	18	greenblue			#15A6C1 (21,166,193)				(0,0,0,0) 	=>	18		
###	19	seaweedgreen		#5EB17F (94,177,127)				(0,0,0,0) 	=>	19		
###	20	yellowgreen			#86B833 (134,184,51)				(0,0,0,0) 	=>	20		
###	21	lightmossgreen		#C5D220 (197,210,32)				(0,0,0,0) 	=>	21		
###	22	mossgreen			#9FC228 (159,194,40)				(0,0,0,0) 	=>	22		or MAF > 0.20 or 0.6 < INFO < 0.8
###	23	lightgreen			#78B113 (120,177,19)				(0,0,0,0) 	=>	23/X
###	24	green				#49A01D (73,160,29)					(0,0,0,0) 	=>	24/Y
###	25	grey				#595A5C (89,90,92)					(0,0,0,0) 	=>	25/XY	or MAF < 0.01 or 0.0 < INFO < 0.2
###	26	lightgrey			#A2A3A4	(162,163,164)				(0,0,0,0) 	=> 	26/MT
### 
### ADDITIONAL COLORS
### 27	midgrey				#D7D8D7
### 28	very lightgrey		#ECECEC
### 29	white				#FFFFFF
### 30	black				#000000
### --------------------------------------------------------------------------------------------------------------------


