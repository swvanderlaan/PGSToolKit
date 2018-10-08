#!/bin/bash

# Created by		Sander W. van der Laan | UMC Utrecht | s.w.vanderlaan[at]gmail[dot]com
# Last edit			2018-10-08
# Version			1.1.0

### Creating display functions
### Setting colouring
NONE='\033[00m'
BOLD='\033[1m'
OPAQUE='\033[2m'
FLASHING='\033[5m'
UNDERLINE='\033[4m'

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

### echo -e "\033[1mbold\033[0m"
### echo -e "\033[3mitalic\033[0m" ### THIS DOESN'T WORK ON MAC!
### echo -e "\033[4munderline\033[0m"
### echo -e "\033[9mstrikethrough\033[0m"
### echo -e "\033[31mHello World\033[0m"
### echo -e "\x1B[31mHello World\033[0m"

function echocyan { #'echobold' is the function name
    echo -e "${CYAN}${1}${NONE}" # this is whatever the function needs to execute.
}
function echobold { #'echobold' is the function name
    echo -e "${BOLD}${1}${NONE}" # this is whatever the function needs to execute, note ${1} is the text for echo
}
function echoitalic { #'echobold' is the function name
    echo -e "\033[3m${1}\033[0m" # this is whatever the function needs to execute.
}

script_copyright_message() {
	echo ""
	THISYEAR=$(date +'%Y')
	echoitalic "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	echoitalic "+ The MIT License (MIT)                                                                                 +"
	echoitalic "+ Copyright (c) 1979-${THISYEAR} Sander W. van der Laan                                                        +"
	echoitalic "+                                                                                                       +"
	echoitalic "+ Permission is hereby granted, free of charge, to any person obtaining a copy of this software and     +"
	echoitalic "+ associated documentation files (the \"Software\"), to deal in the Software without restriction,         +"
	echoitalic "+ including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, +"
	echoitalic "+ and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, +"
	echoitalic "+ subject to the following conditions:                                                                  +"
	echoitalic "+                                                                                                       +"
	echoitalic "+ The above copyright notice and this permission notice shall be included in all copies or substantial  +"
	echoitalic "+ portions of the Software.                                                                             +"
	echoitalic "+                                                                                                       +"
	echoitalic "+ THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT     +"
	echoitalic "+ NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND                +"
	echoitalic "+ NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES  +"
	echoitalic "+ OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN   +"
	echoitalic "+ CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                            +"
	echoitalic "+                                                                                                       +"
	echoitalic "+ Reference: http://opensource.org.                                                                     +"
	echoitalic "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
}


echocyan "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echocyan "                           POLYGENIC SCORE CALCULATIONS"
echocyan ""
echocyan ""
echocyan "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echobold "*Setting the environment."

SOFTWARE="/hpc/local/CentOS7/dhl_ec/software/"
HERCULESTOOLKIT="${SOFTWARE}/HerculesToolKit"

### Project specific
PROJECTDIR="/hpc/dhl_ec/svanderlaan/projects/polygenicscores"
ORIGINALDATA="/hpc/dhl_ec/data/_ae_originals"

echobold "* Making some directories."
### Make directories for script if they do not exist yet (!!!PREREQUISITE!!!)
if [ ! -d ${PROJECTDIR}/PRSICE/ ]; then
	echo "The PRSICE directory does not exist, Mr. Bourne will make it for you!"
	mkdir -v ${PROJECTDIR}/PRSICE/
fi
PRSICEDIR=${PROJECTDIR}/PRSICE

chmod -R a+rwx ${PROJECTDIR}

echobold "* Setting some variables specific for PRSice calculations."
echoitalic " > data files ..."

### Fill in BASE and TARGET data-files
### Note: make sure you have the same variantID-nomenclature in both files...
BASEDATA="Inouye_bioRxiv_2018/metaGRS_hg19_20180205.cardiogramplusc4d_1kg_cad_add.pval.aegs_matched.txt"
TARGETDATA="${ORIGINALDATA}/AEGS_COMBINED_IMPUTE2_1000Gp3_GoNL5/aegs_combo_1kGp3GoNL5_RAW_chr# "

### Making the phenotype-file is project specific
echoitalic " > getting a phenotype file and exclusion list."
echo "FID IID AEGS_type COHORT STUDY_TYPE sex Age AgeSQR OR_year OR_year_C PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 Calcification_bin Collagen_bin Fat10_bin Fat40_bin Macrophages_bin SMC_bin IPH Macrophages_BC Mastcells_BC Neutrophils_BC SMC_BC VesselDensityAvg_BC" > ${PROJECTDIR}/aegscombo_phenocov.pheno
cat ${ORIGINALDATA}/pheno_cov_exclusions/aegscombo_phenocov.sample | \
parseTable --col ID_1,ID_2,AEGS_type,COHORT,STUDY_TYPE,sex,Age,AgeSQR,OR_year,OR_year_C,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,Calcification_bin,Collagen_bin,Fat10_bin,Fat40_bin,Macrophages_bin,SMC_bin,IPH,Macrophages_BC,Mastcells_BC,Neutrophils_BC,SMC_BC,VesselDensityAvg_BC | tail -n +3 | \
sed 's/FEMALE/2/g' | sed 's/MALE/1/g' >> ${PROJECTDIR}/aegscombo_phenocov.pheno

### 
echoitalic " > PRSice general and plotting settings ..."
### Specific PRSice settings
### --no-clump:  don't use clump if you already filtered the data; in case of most GWAS results
###              you do want to use clumping
### --print-snp: print a list of SNPs used in the end in the modeling
### --all-score: you want to print out all the calculated scores for each p-value threshold,
###              handy if you want to do some offline modeling down the road
### --extract:   the first time you run PRSice it will stop and state "you should --extract"
###              to include only the valid variants. You can re-run doing just that...
# PRSICESETTINGS="--no-clump --print-snp --extract PRSice.valid --score sum --missing center --all-score --perm 10000"
PRSICESETTINGS="--print-snp --extract PRSice.valid --score sum --all-score "
PRSICEPLOTTING="--fastscore --bar-col-high #E55738 --bar-col-low #1290D9 --quantile 100 --quant-break 2.5,5,10,20,40,60,80,90,95,97.5,100 --quant-ref 60"
PRSICETHREADS="4"
PRSICESEED="91149214" # just a random number
PRSICEBARLEVELS="0.00000005,0.0000005,0.000005,0.00005,0.0005,0.005,0.05,0.1,0.2,0.3,0.4,0.5,1"
SCORETYPE="PRSsum" # you can use different scoring algorithms, 'sum', or 'average', etc.
DATATYPE="BED" # you can use the Oxford-style BGEN format or the PLINK-style BED format
PERMUTATION="NOPERM_NOCENTER" # be explicit in the output-file name: did you use permutation, center?

echoitalic " > PRSice statistics settings ..."
STATTYPE="--beta" 

### Make sure these are the exact column-names in your BASEDATA
SNPID="rsid"
CHRID="chr"
BPID="position"
A1ID="effect_allele"
A2ID="other_allele"
STATID="beta"
PVALUEID="P_cardiogramplusc4d_1kg_cad_add"

echoitalic " > phenotypes and covariates ..."

PHENOTYPEFILE="${PROJECTDIR}/aegscombo_phenocov.pheno"
COVARIATES="sex,Age,OR_year,@PC[1-4],COHORT"
COVARIATESFACTOR="sex"
COVARIATESFILE="${PROJECTDIR}/aegscombo_phenocov.pheno"
EXCLUSION="${PROJECTDIR}/aegscombo.exclusion_nonCEA.list"
	
### Example binary trait
for PHENO in Calcification_bin Collagen_bin Fat10_bin Fat40_bin Macrophages_bin SMC_bin IPH; do 
	TARGETTYPE="T"
	PHENOTYPE="${PHENO}"
	PRSICEOUTPUTNAME="--out PRSice.${PHENOTYPE}.${SCORETYPE}.${DATATYPE}.${PERMUTATION}"
	
	cd ${PRSICEDIR}

	prsice.R --prsice $(command -v prsice) \
	--dir ${PRSICEDIR} \
	--seed ${PRSICESEED} \
	--bar-levels ${PRSICEBARLEVELS} \
	--base ${PROJECTDIR}/${BASEDATA} \
	--target ${TARGETDATA} \
	--thread ${PRSICETHREADS} \
	${STATTYPE} \
	--binary-target ${TARGETTYPE} \
	--snp ${SNPID} \
	--chr ${CHRID} \
	--bp ${BPID} \
	--A1 ${A1ID} \
	--A2 ${A2ID} \
	--stat ${STATID} \
	--pvalue ${PVALUEID} \
	--cov-file ${COVARIATESFILE} \
	--cov-col ${COVARIATES} \
	--pheno-file ${PHENOTYPEFILE} \
	--pheno-col ${PHENOTYPE} \
	--remove ${EXCLUSION} \
	${PRSICEPLOTTING} \
	${PRSICESETTINGS} \
	${PRSICEOUTPUTNAME}

done

### Example continuous trait
for PHENO in Macrophages_BC Mastcells_BC Neutrophils_BC SMC_BC VesselDensityAvg_BC; do 
	TARGETTYPE="F"
	PHENOTYPE="${PHENO}"
	PRSICEOUTPUTNAME="--out PRSice.${PHENOTYPE}.${SCORETYPE}.${DATATYPE}.${PERMUTATION}"
	
	cd ${PRSICEDIR}

	prsice.R --prsice $(command -v prsice) \
	--dir ${PRSICEDIR} \
	--seed ${PRSICESEED} \
	--bar-levels ${PRSICEBARLEVELS} \
	--base ${PROJECTDIR}/${BASEDATA} \
	--target ${TARGETDATA} \
	--thread ${PRSICETHREADS} \
	${STATTYPE} \
	--binary-target ${TARGETTYPE} \
	--snp ${SNPID} \
	--chr ${CHRID} \
	--bp ${BPID} \
	--A1 ${A1ID} \
	--A2 ${A2ID} \
	--stat ${STATID} \
	--pvalue ${PVALUEID} \
	--cov-file ${COVARIATESFILE} \
	--cov-col ${COVARIATES} \
	--pheno-file ${PHENOTYPEFILE} \
	--pheno-col ${PHENOTYPE} \
	--remove ${EXCLUSION} \
	${PRSICEPLOTTING} \
	${PRSICESETTINGS} \
	${PRSICEOUTPUTNAME}

done

echo ""
echocyan "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echobold "Wow. I'm all done buddy. What a job 😱 ! let's have a 🍻🍻 ... 🖖 "

script_copyright_message
