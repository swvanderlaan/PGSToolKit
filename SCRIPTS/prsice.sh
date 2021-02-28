#!/bin/bash

# Created by		Sander W. van der Laan | UMC Utrecht | s.w.vanderlaan[at]gmail[dot]com
# Last edit			2018-10-08
# Version			1.1.0

echocyan "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echocyan "                                      CALCULATING PRS USING PRSICE2                                      "
echocyan ""


#HERCULESTOOLKIT="${SOFTWARE}/HerculesToolKit"

### Making the phenotype-file is project specific
# echoitalic " > getting a phenotype file and exclusion list."
# echo "FID IID AEGS_type COHORT STUDY_TYPE sex Age AgeSQR OR_year OR_year_C PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 Calcification_bin Collagen_bin Fat10_bin Fat40_bin Macrophages_bin SMC_bin IPH Macrophages_BC Mastcells_BC Neutrophils_BC SMC_BC VesselDensityAvg_BC" > ${PROJECTDIR}/aegscombo_phenocov.pheno
# cat ${ORIGINALDATA}/pheno_cov_exclusions/aegscombo_phenocov.sample | \
# parseTable --col ID_1,ID_2,AEGS_type,COHORT,STUDY_TYPE,sex,Age,AgeSQR,OR_year,OR_year_C,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,Calcification_bin,Collagen_bin,Fat10_bin,Fat40_bin,Macrophages_bin,SMC_bin,IPH,Macrophages_BC,Mastcells_BC,Neutrophils_BC,SMC_BC,VesselDensityAvg_BC | tail -n +3 | \
# sed 's/FEMALE/2/g' | sed 's/MALE/1/g' >> ${PROJECTDIR}/aegscombo_phenocov.pheno

### Be sure to make the proper exclusion list and make sure it has the headers FID IID!


### Binary trait
for PHENO in ${BINARYTRAITS[@]}; do 
	TARGETTYPE="T"
	PHENOTYPE="${PHENO}"
	PRSICEOUTPUTNAME="--out PRSice.${PHENOTYPE}.${SCORETYPE}.${DATATYPE}.${PERMUTATION}"
	
	cd ${PRSDIR}

	# prsice.R --prsice $(command -v prsice) \
	#the --extract option is for avoiding the 'duplicate SNPs' error
	#	--extract PRSice.Collagen_bin.PRSsum.BED.NOPERM_NOCENTER.valid \
	#	--extract PRSice.IPH.PRSsum.BED.NOPERM_NOCENTER.valid \

	${RSCRIPTPATH} ${PRSICE2R} --prsice ${PRSICE2SH} \
	--dir ${PRSDIR} \
	--seed ${PRSICESEED} \
	--bar-levels ${PRSICEBARLEVELS} \
	--base ${BASEDATA} \
	--target ${VALIDATIONDATA} \
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
	--cov-factor ${COVARIATESFACTOR} \
	--pheno-file ${PHENOTYPEFILE} \
	--pheno-col ${PHENOTYPE} \
	--remove ${EXCLUSION} \
	--type bgen \
	--extract PRSice.${PHENOTYPE}.PRSsum.BED.NOPERM_NOCENTER.valid \
	${PRSICEPLOTTING} \
	${PRSICESETTINGS} \
	${PRSICEOUTPUTNAME}

done

### Continuous trait
# for PHENO in Macrophages_BC Mastcells_BC Neutrophils_BC SMC_BC VesselDensityAvg_BC; do 
for PHENO in ${CONTINUOUSTRAITS[@]}; do
	TARGETTYPE="F"
	PHENOTYPE="${PHENO}"
	PRSICEOUTPUTNAME="--out PRSice.${PHENOTYPE}.${SCORETYPE}.${DATATYPE}.${PERMUTATION}"
	
	cd ${PRSDIR}

	prsice.R --prsice $(command -v prsice) \
	--dir ${PRSDIR} \
	--seed ${PRSICESEED} \
	--bar-levels ${PRSICEBARLEVELS} \
	--base ${BASEDATA} \
	--target ${VALIDATIONDATA} \
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
	--cov-factor ${COVARIATESFACTOR} \
	--pheno-file ${PHENOTYPEFILE} \
	--pheno-col ${PHENOTYPE} \
	--remove ${EXCLUSION} \
	${PRSICEPLOTTING} \
	${PRSICESETTINGS} \
	${PRSICEOUTPUTNAME}

done





# echobold "* Setting some variables specific for PRSice calculations."
		# echoitalic " > data files ..."
		# #BASEDATA="Inouye_bioRxiv_2018/metaGRS_hg19_20180205.cardiogramplusc4d_1kg_cad_add.pval.aegs_matched.txt"
		# BASEDATA="/hpc/dhl_ec/data/_gwas_datasets/_CARDIoGRAM/ukbiobank_cardiogramplusc4d/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz"
		# #TARGETDATA="/hpc/dhl_ec/data/_ae_originals/AEGS_COMBINED_IMPUTE2_1000Gp3_GoNL5/aegs_combo_1kGp3GoNL5_RAW_chr#" #" --type bgen "
		# TARGETDATA="/hpc/dhl_ec/data/_ae_originals/AEGS_COMBINED_EAGLE2_1000Gp3v5HRCr11/_chr23_1kg_gonl5/aegs.1kgp3gonl5.idfix.chr#"
		# echoitalic " > PRSice general and plotting settings ..."
		# PRSICESETTINGS="--no-clump --print-snp --all-score --extract PRSice.valid"
		# PRSICEPLOTTING="--bar-col-high \#E55738 --bar-col-low \#1290D9 --quantile 10 --quant-ref 1 --multi-plot 3"
		# PRSICETHREADS="4"
		# PRSICESEED="91149214"
		# PRSICEBARLEVELS="0.00000005,0.0000005,0.000005,0.00005,0.0005,0.005,0.05,0.1,0.2,0.3,0.4,0.5,1"

		# echoitalic " > statistics ..."
		# STATTYPE="--beta" 
		# TARGETTYPE="T,T,T,T,T,T,T,F,F,F,F,F"
		# SNPID="rsid"
		# CHRID="chr"
		# BPID="position"
		# A1ID="effect_allele"
		# A2ID="other_allele"
		# STATID="beta"
		# PVALUEID="P_cardiogramplusc4d_1kg_cad_add"

		# echoitalic " > phenotypes and covariates ..."
		# PHENOTYPE="Calcification_bin,Collagen_bin,Fat10_bin,Fat40_bin,Macrophages_bin,SMC_bin,IPH,Macrophages_BC,Mastcells_BC,Neutrophils_BC,SMC_BC,VesselDensityAvg_BC"
		# PHENOTYPEFILE="${PROJECTDIR}/aegscombo_phenocov.pheno"
		# COVARIATES="sex,Age,OR_year,@PC[1-4]"
		# COVARIATESFILE="${PROJECTDIR}/aegscombo_phenocov.pheno"
		# EXCLUSION="${PROJECTDIR}/aegscombo.exclusion_nonCEA.list"

		# cd ${PRSICEDIR}

		# echo "${SOFTWARE}/PRSice_linux_233/PRSice.R --prsice ${SOFTWARE}/PRSice_linux_233/PRSice_linux"


		# #prsice.R --prsice $(command -v prsice) \
		# #--base ${PROJECTDIR}/${BASEDATA} \

		# ${SOFTWARE}/R-3.6.3/R-3.6.3/bin/Rscript ${SOFTWARE}/PRSice_linux_233/PRSice.R --prsice ${SOFTWARE}/PRSice_linux_233/PRSice_linux \
		# --dir ${PRSICEDIR} \
		# --seed ${PRSICESEED} \
		# --bar-levels ${PRSICEBARLEVELS} \
		# --base ${BASEDATA} \
		# --type bgen
		# --target ${TARGETDATA} \
		# --thread ${PRSICETHREADS} \
		# ${STATTYPE} \
		# --binary-target ${TARGETTYPE} \
		# --snp ${SNPID} \
		# --chr ${CHRID} \
		# --bp ${BPID} \
		# --A1 ${A1ID} \
		# --A2 ${A2ID} \
		# --stat ${STATID} \
		# --pvalue ${PVALUEID} \
		# --cov-file ${COVARIATESFILE} \
		# --cov-col ${COVARIATES} \
		# --pheno-file ${PHENOTYPEFILE} \
		# --pheno-col ${PHENOTYPE} \
		# --remove ${EXCLUSION} \
		# ${PRSICEPLOTTING} \
		# ${PRSICESETTINGS}