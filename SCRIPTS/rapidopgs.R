#!/hpc/local/CentOS7/dhl_ec/software/R-3.6.3/R-3.6.3/bin/Rscript --vanilla

# Wrapper script for RapidoPGS
#
#   Reales, G., Kelemen, M., & Wallace, C. (2020). RápidoPGS: A rapid polygenic score 
#   calculator for summary GWAS data without validation dataset. BioRxiv.
#
# For additional documentation on RapidoPGS see 
# https://cran.r-project.org/web/packages/RapidoPGS/vignettes/Computing_RapidoPGS.html

# Required packages
packages = c("optparse",
            "data.table",
            "remotes",
            "RapidoPGS")

## Load and install packages if necessary
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE, repos='http://cran.us.r-project.org')
      try(library(x), silent=TRUE)
    }
  }
)

#Parse arguments from commandline
args <- commandArgs(trailingOnly = TRUE)
option_list = list(
  make_option(c("-b", "--base"), action="store", default=NA, type='character', help="Path to the base file, e.g. \'/home/documents/GWAS_summary.txt.gz\'. The file may be zipped"),
  make_option(c("-k", "--dir"), action="store", default=NA, type='character', help="Working directory for RapidoPGS"),
  make_option(c("-o", "--out"), action="store", default=NA, type='character', help="Path to the file in which the weights are stored"),
  make_option(c("-d", "--build"), action="store", default=NA, type='character', help="Build of the base file, e.g. \'hg19\' or \'hg38\'"),
  make_option(c("-i", "--idcol"), action="store", default=NA, type='character', help="Name of the SNP ID column in the base file"),
  make_option(c("-c", "--chrcol"), action="store", default=NA, type='character', help="Name of the chromosome column in the base file"),
  make_option(c("-p", "--poscol"), action="store", default=NA, type='character', help="Name of the position column in the base file"),
  make_option(c("-r", "--refcol"), action="store", default=NA, type='character', help="Name of the reference allele column in the base file"),
  make_option(c("-a", "--altcol"), action="store", default=NA, type='character', help="Name of the alternative allele column in the base file"),
  make_option(c("-f", "--freqcol"), action="store", default=NA, type='character', help="Name of the allele frequency column in the base file"),
  make_option(c("-w", "--whichfrq"), action="store", default=NA, type='character', help="'effect' if the --freqcol points to the effect allele frequence, 'alt' if it points to the alternative allele frequency"),
  make_option(c("-m", "--measurecol"), action="store", default=NA, type='character', help="Name of the beta/log(OR)/effect size column in the base file"),
  make_option(c("-e", "--secol"), action="store", default=NA, type='character', help="Name of the column of the standard error of the beta/log(OR)/effect size value in the base file"),
  make_option(c("-s", "--samples"), action="store", default=NA, type='character', help="Sample size of the GWAS"),
  make_option(c("-g", "--fthrs"), action="store", default=NA, type='character', help="scalar indicating the ppi threshold (if filt_threshold < 1) or the number of top SNPs by absolute weights (if filt_threshold >= 1) to filter the dataset after PGS computation. If NULL (DEFAULT), no thresholding will be applied."),
  make_option(c("-j", "--recalc"), action="store", default=NA, type='character', help="logical indicating if weights should be recalculated after thresholding, only relevant if filt_threshold is defined"),
  make_option(c("-l", "--trait"), action="store", default=NA, type='character', help="string specifying if the dataset corresponds to a case-control (\"cc\") or a quantitative trait (\"quant\") GWAS. If trait==quant, an ALT_FREQ column is required"),
  make_option(c("-v", "--pii"), action="store", default=NA, type='character', help="scalar representing the prior probability"),
  make_option(c("-x", "--prior"), action="store", default=NA, type='character', help=""),
  make_option(c("-z", "--ref"), action="store", default=NA, type='character', help="")
  )
opt = parse_args(OptionParser(option_list=option_list))

# Read the GWAS summary file to a table
sumstats <- read.table(opt$base, sep='\t', header=TRUE)

# Change the column names to the required format of RapidoPGS
colnames(sumstats)[colnames(sumstats) == opt$idcol] <- "SNPID"
colnames(sumstats)[colnames(sumstats) == opt$chrcol] <- "CHR"
colnames(sumstats)[colnames(sumstats) == opt$poscol] <- "BP"
colnames(sumstats)[colnames(sumstats) == opt$refcol] <- "REF"
colnames(sumstats)[colnames(sumstats) == opt$altcol] <- "ALT"
colnames(sumstats)[colnames(sumstats) == opt$freqcol] <- "ALT_FREQ"
colnames(sumstats)[colnames(sumstats) == opt$measurecol] <- "BETA"
colnames(sumstats)[colnames(sumstats) == opt$secol] <- "SE"

# RapidoPGS requires the effect allele frequency. If the frequency column points to the non-effect allele frequency
# we will do '1 - non-effect allele frequency' to calculate the effect allele frequency
if (opt$whichfrq == 'noneffect') {
  sumstats$ALT_FREQ <- (1-sumstats$ALT_FREQ)  
} else if (opt$whichfrq == 'effect') {
  print('Frequency column points to the effect allele')
} else {
  print('ERROR: WFRQ parameter not recognized')
  quit(save="no", status=1)
}

# RapidoPGS can only handle autosomes, therefore we remove SNPs on the X and Y chromosomes
sumstats <- subset(sumstats, CHR!="x" & CHR!="X" & CHR!="y" & CHR!="Y")

# Parse some parameters to the required formats
if (is.null(opt$samples) || all(opt$samples == "")) {
  opt$samples <- NULL
}
if (is.null(opt$fthrs) || all(opt$fthrs == "")) {
  opt$fthrs <- NULL
}
if (is.null(opt$sbjcol)) {
  opt$recalc <- NULL
} 
else if (all(opt$recalc == "TRUE")) {
  opt$recalc <- TRUE
} 
else if (all(opt$recalc == "FALSE")) {
  opt$recalc <- FALSE
}
if (is.null(opt$ref) || all(opt$ref == "")) {
  opt$ref <- NULL
}
#  pi_i = as.numeric(opt$pii),
#   sd.prior = as.numeric(opt$prior),
#   pi_i = NULL,
#   sd.prior = NULL,

# Compute PGS using GWAS summary statistics
# PGS <- rapidopgs_single(
#   data = sumstats,
#   N = opt$sbjcol,
#   trait = opt$trait,
#   build = opt$build,

#   filt_threshold = opt$fthrs,
#   recalc = opt$recalc,
#   reference = opt$ref
# )


PGS <- rapidopgs_single(
  data = sumstats,
  N = opt$samples,
  trait = opt$trait,
  build = opt$build,
  pi_i = 1e-04,
  sd.prior = 0.15,
  filt_threshold = NULL,
  recalc = TRUE,
  reference = NULL
)

# Write the PGS table to a file so it can be read by PLINK
export <- PGS[,c("SNPID", "REF", "WEIGHT")]
write.table(export, file=opt$out, row.names=FALSE, sep="\t", quote=FALSE)
