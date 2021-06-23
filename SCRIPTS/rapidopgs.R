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
  make_option(c("-f", "--freqcol"), action="store", default=NA, type='character', help="Name of the effect allele frequency column in the base file"),
  make_option(c("-m", "--measurecol"), action="store", default=NA, type='character', help="Name of the beta/log(OR)/effect size column in the base file"),
  make_option(c("-e", "--secol"), action="store", default=NA, type='character', help="Name of the column of the standard error of the beta/log(OR)/effect size value in the base file"),
  make_option(c("-s", "--samples"), action="store", default=NA, type='character', help="Sample size of the GWAS"),
  make_option(c("-n", "--sbjcol"), action="store", default=NA, type='character', help="Name of the column containing the sample size for each variant"),
  make_option(c("-g", "--fthrs"), action="store", default=NA, type='character', help="scalar indicating the ppi threshold (if filt_threshold < 1) or the number of top SNPs by absolute weights (if filt_threshold >= 1) to filter the dataset after PGS computation. If NULL (DEFAULT), no thresholding will be applied."),
  make_option(c("-j", "--recalc"), action="store", default=NA, type='character', help="logical indicating if weights should be recalculated after thresholding, only relevant if filt_threshold is defined"),
  make_option(c("-l", "--trait"), action="store", default=NA, type='character', help="string specifying if the dataset corresponds to a case-control (\"cc\") or a quantitative trait (\"quant\") GWAS. If trait==quant, the ALT_FREQ column and sample size are required"),
  make_option(c("-v", "--pii"), action="store", default=NA, type='character', help="scalar representing the prior probability"),
  make_option(c("-x", "--prior"), action="store", default=NA, type='character', help="the prior specifies that BETA at causal SNPs follows a centred normal distribution with standard deviation sd.prior, sensible and widely used DEFAULTs are 0.2 for case control traits, and 0.15 * var(trait) for quantitative (selected if trait == \"quant\")"),
  make_option(c("-z", "--reffile"), action="store", default=NA, type='character', help="path to the reference file the SNPs should be filtered and aligned to, this file should have 5 columns (CHR, BP, SNPID, REF and ALT) and should be in the same build as the summary statistics")
  )
opt = parse_args(OptionParser(option_list=option_list))

print(opt)

# Parse some parameters to the required formats
if (!opt$sbjcol == "") {
  sample_par <- opt$sbjcol
} else if (!opt$samples == "") {
  sample_par <- opt$samples
} else {
  sample_par <- NULL
}

if (opt$fthrs == "") {
  opt$fthrs <- NULL
}

if (opt$recalc == "") {
  opt$recalc <- NULL
} else if (opt$recalc == "TRUE") {
  opt$recalc <- TRUE
} else if (opt$recalc == "FALSE") {
  opt$recalc <- FALSE
}

if (opt$pii == "") {
  opt$pii <- NULL
} else {
  opt$pii <- as.numeric(opt$pii)
}

if (opt$prior == "") {
  opt$prior <- NULL
} else {
  opt$prior <- as.numeric(opt$prior)
}

if (opt$reffile == "") {
  opt$reffile <- NULL
}

# Read the GWAS summary file to a table
sumstats <- read.table(opt$base, sep='\t', header=TRUE)

# Change the column names to the required format of RapidoPGS
colnames(sumstats)[colnames(sumstats) == opt$idcol] <- "SNPID"
colnames(sumstats)[colnames(sumstats) == opt$chrcol] <- "CHR"
colnames(sumstats)[colnames(sumstats) == opt$poscol] <- "BP"
colnames(sumstats)[colnames(sumstats) == opt$refcol] <- "REF"
colnames(sumstats)[colnames(sumstats) == opt$altcol] <- "ALT"
colnames(sumstats)[colnames(sumstats) == opt$measurecol] <- "BETA"
colnames(sumstats)[colnames(sumstats) == opt$secol] <- "SE"

if (!opt$freqcol == "") {
  colnames(sumstats)[colnames(sumstats) == opt$freqcol] <- "ALT_FREQ"
}

# RapidoPGS can only handle autosomes, therefore we remove SNPs on the X and Y chromosomes
sumstats <- subset(sumstats, CHR!="x" & CHR!="X" & CHR!="y" & CHR!="Y")

str(sumstats)

# Compute PGS using GWAS summary statistics
PGS <- rapidopgs_single(
  data = sumstats,
  N = opt$sample_par,
  trait = opt$trait,
  build = opt$build,
  pi_i = opt$pii,
  sd.prior = opt$prior,
  filt_threshold = opt$fthrs,
  recalc = opt$recalc,
  reference = opt$reffile
)

str(PGS)

# Write the PGS table to a file so it can be read by PLINK
export <- PGS[,c("SNPID", "ALT", "weight")]
write.table(export, file=opt$out, row.names=FALSE, sep="\t", quote=FALSE)
