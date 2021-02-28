# Wrapper script for RapidoPGS
#
#   Reales, G., Kelemen, M., & Wallace, C. (2020). RápidoPGS: A rapid polygenic score 
#   calculator for summary GWAS data without validation dataset. BioRxiv.
#
# For additional documentation on RapidoPGS see 
# https://cran.r-project.org/web/packages/RapidoPGS/vignettes/Computing_RapidoPGS.html

# # Required packages
# packages = c("optparse",
#             "data.table")

## Load and install packages if necessary
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

#Parse arguments from commandline
args <- commandArgs(trailingOnly = TRUE)
option_list <- list(
    make_option(c("-b", "--basefile"), type = "character", dest = "base_file"),
    make_option(c("--sep"), type = "character", dest = "sep"),
    make_option(c("--build"), type = "character", dest = "build"), 
    make_option(c("--chr"), type = "character", dest = "chr_column"),
    make_option(c("--pos"), type = "character", dest = "pos_column"),
    make_option(c("--chr"), type = "character", dest = "snpID_column"),
    make_option(c("--ref"), type = "character", dest = "ref_column"),
    make_option(c("--alt"), type = "character", dest = "alt_column"),
    make_option(c("--measure"), type = "character", dest = "measure_column"), # This can be beta values, Log(OR)
    make_option(c("--SE"), type = "character", dest = "SE_column"),
    make_option(c("--pvalue"), type = "character", dest = "pvalue_column"),
    make_option(c("--altfreq"), type = "character", dest = "altfreq_column"))

# argv <- parse_args(OptionParser(option_list = option_list))

# Temp commandline args
# base_file <- "/home/dhl_ec/aligterink/ProjectFiles/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz"
# sep <- "\t"
# build <- "hg19" 
# chr_column <- "chr"
# pos_column <- "bp_hg19"
# snpID_column <- "snptestid"
# ref_column <- "effect_allele"
# alt_column <- "noneffect_allele"
# measure_column <- "logOR"
# SE_column <- "se_gc"
# pvalue_column <- "p-value_gc"
# altfreq_column <- "info_ukbb" # is this correct?

sumstats <- read.table(base_file, sep=sep, header=TRUE)

setnames(sumstats, old = c(snpID_column, chr_column, pos_column, ref_column, alt_column, measure_column, altfreq_column, SE_column, pvalue_column), 
    new=c("SNPID","CHR", "BP","REF", "ALT", "BETA", "ALT_FREQ", "SE", "P"))

# RapidoPGS can only handle autosomes, therefore we remove SNPs on the X or Y chromosomes
sumstats <- subset(sumstats, CHR!="x" & CHR!="X" & CHR!="y" & CHR!="Y")

# Compute PGS using GWAS summary statistics
PGS <- computePGS(
  data = sumstats,
  N0 = ,
  N1 = NULL,
  build = build,
  pi_i = 1e-04,
  sd.prior = if (is.null(N1)) {     0.15 } else {     0.2 },
  log.p = FALSE,
  filt_threshold = NULL,
  recalc = TRUE,
  reference = NULL,
  forsAUC = FALSE,
  altformat = FALSE
)
