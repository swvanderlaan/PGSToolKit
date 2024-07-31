#!/hpc/local/CentOS7/dhl_ec/software/R-3.6.3/R-3.6.3/bin/Rscript --vanilla

# R script for parsing the GWAS summary file to the default PRS-CS format.
# More on the PRS-CS format can be found here: https://github.com/getian107/PRScs

# Required packages
packages = c("optparse",
            "data.table")

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
option_list <- list(
    make_option(c("-i", "--input"), action="store", default=NA, type='character', help="Path to the base file"),
    make_option(c("-o", "--output"), action="store", default=NA, type='character', help="Path to where the output file should be written to"),
    make_option(c("-d", "--idcol"), action="store", default=NA, type='character', help="Name of the SNP ID column in the base file"),
    make_option(c("-r", "--refcol"), action="store", default=NA, type='character', help="Name of the reference allele column in the base file"),
    make_option(c("-a", "--altcol"), action="store", default=NA, type='character', help="Name of the alternative allele column in the base file"),
    make_option(c("-z", "--measure"), action="store", default=NA, type='character', help="Type of measure: either 'OR' or 'BETA'"),
    make_option(c("-m", "--measurecol"), action="store", default=NA, type='character', help="Name of the beta/log(OR)/effect size column in the base file"),
    make_option(c("-v", "--pvaluecol"), action="store", default=NA, type='character', help="Name of the p-value column in the base file"))

opt = parse_args(OptionParser(option_list=option_list))

print(opt)

# Read the GWAS summary file to a table
sumstats <- read.table(opt$input, sep='\t', header=TRUE, check.names=FALSE)

# Change the column names to the required format of PRS-CS
colnames(sumstats)[colnames(sumstats) == opt$idcol] <- "SNP"
colnames(sumstats)[colnames(sumstats) == opt$refcol] <- "A1"
colnames(sumstats)[colnames(sumstats) == opt$altcol] <- "A2"
colnames(sumstats)[colnames(sumstats) == opt$measurecol] <- toupper(opt$measure)
colnames(sumstats)[colnames(sumstats) == opt$pvaluecol] <- "P"

str(sumstats)

# Write the PGS table to a file so it can be read by PRS-CS
export <- sumstats[,c("SNP", "A1", "A2", toupper(opt$measure), "P")]
write.table(export, file=opt$output, row.names=FALSE, sep="\t", quote=FALSE)
