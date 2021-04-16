#!/hpc/local/CentOS7/dhl_ec/software/R-3.6.3/R-3.6.3/bin/Rscript --vanilla
#
# R script for applying quality control to the dataset. 
#
# We will extract variants that pass the threshold from the stats file. 
# The IDs of these variants will be used to create a subset out of the variants in the GWAS summary file.
# The subset can then be used for further analysis.

# Required packages
packages = c("optparse",
            "data.table"
            )

## Load and install packages if necessary
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      try(library(x), silent=TRUE)
    }
  }
)

#Parse arguments from commandline
args <- commandArgs(trailingOnly = TRUE)
option_list = list(
    make_option(c("-b", "--base"), action="store", default=NA, type='character', help="Path to the GwAS summary file"),
    make_option(c("-s", "--stats"), action="store", default=NA, type='character', help="Path to the stats file"),
    make_option(c("-o", "--out"), action="store", default=NA, type='character', help="Path to the output file"),
    make_option(c("-m", "--mafthres"), action="store", default=NA, type='character', help="Minimum MAF"),
    make_option(c("-i", "--infothres"), action="store", default=NA, type='character', help="Minimum info (imputation) score"),
    make_option(c("-a", "--baseidcol"), action="store", default=NA, type='character', help="Name of the column in the base file containing the ID"),
    make_option(c("-d", "--statsidcolname"), action="store", default=NA, type='character', help="Name of the column in the stats file containing the ID"),
    make_option(c("-c", "--statsmafcolname"), action="store", default=NA, type='character', help="Name of the column in the stats file containing the MAF"),
    make_option(c("-t", "--statsinfocolname"), action="store", default=NA, type='character', help="Name of the column in the stats file containing the INFO score"),
    make_option(c("-r", "--result"), action="store", default=NA, type='character', help="Path to the file to which the QC results will be written")
  )
opt = parse_args(OptionParser(option_list=option_list))

# Read the statsfile to a table
snpstats <- read.table(opt$stats, sep='', header=TRUE)
snpstats <- snpstats[ , c(opt$statsidcol, opt$statsmafcol, opt$statsmafcol)]
colnames(snpstats)<-c("ID", "MAF", "INFO")

# Extract the IDs of the variants that meet the threshold
snpstatsIDs <- snpstats[snpstats$MAF>=opt$mafthres & snpstats$INFO>=opt$infothres , ][ ,1]

# Read the base file to a table
basetable <- read.table(opt$base, sep='', header=TRUE)

# Extract the base file variants of which the ID occurs in the list of statsfile variants that meet the thresholds
QCd_variants <- basetable[basetable$baseidcol %in% snpstatsIDs, ]

# These two lines can be used in case the stats file has chr:bp IDs and the base file has rs IDs 
# basetable[,"newid"] = paste(basetable$chr,":",basetable$bp_hg19,sep="")
# QCd_variants <- basetable[basetable$newid %in% snpstatsIDs, ]

# Write the extracted variants to a new file
fwrite(QCd_variants, file=opt$out, row.names=FALSE, sep="\t", quote=FALSE)

cat(paste(nrow(QCd_variants), " out of ", nrow(snpstats), " variants have met the thresholds\n", nrow(snpstats)-nrow(QCd_variants), " out of ", nrow(snpstats), "variants were removed from the data", sep=""), file=opt$result)
