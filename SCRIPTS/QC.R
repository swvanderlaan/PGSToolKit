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
    make_option(c("-a", "--baseidcol"), action="store", default=NA, type='character', help="Index of the column in the base file containing the ID (start from 1)"),
    make_option(c("-d", "--statsidcolname"), action="store", default=NA, type='character', help="Index of the column in the stats file containing the ID (start from 1)"),
    make_option(c("-c", "--statsmafcolname"), action="store", default=NA, type='character', help="Index of the column in the stats file containing the MAF (start from 1)"),
    make_option(c("-t", "--statsinfocolname"), action="store", default=NA, type='character', help="Index of the column in the stats file containing the INFO score (start from 1)"),
    make_option(c("-r", "--result"), action="store", default=NA, type='character', help="Path to the file to which the QC results will be written")
  )
opt = parse_args(OptionParser(option_list=option_list))

# noquote("Quality control will now commence with the following parameters")
# noquote(paste("", "Minimum minor allele frequency:", opt$mafthres, sep="   "))
# noquote(paste("", "Minimum imputation (info) score:", opt$infothres, sep="   "))
# noquote("")

# Read the statsfile to a table
# opt$stats <- "/hpc/dhl_ec/aligterink/ProjectFiles/sumstatssample.txt"
snpstats <- read.table(opt$stats, sep='', header=TRUE, nrows=500)
snpstats <- snpstats[ , c(opt$statsidcol, opt$statsmafcol, opt$statsmafcol)]
colnames(snpstats)<-c("ID", "MAF", "INFO")

# Extract the IDs of the variants that meet the threshold
# snpstats[1] <- lapply(snpstats[1], as.character)
snpstatsIDs <- snpstats[snpstats$MAF>=opt$mafthres & snpstats$INFO>=opt$infothres , ][ ,1]

# Read the base file to a table
basetable <- read.table(opt$base, sep='', header=TRUE, nrows=500)

# Extract the base file variants of which the ID occurs in the list of statsfile variants that meet the thresholds
basetable[,"newid"] = paste(basetable$chr,":",basetable$bp_hg19,sep="")
QCd_variants <- basetable[basetable$newid %in% snpstatsIDs, ]

# Write the extracted variants to a new file
fwrite(QCd_variants, file=opt$out, row.names=FALSE, sep="\t", quote=FALSE)

# noquote(paste(nrow(QCd_variants), "out of", nrow(snpstats), "variants have met the thresholds"))
# noquote(paste(nrow(snpstats)-nrow(QCd_variants), "out of", nrow(snpstats), "variants were removed from the data"))

cat(paste(nrow(QCd_variants), " out of ", nrow(snpstats), " variants have met the thresholds\n", nrow(snpstats)-nrow(QCd_variants), " out of ", nrow(snpstats), "variants were removed from the data", sep=""), file=opt$result)




# alternate_ids rsid chromosome position alleleA alleleB index average_maximum_posterior_call info cohort_1_AA cohort_1_AB cohort_1_BB cohort_1_NULL all_AA all_AB all_BB all_NULL all_total all_maf missing_data_proportion cohort_1_hwe comment
# . 1:10177 01 10177 A AC 1 0.486761 0.00417734 718.628 1033.88 371.652 0 718.628 1033.88 371.652 0 2124 0.418326 -3.79002e-05 1 NA
# . 1:10235 01 10235 T TA 2 0.997695 0.00169929 2119.1 4.893 0 0.003 2119.1 4.893 0 0.003 2124 0.00115184 7.06215e-07 1 NA

# Markername      snptestid       chr     bp_hg19 effect_allele   noneffect_alleleeffect_allele_freq      logOR   se_gc   p-value_gc      n_samples       exome  info_ukbb
# 1:569406_G_A    rs561255355     1       569406  G       A       0.99858 0.051910.27358  0.849496        154959  yes     0.41
# 1:751756_T_C    rs143225517     1       751756  C       T       0.14418 0.005460.01416  0.699604        254199  no      0.98
# 1:753405_C_A    rs3115860       1       753405  C       A       0.17332 0.003160.0138   0.818922        255643  no      0.98


# v <- "MAF,19,>,0.005|info,9,>,0.3"
# l <- strsplit(v, "\\|")[[1]]

# x <- c()
# y <- "snpstatsIDs <- snpstats["

# for(i in 1:length(l)) {
#     # x[i] <- strsplit(l[i], ",")
#     a <- strsplit(l[i], ",")
#     y <- paste(y, "snpstats[, ", )
# }
