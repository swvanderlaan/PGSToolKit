#!/hpc/local/CentOS7/dhl_ec/software/R-3.6.3/R-3.6.3/bin/Rscript --vanilla

### This script is still a WIP.

# Required packages
packages = c("optparse",
            "data.table",
            "remotes",
            "bigsnpr",
            "doParallel",
            "RSQLite",
            "dbplyr",
            "stringr",
            "rlist")

## Load and install packages if necessary
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only=TRUE)) {
      install.packages(x, dependencies=TRUE, repos='http://cran.us.r-project.org', quiet=TRUE)
      library(x, quietly=TRUE)
    }
  }
)

# tmp pars
generate_notfound <- "FALSE"
sumstats_path <- "/hpc/dhl_ec/aligterink/ProjectFiles/base_data/CAD_UKBB_QCd_maf0.01_info0.8.txt.gz"
prsdir <- "/hpc/dhl_ec/aligterink/ProjectFiles/tmp/ldpredtmp"

bf_chr_col <- "chr"
bf_pos_col <- "bp_hg19"
bf_effect_allele_col <- "effect_allele"
bf_non_effect_allele_col <- "noneffect_allele"
bf_stat_col <- "logOR"
bf_se_col <- "se_gc"
bf_samples_col <- "n_samples"

bgen_object_available <- "TRUE"
bgen_object_path <- "/hpc/dhl_ec/aligterink/ProjectFiles/backingfile" #paste0(prsdir, "/backingfile_all_bgen.rds)"

exclude_snps <- "TRUE" #/hpc/dhl_ec/aligterink/ProjectFiles/aegs.qc.1kgp3hrcr11.idfix.rsid.chr21_8bit_not_found.rds
exclude_snps_dir <- "/hpc/dhl_ec/aligterink/ProjectFiles/bgen_notfound_info08_maf001"

bgen_folder <- "/hpc/dhl_ec/data/_ae_originals/AEGS_COMBINED_EAGLE2_1000Gp3v5HRCr11"
bgen_prefix <- "aegs.qc.1kgp3hrcr11.idfix.rsid.8bit.chr"

reference_data <- "/hpc/dhl_ec/aligterink/ProjectFiles/reference_data/PRScs/ldblk_1kg_eur/snpinfo_1kg_hm3"

cores <- nb_cores()
if (cores < 1) {
  cores <- 1
}

# Read sumstats
sumstats <- bigreadr::fread2(sumstats_path)

# Change column names to the required format
colnames(sumstats)[colnames(sumstats) == bf_chr_col] <- "chr"
colnames(sumstats)[colnames(sumstats) == bf_pos_col] <- "pos"
colnames(sumstats)[colnames(sumstats) == bf_effect_allele_col] <- "a0"
colnames(sumstats)[colnames(sumstats) == bf_non_effect_allele_col] <- "a1"
colnames(sumstats)[colnames(sumstats) == bf_stat_col] <- "beta"
colnames(sumstats)[colnames(sumstats) == bf_se_col] <- "beta_se"
colnames(sumstats)[colnames(sumstats) == bf_samples_col] <- "n_eff"

sumstats <- subset(sumstats, select=c("chr", "pos", "a0", "a1", "beta", "beta_se", "n_eff"))

# Attach bgen object
bgen_object <- snp_attach(paste(bgen_object_path, ".rds", sep=""))

# Set aliases
G <- bgen_object$genotypes
CHR <- as.integer(bgen_object$map$chromosome)
POS <- as.integer(bgen_object$map$physical.pos)
# POS2 <- snp_asGeneticPos(CHR, POS, dir=prsdir, ncores=cores)

# CHR_test <- CHR[270937:270946]
# POS_test <- POS[270937:270946]
# print(CHR_test)
# print(POS_test)

# POS2_test <- snp_asGeneticPos(CHR_test, POS_test, dir="/hpc/dhl_ec/aligterink/ProjectFiles", ncores=cores)
# print(POS2_test)


# Map object
map <- bgen_object$map[-(2:3)]
names(map) <- c("chr", "pos", "a0", "a1")
map$chr <- as.integer(map$chr)

# Match variants
info_snp <- snp_match(sumstats, map)

# Generate LD matrix
for (chr in 1:22) {

  print(chr)

  ind.chr <- which(info_snp$chr == chr)
  ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
  ind.chr3 <- match(ind.chr2, which(CHR == chr))

  corr0 <- snp_cor(G, ind.row=1:2124, ind.col=ind.chr, alpha=1, infos.pos=POS2[ind.chr], size=3/1000, ncores=cores)

  if (chr == 1) {
    df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")]
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, sfbm_backingfile)
  } else {
    df_beta <- rbind(df_beta, info_snp[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")])
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}

str(corr)


# # Read sumstats
# sumstats <- bigreadr::fread2(sumstats_path)

# # Change column names to the required format
# colnames(sumstats)[colnames(sumstats) == bf_chr_col] <- "chr"
# colnames(sumstats)[colnames(sumstats) == bf_pos_col] <- "pos"
# colnames(sumstats)[colnames(sumstats) == bf_effect_allele_col] <- "a0"
# colnames(sumstats)[colnames(sumstats) == bf_non_effect_allele_col] <- "a1"
# colnames(sumstats)[colnames(sumstats) == bf_stat_col] <- "beta"
# colnames(sumstats)[colnames(sumstats) == bf_se_col] <- "beta_se"
# colnames(sumstats)[colnames(sumstats) == bf_samples_col] <- "n_eff"

# sumstats <- subset(sumstats, select=c("chr", "pos", "a0", "a1", "beta", "beta_se", "n_eff"))

# # Attach bgen object
# bgen_object <- snp_attach(paste(bgen_object_path, ".rds", sep=""))

# # Set aliases
# G <- bgen_object$genotypes
# CHR <- as.integer(bgen_object$map$chromosome)
# POS <- as.integer(bgen_object$map$physical.pos)
# POS2 <- snp_asGeneticPos(CHR, POS, dir=prsdir, ncores=cores)

# # Map object
# map <- bgen_object$map[-(2:3)]
# names(map) <- c("chr", "pos", "a0", "a1")
# map$chr <- as.integer(map$chr)

# # Match variants
# info_snp <- snp_match(sumstats, map)

# # Generate LD matrix
# for (chr in 1:22) {

#   print(chr)

#   ind.chr <- which(info_snp$chr == chr)
#   ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
#   ind.chr3 <- match(ind.chr2, which(CHR == chr))

#   str(ind.chr)
#   str(ind.chr2)
#   str(ind.chr3)

#   df <- data.frame(POS[which(CHR == chr)], POS2[which(CHR == chr)])
#   write.table(df, file=paste0("/hpc/dhl_ec/aligterink/ProjectFiles/poging2_POS_POS2_chr", as.character(chr), ".txt"), row.names=FALSE, sep="\t", quote=FALSE)

#   corr0 <- snp_cor(G, ind.row=1:2124, ind.col=ind.chr, alpha=1, infos.pos=POS2[ind.chr], size=3/1000, ncores=cores)
#   str(corr0)
#   # corr0 <- corr0[ind.chr3, ind.chr3]
#   # str(corr0)

#   if (chr == 1) {
#     df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")]
#     ld <- Matrix::colSums(corr0^2)
#     corr <- as_SFBM(corr0, paste0(prsdir, "/backingfile_testx"))
#   } else {
#     df_beta <- rbind(df_beta, info_snp[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")])
#     ld <- c(ld, Matrix::colSums(corr0^2))
#     corr$add_columns(corr0, nrow(corr))
#   }
# }

# str(corr)






# # Read base file
# sumstats <- bigreadr::fread2(sumstats_path)

# # Change column names to the required format
# colnames(sumstats)[colnames(sumstats) == bf_chr_col] <- "chr"
# colnames(sumstats)[colnames(sumstats) == bf_pos_col] <- "pos"
# colnames(sumstats)[colnames(sumstats) == bf_effect_allele_col] <- "a0" # should effect col be a0?
# colnames(sumstats)[colnames(sumstats) == bf_non_effect_allele_col] <- "a1" # should non effect be a1?
# colnames(sumstats)[colnames(sumstats) == bf_stat_col] <- "beta"
# colnames(sumstats)[colnames(sumstats) == bf_se_col] <- "beta_se"
# colnames(sumstats)[colnames(sumstats) == "n_samples"] <- "n_eff" # not sure if this is ok

# sumstats <- subset(sumstats, select=c("chr", "pos", "a0", "a1", "beta", "beta_se", "n_eff"))
# # sumstats$chr <- as.character(sumstats$chr)


# # Unfortunately, the bigsnpr::snp_readBGEN function errors out when it finds base file variants that do not exist in 
# # the bgen files. When this happens it generates a file for the bgen file it is currently reading containing the variants 
# # that were not found. This file can be found in the same directory as the bgen file with "_not_found.rds" at the end. 
# # In order to generate a "_not_found.rds" file for each bgen file, we need to try to read each bgen file separately, 
# # which is done if "generate_notfound" is set to true. These files will then be used to filter out base file variants
# # on two conditions: 1) "exclude_snps_path" points to the directory containing these files (they are generated in the 
# # same directory as the bgen files) and 2) the names of these files are unaltered, meaning they are identical except 
# # for having "_not_found.rds" at the end instead of ".bgen".
# if (generate_notfound == "TRUE") {
#   gen_not_found <- function(bgenfile) {
    
#     # Extract the chromosome from the bgen file name
#     bgen_chr <- str_match(bgenfile, paste0(bgen_prefix, "(.*).bgen"))[2]
#     snp_list_index <- which(sapply(snp_list, function(x) x[[1]]) == bgen_chr)

#     # Path for the temporary backingfile
#     tmp_backingfile <- paste0(prsdir, "/chr", bgen_chr, "_tmp_backingfile")

#     print(paste("Generating a _not_found.rds file for chromosome", bgen_chr))

#     # Run the readBGEN function in order to generate a "_not_found.rds" file for later variant removal
#     tryCatch(                      
#       expr = {                     
#         bigsnpr::snp_readBGEN(
#           bgenfiles   = c(bgenfile),
#           backingfile = tmp_backingfile,
#           list_snp_id = snp_list[[snp_list_index]][2],
#           ncores      = cores
#         )
#         print("_not_found.rds file was not generated, likely because all base file variants occur in the bgen file")
#       },
#       error = function(e){          
#         print(e)
#     })
#     print("")    

#     # Remove the backingfile that was temporarily generated
#     if (file.exists(tmp_backingfile)) {
#       file.remove(tmp_backingfile)
#     }
#   }

#   lapply(bgenfiles, gen_not_found)
# }


# # Check if the bgen files are already stored in an object
# if (bgen_object_available == "FALSE") {

#   # The goal here is to make a list of lists where each of the inner lists contains the variants for a 
#   # single chromosome in the following format: "chr_pos_a0_a1", e.g. "1_753405_C_A". Important is that 
#   # the amount of inner lists is equal to the amount of bgen files, otherwise we will get an 'incompatibility 
#   # between dimensions' error. I.e, make sure each bgen file has matching variants in the sumstats file. 
#   snp_list <- lapply(split(sumstats, sumstats$chr), function(x) x <- list(as.character(x$chr[1]), paste(gsub("^0", "", x$chr), x$pos, x$a0, x$a1, sep="_")))
#   print("Extracted variants from the base file:")
#   str(snp_list)
#   print("")

#   # Filter out variants that do not occur in the bgen files using the "_not_found.rds" files
#   if (exclude_snps == "TRUE") {
    
#     for (bgenfile in bgenfiles) {
#       bgen_chr <- str_match(bgenfile, paste0(bgen_prefix, "(.*).bgen"))[2]
#       exclusion_file <- paste0(exclude_snps_dir, "/", bgen_prefix, bgen_chr, "_not_found.rds")
#       snp_list_index <- which(sapply(snp_list, function(x) x[[1]]) == bgen_chr)

#       if (file.exists(exclusion_file)) {
#         exclusion_snps <- readRDS(exclusion_file)

#         # Remove the '0' that LDpred likes to write in front of chromosomes
#         exclusion_snps <- gsub("^0", "", exclusion_snps)

#         print(paste(as.character(length(exclusion_snps)), "variants from", exclusion_file, "will be excluded for chromosome", bgen_chr))

#         filtered_snps[[length(filtered_snps)+1]] <- snp_list[[snp_list_index]][[2]][! snp_list[[snp_list_index]][[2]] %in% exclusion_snps]

#       } else {
#         print(paste("No _not_found.rds file found for ", bgenfile))
#         filtered_snps[[length(filtered_snps)+1]] <- snp_list[[snp_list_index]][[2]]
#       }
#     }
#     snp_list <- filtered_snps

#     str(snp_list)
#   }
#   # Parse the bgen files to a .rds file
#   bigsnpr::snp_readBGEN(
#     bgenfiles   = bgenfiles,
#     backingfile = bgen_objects_path,
#     list_snp_id = snp_list,
#     ncores      = cores
#   )
  
#   print(paste("Bgen file were stored in ", bgen_object_path, ".rds", sep=""))

# } else if (bgen_object_available == "TRUE") {
#   print(paste("Using the bgen file stored in ", bgen_object_path, ".rds", sep=""))

# } else {
#   print("BGEN_OBJECT_AVAILABLE parameter not recognized, exiting ...")
#   quit(status=1)
# }

# # Read bgen object
# bgen_object <- snp_attach(paste(bgen_object_path, ".rds", sep=""))
# G <- bgen_object$genotypes
# CHR <- as.integer(bgen_object$map$chromosome)
# POS <- as.integer(bgen_object$map$physical.pos)
# POS2 <- snp_asGeneticPos(CHR, POS, dir=prsdir, ncores=cores)

# # map contains a list of variants present in the bgen files
# map <- bgen_object$map[-(2:3)]
# names(map) <- c("chr", "pos", "a0", "a1")
# map$chr <- as.integer(map$chr)

# # map_chr2 <- map[which(CHR == 2), ]
# # write.table(map_chr2, file="/hpc/dhl_ec/aligterink/ProjectFiles/map_chr2.txt", row.names=FALSE, sep="\t", quote=FALSE)

# # Match variants between summary statistics and genotype data
# info_snp <- snp_match(sumstats, map)

# # prt <- which(CHR == 2)
# # write.table(prt, file="/hpc/dhl_ec/aligterink/ProjectFiles/ind.chr_chr2.txt", row.names=FALSE, sep="\t", quote=FALSE)



# # chr <- 2
# # print(chr)

# # ind.chr <- which(CHR == chr)
# # corr0 <- snp_cor(G, ind.row=1:2124, ind.col=ind.chr, alpha=1, infos.pos=POS2[ind.chr], size=3/1000, ncores=cores)
# # str(corr0)

# for (chr in 1:22) {

#   print(chr)

#   ind.chr <- which(info_snp$chr == chr)
#   ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
#   ind.chr3 <- match(ind.chr2, which(CHR == chr))

#   str(ind.chr)
#   str(ind.chr2)
#   str(ind.chr3)

#   df <- data.frame(POS[which(CHR == chr)], POS2[which(CHR == chr)])
#   write.table(df, file=paste0("/hpc/dhl_ec/aligterink/ProjectFiles/poging2_POS_POS2_chr", as.character(chr), ".txt"), row.names=FALSE, sep="\t", quote=FALSE)

#   corr0 <- snp_cor(G, ind.row=1:2124, ind.col=ind.chr2, alpha=1, infos.pos=POS2[ind.chr], size=3/1000, ncores=cores)
#   str(corr0)
#   corr0 <- corr0[ind.chr3, ind.chr3]
#   str(corr0)

#   if (chr == 1) {
#     df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")]
#     ld <- Matrix::colSums(corr0^2)
#     corr <- as_SFBM(corr0, paste0(prsdir, "/backingfile_test200"))
#   } else {
#     df_beta <- rbind(df_beta, info_snp[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")])
#     ld <- c(ld, Matrix::colSums(corr0^2))
#     corr$add_columns(corr0, nrow(corr))
#   }
# }

# str(corr)















# ## Information for the variants provided in the LD reference
# # map_ldref <- readRDS("/hpc/dhl_ec/aligterink/ProjectFiles/reference_data/LDpred/map.rds")
# map_ldref <- bigreadr::fread2("/hpc/dhl_ec/aligterink/ProjectFiles/reference_data/PRScs/ldblk_1kg_eur/snpinfo_1kg_hm3")
# colnames(map_ldref)[colnames(map_ldref) == "CHR"] <- "chr"
# colnames(map_ldref)[colnames(map_ldref) == "BP"] <- "pos"
# colnames(map_ldref)[colnames(map_ldref) == "A1"] <- "a0"
# colnames(map_ldref)[colnames(map_ldref) == "A2"] <- "a1"
# str(map_ldref)

# info_snp <- snp_match(sumstats, map_ldref)

# df_beta <- info_snp


# for (chr in 1:22) {

#   cat(chr, ".. ", sep = "")

#   ## indices in 'df_beta'
#   ind.chr <- which(df_beta$chr == chr)
#   ## indices in 'map_ldref'
#   ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
#   ## indices in 'corr_chr'
#   ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))

#   corr_file <- paste0("/hpc/dhl_ec/aligterink/ProjectFiles/reference_data/PRScs/ldblk_1kg_eur/ldblk_1kg_chr", chr, ".hdf5")
#   # corr_chr <- readRDS(paste0("/hpc/dhl_ec/aligterink/ProjectFiles/reference_data/LDpred/LD_chr", chr, ".rds"))[ind.chr3, ind.chr3]
#   # corr_chr <- readRDS("/hpc/dhl_ec/aligterink/ProjectFiles/reference_data/LDpred/LD_chr1.rds")
#   # corr_chr <- bigreadr::fread2(paste0("/hpc/dhl_ec/aligterink/ProjectFiles/reference_data/PRScs/ldblk_1kg_eur/ldblk_1kg_chr", chr, ".hdf5"))[ind.chr3, ind.chr3]
#   # corr_chr <- h5file(paste0("/hpc/dhl_ec/aligterink/ProjectFiles/reference_data/PRScs/ldblk_1kg_eur/ldblk_1kg_chr", chr, ".hdf5"), mode = "r")[ind.chr3, ind.chr3]
#   # str(corr_chr)



#   str(h5ls(corr_file))
  
#   str(h5read(corr_file))


#   # if (chr == 1) {
#   #   corr <- as_SFBM(corr_chr, "/hpc/dhl_ec/aligterink/ProjectFiles/tmp/ldpredtmp/tmp_backingfile_sfbm")
#   # } else {
#   #   corr$add_columns(corr_chr, nrow(corr))
#   # }
# }

# str(corr)


















# # Read bgen object
# bgen_object <- snp_attach(paste(bgen_object_path, ".rds", sep=""))
# G <- bgen_object$genotypes
# CHR <- as.integer(bgen_object$map$chromosome)
# POS <- as.integer(bgen_object$map$physical.pos)
# POS2 <- snp_asGeneticPos(CHR, POS, dir=prsdir, ncores=cores)

# # map contains a list of variants present in the bgen files
# map <- bgen_object$map[-(2:3)]
# names(map) <- c("chr", "pos", "a0", "a1")
# map$chr <- as.integer(map$chr)

# # Match variants between summary statistics and genotype data
# info_snp <- snp_match(sumstats, map)

# head(info_snp[which(chr == 2)], )
# tail(info_snp[which(chr == 2)], )

# info_snp <- info_snp[ with(info_snp, order(info_snp$pos)), ]

# head(info_snp[which(chr == 2)], )
# tail(info_snp[which(chr == 2)], )

# chr <- 2
# ind.chr <- which(info_snp$chr == chr)
# ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
# ind.chr3 <- match(ind.chr2, which(CHR == chr))

# options("scipen"=100, "digits"=4)
# print(ind.chr)
# print(POS[ind.chr])
# print(POS2[ind.chr])


# for (chr in 1:22) {
#     chr <- 2

#     # print(chr)

#     # indices in 'sumstats'
#     ind.chr <- which(info_snp$chr == chr)

#     # indices in 'G'
#     ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]

#     # indices in 'corr'
#     ind.chr3 <- match(ind.chr2, which(CHR == chr))

#     # str(ind.chr)
#     # str(ind.chr2)
#     # str(ind.chr3)

#     # str(G)
#     # str(POS2)

#     options("scipen"=100, "digits"=4)

#     print(ind.chr)
#     print(POS[ind.chr])
#     print(POS2[ind.chr])

#     # Get the Pearson correlations of nearby SNPs on the same chromosome
#     # snp_cor(G, ind.row=ind.val, ind.col=ind.chr, alpha=1, infos.pos=POS2[ind.chr], size=3/1000, ncores=cores),
#     cor_matrix <- snp_cor(G, ind.row=1:2124, ind.col=ind.chr, alpha=1, infos.pos=POS2[ind.chr], size=3/1000, ncores=cores)
#     # cor_matrix <- snp_cor(G, ind.row=1:20, ind.col=1:19, alpha=1, infos.pos=POS2[1:19], size=3/1000, ncores=cores)

#     str(cor_matrix) # is 233301 x 233301
    
#     print('1')
#     # corr0 <- cor_matrix[ind.chr3, ind.chr3]
#     # corr0 <- cor_matrix[1:19, 1:19]
#     print('1a')

#     corr0 <- cor_matrix
    
#     print('2')

#     if (chr == 1) {
#       df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")]
#       print('3')
#       ld <- Matrix::colSums(corr0^2)
#       print('4')
#       corr <- as_SFBM(corr0, paste0(prsdir, "/SFBM_backingfile"))
#       print('5')
#     } else {
#       df_beta <- rbind(df_beta, info_snp[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")])
#       ld <- c(ld, Matrix::colSums(corr0^2))
#       corr$add_columns(corr0, nrow(corr))
#     }

#     str(corr)
#     saveRDS(corr, file="/hpc/dhl_ec/aligterink/ProjectFiles/corr_matrix.rds")
# }






# map <- obj.bigSNP$map[-(2:3)]
# names(map) <- c("chr", "pos", "a0", "a1")
# map$chr <- as.numeric(map$chr)

# # Match variants between summary statistics and genotype data
# info_snp <- snp_match(sumstats, map)

# G <- obj.bigSNP$genotypes
# CHR <- map$chr
# POS <- map$pos
# POS2 <- snp_asGeneticPos(CHR, POS, dir=prsdir, ncores=cores)
# chr_set <- as.set(CHR)
# print(chr_set)






# # Check if the bgen files are already stored in objects
# if (bgen_objects_available == "FALSE") {

#   bgen_objects_path <- prsdir
#   bgenfiles <- list.files(path=bgen_folder, pattern=paste0("^", bgen_prefix, "(.*)\\.bgen$"), full.names=TRUE)
#   print("Found .bgen files:")
#   print(bgenfiles)
#   print("")

#   # The goal here is to make a list of lists where each of the inner lists contains the variants for a 
#   # single chromosome in the following format: "chr_pos_a0_a1", e.g. "1_753405_C_A". Important is that 
#   # the amount of inner lists is equal to the amount of bgen files, otherwise we will get an 'incompatibility 
#   # between dimensions' error. I.e, make sure each bgen file has matching variants in the sumstats file. 
#   snp_list <- lapply(split(sumstats, sumstats$chr), function(x) x <- list(as.character(x$chr[1]), paste(gsub("^0", "", x$chr), x$pos, x$a0, x$a1, sep="_")))
#   print("Extracted variants from the base file:")
#   str(snp_list)
#   print("")

#   for (bgenfile in bgenfiles) {
#     bgen_chr <- str_match(bgenfile, paste0(bgen_prefix, "(.*).bgen"))[2]
#     chr_snps <- snp_list[[which(sapply(snp_list, function(x) x[[1]]) == bgen_chr)]][[2]]
#     chr_snps_length <- length(chr_snps)
#     chr_backingfile <- paste0(bgen_objects_path, "/backingfile_chr", bgen_chr)


#     # Filter out variants that do not occur in the bgen files using the "_not_found.rds" files
#     if (exclude_snps == "TRUE") {

#       exclusion_file <- paste0(exclude_snps_dir, "/", bgen_prefix, bgen_chr, "_not_found.rds")

#       if (file.exists(exclusion_file)) {
#         exclusion_snps <- readRDS(exclusion_file)

#         # Remove the '0' that LDpred likes to write in front of chromosomes
#         exclusion_snps <- gsub("^0", "", exclusion_snps)

#         chr_snps <- list(chr_snps[! chr_snps %in% exclusion_snps])

#       } else {
#         print(paste("No _not_found.rds file found for ", bgenfile))
#       }
#     }
  
#     # Parse the bgen file to a .rds file
#     bigsnpr::snp_readBGEN(
#       bgenfiles   = c(bgenfile),
#       backingfile = chr_backingfile,
#       list_snp_id = chr_snps,
#       ncores      = cores
#     )

#     print(paste(chr_snps, "out of", chr_snps_length, "base file variants for chromosome", bgen_chr, "have been matched with the bgen variants and were written to", chr_backingfile))
#   }

# } else if (bgen_objects_available == "TRUE") {
#   print(paste("Using the bgen objects stored in", bgen_objects_path))

# } else {
#   print("BGEN_OBJECT_AVAILABLE parameter not recognized, exiting ...")
#   quit(status=1)
# }

# # Read reference data
# # ref <- bigreadr::fread2(reference_data)
# # CHR <- ref$chr
# # POS <- ref$BP

# # Match variants between summary statistics and genotype data
# # info_snp <- snp_match(sumstats, map)

# # Get list of parsed bgen (.rds) files
# # bgen_object_files <- list.files(path=bgen_objects_path, pattern=paste0("^", "backingfile_chr", "(.*)\\.rds$"), full.names=TRUE)
# # print("Found bgen objects:")
# # print(bgen_object_files)

# # for (bgen_object in bgen_object_files) { 
  
# #   oi <- readRDS(bgen_object)
# #   G <- oi$genotypes
# #   CHR <- as.integer(oi$map$chromosome)
# #   POS <- oi$map$physical.pos
# #   POS2 <- snp_asGeneticPos(CHR, POS, dir=prsdir, ncores=cores)

# #   chr <- as.integer(str_match(bgen_object, paste0("backingfile_chr", "(.*).rds"))[2])

# #   ind.chr <- which(CHR == chr)

# #   corr <- snp_cor(G, ind.col=ind.chr2, ncores=cores, infos.pos=POS2[ind.chr2], size=3/1000)








#   map <- oi$map[-(2:3)]
#   names(map) <- c("chr", "pos", "a0", "a1")
#   info_snp <- snp_match(sumstats, map)


#   ## indices in 'sumstats'
#   ind.chr <- which(sumstats$chr == chr)
#   ## indices in 'G'
#   ind.chr2 <- sumstats$`_NUM_ID_`[ind.chr]
#   ## indices in 'corr'
#   ind.chr3 <- match(ind.chr2, which(CHR == chr))
  
#   str(sumstats)
#   str(chr)
#   str(ind.chr)
#   str(ind.chr2)
#   str(ind.chr3)

#   # corr0 <- readRDS(paste0("data/corr/chr", chr, ".rds"))[ind.chr3, ind.chr3]
#   # corr0 <- readRDS(bgen_object)[ind.chr3, ind.chr3]

#   # if (chr == 1) {
#   if (which(bgen_object %in% bgen_object_files) == 1) {
#     df_beta <- sumstats[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")]
#     ld <- Matrix::colSums(corr0^2)
#     corr <- as_SFBM(corr0, tmp)
#   } else {
#     df_beta <- rbind(df_beta, sumstats[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")])
#     ld <- c(ld, Matrix::colSums(corr0^2))
#     corr$add_columns(corr0, nrow(corr))
#   }
# }

# str(corr)


# ## indices in info_snp
# ind.chr <- which(info_snp$chr==info_snp$chr)

# df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff")]

# ## indices in G
# ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]

# # str(ind.chr)
# # str(info_snp)
# # str(POS2)
# # str(df_beta)

# print("G")
# str(G)
# print("ind.chr2")
# str(ind.chr2)
# print("POS2")
# str(POS2[ind.chr2])
# print("ind.chr")
# str(ind.chr)

# print('1')
# ldsc <- snp_ldsc2(corr, df_beta)
# print('3')
# h2_est <- ldsc[["h2"]]
# print('4')
# corr <- as_SFBM(corr, backingfile="/hpc/dhl_ec/aligterink/ProjectFiles/SFBM_backingfile")
# print('5')
# beta_inf <- snp_ldpred2_inf(corr, df_beta, h2=h2_est)
# print('6')
# pred_inf <- big_prodVec(G, beta_inf, ind.col=ind.chr2)
# print('7')



# scores_table <- data.frame(scores=pred_inf)
# write.table(scores_table, file="/hpc/dhl_ec/aligterink/ProjectFiles/LDPRED_scores_chr21.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)





















# ## indices in info_snp
# ind.chr <- which(info_snp$chr==info_snp$chr)

# df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff")]

# ## indices in G
# ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]

# # str(ind.chr)
# # str(info_snp)
# # str(POS2)
# # str(df_beta)

# print("G")
# str(G)
# print("ind.chr2")
# str(ind.chr2)
# print("POS2")
# str(POS2[ind.chr2])
# print("ind.chr")
# str(ind.chr)

# # corr <- snp_cor(G, ind.col=ind.chr2, ncores=cores, infos.pos=POS2[ind.chr2], size=3/1000)
# print('1')
# corr <- snp_cor(G, ind.col=ind.chr2, ncores=cores, infos.pos=POS2[ind.chr2], size=3/1000)
# print('2')
# ldsc <- snp_ldsc2(corr, df_beta)
# print('3')
# h2_est <- ldsc[["h2"]]
# print('4')
# corr <- as_SFBM(corr, backingfile="/hpc/dhl_ec/aligterink/ProjectFiles/SFBM_backingfile")
# print('5')
# beta_inf <- snp_ldpred2_inf(corr, df_beta, h2=h2_est)
# print('6')
# pred_inf <- big_prodVec(G, beta_inf, ind.col=ind.chr2)
# print('7')

# str(pred_inf)
# str(h2_est)
# str(beta_inf)
# str(ind.chr2)
# str(df_beta)

# scores_table <- data.frame(scores=pred_inf)
# write.table(scores_table, file="/hpc/dhl_ec/aligterink/ProjectFiles/LDPRED_scores_chr21.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)
















# ### This is for only chromosome 21
# # Read bgen object
# obj.bigSNP <- snp_attach(paste(bgen_objects_path, ".rds", sep=""))
# map <- obj.bigSNP$map[-(2:3)]
# names(map) <- c("chr", "pos", "a0", "a1")
# map$chr <- as.numeric(map$chr)

# # Match variants between summary statistics and genotype data
# info_snp <- snp_match(sumstats, map)

# G <- obj.bigSNP$genotypes
# CHR <- map$chr
# POS <- map$pos

# POS2 <- snp_asGeneticPos(CHR, POS, dir=prsdir, ncores=cores)

# ## indices in info_snp
# ind.chr <- which(info_snp$chr == 21)

# df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff")]

# ## indices in G
# ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]

# corr <- snp_cor(G, ind.col=ind.chr2, ncores=cores, infos.pos=POS2[ind.chr2], size=3/1000)
# ldsc <- snp_ldsc2(corr, df_beta))
# h2_est <- ldsc[["h2"]]
# corr <- as_SFBM(corr, backingfile="/hpc/dhl_ec/aligterink/ProjectFiles/SFBM_backingfile")
# beta_inf <- snp_ldpred2_inf(corr, df_beta, h2=h2_est)
# pred_inf <- big_prodVec(G, beta_inf, ind.col=ind.chr2)
# str(pred_inf)

