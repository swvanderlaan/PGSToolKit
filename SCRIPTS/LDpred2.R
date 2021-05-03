#!/hpc/local/CentOS7/dhl_ec/software/R-3.6.3/R-3.6.3/bin/Rscript --vanilla

# Required packages
packages = c("optparse",
            "data.table",
            "remotes",
            "bigsnpr",
            "doParallel",
            "RSQLite",
            "dbplyr")

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

#parameters

rds <- bigsnpr::snp_readBGEN(
  bgenfiles   = c("/hpc/dhl_ec/aligterink/ProjectFiles/chr1_validation_1000lines_8bit.bgen"),
  list_snp_id = list(c("1_10177_A_AC", "1_10235_T_TA", "1_10352_T_TA", "1_10539_C_A", "1_10642_G_A")),
  backingfile = "/hpc/dhl_ec/aligterink/ProjectFiles/tmp/ldpredtmp/backingfile",
  ncores      = 1
)

ukb <- snp_attach(rds)

G <- ukb$genotypes
CHR <- as.numeric(ukb$map$chromosome)
POS <- ukb$map$physical.pos

print(G)
print(CHR)
print(POS)


###### Read bgen files

# map_hapmap3 <- bigreadr::fread2("ldblk_1kg_eur/snpinfo_1kg_hm3")

# cl <- makeCluster(22)
# parallel::clusterExport(cl, "map_hapmap3")
# list_snp_id <- parLapply(cl, 1:22, function(chr) {
#   mfi <- paste0("UKBB/mfi/ukb_mfi_chr", chr, "_v3.txt")
#   infos_chr <- bigreadr::fread2(mfi, showProgress = FALSE)
#   joined <- dplyr::inner_join(cbind(chr = chr, infos_chr), map_hapmap3[1:2],
#                               by = c("chr" = "CHR", "V2" = "SNP"))
#   with(joined[!vctrs::vec_duplicate_detect(joined$V2), ],
#        paste(chr, V3, V4, V5, sep = "_"))
# })
# stopCluster(cl)

# sum(lengths(list_snp_id))  # 1,117,493


# sample <- bigreadr::fread2("UKBB/ukb58024_imp_chr1_v3_s487296.sample")
# str(sample)
# sample <- sample[-1, ]


# csv <- "UKBB/ukb41181.csv"
# df0 <- bigreadr::fread2(
#   csv,
#   select = c("eid", "22020-0.0", paste0("22009-0.", 1:16)),
#   col.names = c("eid", "used_in_pca", paste0("PC", 1:16))
# )

# # not removed & quality controled / unrelated & white-British
# ind.indiv <- match(df0$eid, sample$ID_2)
# sub <- which(!is.na(ind.indiv) & df0$used_in_pca)

# # Genetically homogeneous
# dist <- bigutilsr::dist_ogk(as.matrix(df0[sub, -(1:2)]))
# sub2 <- sub[log(dist) < 5]
# length(sub2) # 362,320
# saveRDS(sub2, "data/ind_sub_csv.rds")


# (NCORES <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) - 1L)

# system.time(
#   rds <- bigsnpr::snp_readBGEN(
#     bgenfiles   = glue::glue("UKBB/bgen/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
#     list_snp_id = list_snp_id,
#     backingfile = "data/UKBB_imp_HM3",
#     ind_row     = ind.indiv[sub2],
#     ncores      = NCORES
#   )
# ) # 43 min with 23 cores

# library(bigsnpr)
# ukb <- snp_attach("data/UKBB_imp_HM3.rds")
# G <- ukb$genotypes
# CHR <- as.numeric(ukb$map$chromosome)
# POS <- ukb$map$physical.pos
# dim(G) # 362,320 x 1,117,493
# file.size(G$backingfile) / 1024^3  # 377 GB


####### Actual calculation

# # Read in the summary statistic file
# sumstats <- bigreadr::fread2("Height.QC.gz") 

# # LDpred 2 requires the header to follow the exact naming
# names(sumstats) <-
#     c("chr",
#     "pos",
#     "rsid",
#     "a1",
#     "a0",
#     "n_eff",
#     "beta_se",
#     "p",
#     "OR",
#     "INFO",
#     "MAF")

# # Transform the OR into log(OR)
# sumstats$beta <- log(sumstats$OR)

# # Filter out hapmap SNPs
# sumstats <- sumstats[sumstats$rsid%in% info$rsid,]


# # Get maximum amount of cores
# NCORES <- nb_cores()

# # Open a temporary file
# tmp <- tempfile(tmpdir = "tmp-data")
# on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

# # Initialize variables for storing the LD score and LD matrix
# corr <- NULL
# ld <- NULL

# # We want to know the ordering of samples in the bed file 
# info_snp <- NULL
# fam.order <- NULL
# for (chr in 1:22) {
  
#     # preprocess the bed file (only need to do once for each data set)
#     # Assuming the file naming is EUR_chr#.bed

#     snp_readBed(paste0("EUR_chr",chr,".bed"))
#     # now attach the genotype object
#     obj.bigSNP <- snp_attach(paste0("EUR_chr",chr,".rds"))
    
#     # extract the SNP information from the genotype
#     map <- obj.bigSNP$map[-3]
#     names(map) <- c("chr", "rsid", "pos", "a1", "a0")

#     # perform SNP matching
#     tmp_snp <- snp_match(sumstats[sumstats$chr==chr,], map)
#     info_snp <- rbind(info_snp, tmp_snp)

#     # Assign the genotype to a variable for easier downstream analysis
#     genotype <- obj.bigSNP$genotypes

#     # Rename the data structures
#     CHR <- map$chr
#     POS <- map$pos

#     # get the CM information from 1000 Genome
#     # will download the 1000G file to the current directory (".")
#     POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")

#     # calculate LD
#     # Extract SNPs that are included in the chromosome
#     ind.chr <- which(tmp_snp$chr == chr)
#     ind.chr2 <- tmp_snp$`_NUM_ID_`[ind.chr]

#     # Calculate the LD
#     corr0 <- snp_cor(
#             genotype,
#             ind.col = ind.chr2,
#             ncores = NCORES,
#             infos.pos = POS2[ind.chr2],
#             size = 3 / 1000
#         )
#     if (chr == 1) {
#         ld <- Matrix::colSums(corr0^2)
#         corr <- as_SFBM(corr0, tmp)
#     } else {
#         ld <- c(ld, Matrix::colSums(corr0^2))
#         corr$add_columns(corr0, nrow(corr))
#     }

#     # We assume the fam order is the same across different chromosomes
#     if(is.null(fam.order)){
#         fam.order <- as.data.table(obj.bigSNP$fam)
#     }
# }

# # Rename fam order
# setnames(fam.order,
#         c("family.ID", "sample.ID"),
#         c("FID", "IID"))






