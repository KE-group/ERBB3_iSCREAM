## Description: Combining readcounts and calculating FC for 
##              bootstrap of controls (10 rounds x 18 mutations)
# ---- Fri, Feb 03, 2023 @ 07:39 ----

rm(list = ls())
gc()
library(data.table)
#--------> Combining files <------------
setwd("/Volumes/UTU_KE_2TB/Seafile/NGS-Data-DC/PacBio/20210909 iSCREAM Library and samples/analysis/Figures/Bootstrap/Data/")
dir.create("./combined_hits", showWarnings = F)

combine <- function(dirname){
  setwd(dirname)
  data.files <- list.files(".",pattern = "csv")
  
  # Finding unique genotypes
  dat <- data.table()
  for(FILE in data.files){
    dat <- rbind(dat,fread(FILE))
  }
  genotypes <- unique(dat$muts)
  
  # Initializing data table for the genotypes vs Time-points
  dat <- data.table(genotypes,
                    V1_WEHI = NA,
                    V1_NRG = NA,
                    V1_NRG_to_NoLig = NA,
                    V2_WEHI = NA,
                    V2_NRG = NA,
                    V2_NRG_to_NoLig = NA,
                    V2_WEHI_to_NoLig = NA)
  
  # V1
  tmp <- fread("ERBB3_NRG_to_NoLig_V1.csv")
  idx <- match(x = dat$genotypes, table = tmp$muts)
  dat$V1_NRG_to_NoLig <- tmp$N[idx]
  
  tmp <- fread("ERBB3_NRG_V1.csv")
  idx <- match(x = dat$genotypes,table = tmp$muts)
  dat$V1_NRG <- tmp$N[idx]
  
  tmp <- fread("ERBB3_WEHI_V1.csv")
  idx <- match(x = dat$genotypes,table = tmp$muts)
  dat$V1_WEHI <- tmp$N[idx]
  rm(tmp)
  
  # V1
  tmp <- fread("ERBB3_NoLigand_V2.csv")
  idx <- match(x = dat$genotypes, table = tmp$muts)
  dat$V2_WEHI_to_NoLig <- tmp$N[idx]
  
  tmp <- fread("ERBB3_NRG_to_NoLig_V2.csv")
  idx <- match(x = dat$genotypes, table = tmp$muts)
  dat$V2_NRG_to_NoLig <- tmp$N[idx]
  
  tmp <- fread("ERBB3_NRG_V2.csv")
  idx <- match(x = dat$genotypes,table = tmp$muts)
  dat$V2_NRG <- tmp$N[idx]
  
  tmp <- fread("ERBB3_WEHI_V2.csv")
  idx <- match(x = dat$genotypes,table = tmp$muts)
  dat$V2_WEHI <- tmp$N[idx]
  rm(tmp)
  
  setorder(dat, -V2_WEHI_to_NoLig, na.last = T)
  
  # convert to Allele Fractions by dividing with total read count in that sample
  # The divisor is coming from the output of samtools flagstat stored next to the BAM file of the corresponding sample in the directory
  # "Volumes/DATA/Seafile/NGS-Data-DC/PacBio/20210909 iSCREAM Run1/analysis/BAM"
  
  dat[,V1_NRG_to_NoLig := V1_NRG_to_NoLig / 291190] # total reads (from flagstat)
  dat[,V1_NRG := V1_NRG / 217694] # total reads (from flagstat)
  dat[,V1_WEHI := V1_WEHI / 144629]   # total reads (from flagstat)
  
  dat[,V2_NRG_to_NoLig := V2_NRG_to_NoLig / 230103] # total reads (from flagstat)
  dat[,V2_NRG := V2_NRG / 246725] # total reads (from flagstat)
  dat[,V2_WEHI := V2_WEHI / 99808]   # total reads (from flagstat)
  dat[,V2_WEHI_to_NoLig := V2_WEHI_to_NoLig / 182885]  # total reads (from flagstat)
  setwd("..") # important as it returns to the parent directory
  return(setDT(dat))
}

# Control mutations (FC 1 - 25)
calculateSummary <- function(dirname) {
  ## Description:   This function runs the combine genotypes across samples
  ##                and computes the number of genotypes with 1, 2 or 3+ mutations
  control <- combine(dirname = dirname)
  control[, FC1 := V1_NRG_to_NoLig/V1_NRG]
  control[, FC2 := V2_NRG_to_NoLig/V2_NRG]
  control[, FC3 := V2_WEHI_to_NoLig/V2_WEHI]
  control[, maxFC := max(FC1, FC2, FC3, na.rm = T), by = 1:nrow(control)]
  
  control <- control[maxFC != -Inf,.(genotypes, maxFC)]
  control[, N_mutants := stringi::stri_count_fixed(str = genotypes,
                                                   pattern = "/") + 1]
  
  # Setting Number of muts to 1, 2 and many
  control$N_mutants[control$N_mutants>2] <- 3
  control$N_mutants <- factor(control$N_mutants)
  levels(control$N_mutants) <- c("1","2","3+")
  return(table(control$N_mutants)/nrow(control))
  
}
library(parallel)

output <- parallel::mclapply(
  X = dir(path = "/Volumes/UTU_KE_2TB/Seafile/NGS-Data-DC/PacBio/20210909 iSCREAM Library and samples/analysis/Figures/Bootstrap/Data", pattern = "20230202", recursive = F),
  FUN = function(X)
    calculateSummary(dirname = X),
  mc.cores = parallel::detectCores()
)

unlist(lapply(output, `[[`, 3))

resultsDT <-
  data.table(
    Replicate = paste0("control-run-", seq(1:10)),
    OneMut = unlist(lapply(output, `[[`, 1)),
    TwoMut = unlist(lapply(output, `[[`, 2)),
    More = unlist(lapply(output, `[[`, 3))
  )

# write.table(
#   x = control,
#   file = "~/Desktop/control.tsv",
#   quote = F,
#   sep = "\t",
#   row.names = F,
#   col.names = T
# )
