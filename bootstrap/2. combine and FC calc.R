## Description: Combining readcounts for genotypes and calculate FC
# ---- Thu, Nov 10, 2022 @ 22:19 ----
rm(list = ls())
gc()
library(data.table)
#--------> Combining files <------------
rm(list=ls())
# True Hits
setwd("/Volumes/DATA/Seafile/NGS-Data-DC/PacBio/20210909 iSCREAM Run1/analysis/Figures/Bootstrap/Data/")
dir.create("./combined_hits", showWarnings = F)

combine <- function(dirname){
  setwd(dirname)
  data.files <- list.files(".",pattern = "V1")
  
  dat <- data.table()
  for(FILE in data.files){
    dat <- rbind(dat,fread(FILE))
  }
  
  genotypes <- unique(dat$muts)
  dat <- data.table(genotypes,
                    WEHI = NA,
                    NRG = NA,
                    NRG_to_NoLig = NA)
  
  # NRG_to_NoLig
  tmp <-
    fread(data.files[grep(pattern = "NRG_to_NoLig", 
                          x = data.files)])
  idx <- match(x = dat$genotypes,table = tmp$muts)
  dat$NRG_to_NoLig <- tmp$N[idx]
  
  tmp <- fread(data.files[grep(pattern = "NRG_[^to]", 
                               x = data.files)])
  idx <- match(x = dat$genotypes,table = tmp$muts)
  dat$NRG <- tmp$N[idx]
  
  tmp <- fread(data.files[grep(pattern = "WEHI", 
                               x = data.files)])
  idx <- match(x = dat$genotypes,table = tmp$muts)
  dat$WEHI <- tmp$N[idx]
  rm(tmp)
  setorder(dat, -NRG_to_NoLig,na.last = T)
  
  
  dat[,NRG_to_NoLig := NRG_to_NoLig / 291190] # total reads (from flagstat)
  dat[,NRG := NRG / 217694] # total reads (from flagstat)
  dat[,WEHI := WEHI / 144629]   # total reads (from flagstat)
  
  dat[, FC := NRG_to_NoLig/NRG]
  dat <- dat[!is.na(FC),]
  dat[, FC_text := format(FC, scientific = F,digits = 2)]
  setorder(dat, -FC, na.last = T)
  setwd("..")
  return(dat)
}

hits <- combine(dirname = "20221110_214857/")
set1 <- combine(dirname = "20221110_220920/")

