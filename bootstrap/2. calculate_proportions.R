## Description: Combining readcounts for genotypes and calculate FC
# ---- Thu, Nov 10, 2022 @ 22:19 ----
rm(list = ls())
gc()
library(data.table)
#--------> Combining files <------------
# True Hits
setwd("/Volumes/DATA/Seafile/NGS-Data-DC/PacBio/20210909 iSCREAM Run1/analysis/Figures/Bootstrap/Data/")
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

hits <- combine(dirname = "20221219_084028/")
# set1 <- combine(dirname = "20221219_083754/")
set1 <- combine(dirname = "20221219_094048/")


hits <- data.table::melt.data.table(hits, id.vars="genotypes")
hits <- hits[-which(is.na(hits$value)),]
hits[, N_mutants := stringi::stri_count_fixed(str = genotypes,
                                              pattern = "/") + 1]

set1 <- data.table::melt.data.table(set1, id.vars="genotypes")
set1 <- set1[-which(is.na(set1$value)),]
set1[, N_mutants := stringi::stri_count_fixed(str = genotypes,
                                              pattern = "/") + 1]

library(ggplot2)
source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type = "L")

ggplot(data = set1,
       aes(x = N_mutants, y=value, color = N_mutants))+
  facet_wrap(.~variable)+
  geom_point()+
  coord_trans(y="log10")+
  ggtitle("Selected controls")+
  xlab("Number of mutations in cDNA")+
  ylab("Allele frequency")+
  viridis::scale_color_viridis(direction = -1,end = 0.9,option = "magma")+
  customtheme

ggplot(data = hits,
       aes(x = N_mutants, y=value, color = N_mutants))+
  facet_wrap(.~variable)+
  geom_point()+
  coord_trans(y="log10")+
  ggtitle("Hits from the screen")+
  xlab("Number of mutations in cDNA")+
  ylab("Allele frequency")+
  scale_y_continuous(breaks = c(1,0.5,0.25,0.1,0.01,0.001, 0.001))+
  viridis::scale_color_viridis(direction = -1,end = 0.9,option = "magma")+
  customtheme

#-------------

hits <- na.exclude(combine(dirname = "20221219_084028/"))
hits$FC <- hits$NRG_to_NoLig/hits$NRG

# Setting for 1, 2 or multiple mutants. Multiple = 3
hits$N_mutants[hits$N_mutants>2] = 3
set1$N_mutants[set1$N_mutants>2] = 3

