## Description: Combining readcounts for genotypes and calculate FC
# ---- Thu, Nov 10, 2022 @ 22:19 ----
rm(list = ls())
gc()
library(data.table)
#--------> Combining files <------------
# True Hits
setwd("/Volumes/DATA/Seafile/NGS-Data-DC/PacBio/20210909 iSCREAM Library and samples/analysis/Figures/Bootstrap/Data/")
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

# hits <- combine(dirname = "20221219_084028/")
# set1 <- combine(dirname = "20221219_083754/")

hits <- combine(dirname = "20221219_161723/")
set1 <- combine(dirname = "20221219_094048/")
control <- combine(dirname = "20221224_041520/")


# Calculating FC and taking max (as done for the scatter plot)
hits <- combine(dirname = "20221219_161723/")
hits[, FC1 := V1_NRG_to_NoLig/V1_NRG]
hits[, FC2 := V2_NRG_to_NoLig/V2_NRG]
hits[, FC3 := V2_WEHI_to_NoLig/V2_WEHI]
hits[, maxFC := max(FC1, FC2, FC3, na.rm = T), by = 1:nrow(hits)]

hits <- hits[maxFC != -Inf,.(genotypes, maxFC)]
hits$FC_without_exponent <- format(hits$maxFC,scientific = F,digits = 2)
hits[, N_mutants := stringi::stri_count_fixed(str = genotypes,
                                              pattern = "/") + 1]

# Setting Number of muts to 1, 2 and many
hits$N_mutants[hits$N_mutants>2] <- 3
hits$N_mutants <- factor(hits$N_mutants)
levels(hits$N_mutants) <- c("1","2","3+")

FC_above_25 <- hits[maxFC >= 25,]
FC_below_25 <- hits[maxFC < 25,]
table(FC_above_25$N_mutants)/nrow(FC_above_25)
table(FC_below_25$N_mutants)/nrow(FC_below_25)


# Set 1
# Calculating FC and taking max (as done for the scatter plot)
set1 <- combine(dirname = "20221219_094048/")
set1[, FC1 := V1_NRG_to_NoLig/V1_NRG]
set1[, FC2 := V2_NRG_to_NoLig/V2_NRG]
set1[, FC3 := V2_WEHI_to_NoLig/V2_WEHI]
set1[, maxFC := max(FC1, FC2, FC3, na.rm = T), by = 1:nrow(set1)]

set1 <- set1[maxFC != -Inf,.(genotypes, maxFC)]
set1$FC_without_exponent <- format(set1$maxFC,scientific = F,digits = 2)
set1[, N_mutants := stringi::stri_count_fixed(str = genotypes,
                                              pattern = "/") + 1]

# Setting Number of muts to 1, 2 and many
set1$N_mutants[set1$N_mutants>2] <- 3
set1$N_mutants <- factor(set1$N_mutants)
levels(set1$N_mutants) <- c("1","2","3+")
table(set1$N_mutants)/nrow(set1)

# Control mutations (FC 1 - 25)
# Calculating FC and taking max (as done for the scatter plot)
control <- combine(dirname = "20221224_041520/")
control[, FC1 := V1_NRG_to_NoLig/V1_NRG]
control[, FC2 := V2_NRG_to_NoLig/V2_NRG]
control[, FC3 := V2_WEHI_to_NoLig/V2_WEHI]
control[, maxFC := max(FC1, FC2, FC3, na.rm = T), by = 1:nrow(control)]

control <- control[maxFC != -Inf,.(genotypes, maxFC)]
control$FC_without_exponent <- format(control$maxFC,scientific = F,digits = 2)
control[, N_mutants := stringi::stri_count_fixed(str = genotypes,
                                              pattern = "/") + 1]

# Setting Number of muts to 1, 2 and many
control$N_mutants[control$N_mutants>2] <- 3
control$N_mutants <- factor(control$N_mutants)
levels(control$N_mutants) <- c("1","2","3+")
table(control$N_mutants)/nrow(control)

write.table(
  x = control,
  file = "~/Desktop/control.tsv",
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = T
)

write.table(
  x = FC_above_25,
  file = "~/Desktop/FC_above_25.tsv",
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = T
)

write.table(
  x = FC_below_25,
  file = "~/Desktop/FC_below_25.tsv",
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = T
)


# Summary
message("FC >= 25")
table(FC_above_25$N_mutants)/nrow(FC_above_25)
message("FC < 25")
table(FC_below_25$N_mutants)/nrow(FC_below_25)
message("Control")
table(control$N_mutants)/nrow(control)
 
# #------ Attempting plots ---------
# hits <- data.table::melt.data.table(hits, id.vars="genotypes")
# hits <- hits[-which(is.na(hits$value)),]
# hits[, N_mutants := stringi::stri_count_fixed(str = genotypes,
#                                               pattern = "/") + 1]
# 
# set1 <- data.table::melt.data.table(set1, id.vars="genotypes")
# set1 <- set1[-which(is.na(set1$value)),]
# set1[, N_mutants := stringi::stri_count_fixed(str = genotypes,
#                                               pattern = "/") + 1]
# 
# library(ggplot2)
# source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
# customtheme <- DC_theme_generator(type = "L")
# 
# ggplot(data = set1,
#        aes(x = N_mutants, y=value, color = N_mutants))+
#   facet_wrap(.~variable)+
#   geom_point()+
#   coord_trans(y="log10")+
#   ggtitle("Selected controls")+
#   xlab("Number of mutations in cDNA")+
#   ylab("Allele frequency")+
#   viridis::scale_color_viridis(direction = -1,end = 0.9,option = "magma")+
#   customtheme
# 
# ggplot(data = hits,
#        aes(x = N_mutants, y=value, color = N_mutants))+
#   facet_wrap(.~variable)+
#   geom_point()+
#   coord_trans(y="log10")+
#   ggtitle("Hits from the screen")+
#   xlab("Number of mutations in cDNA")+
#   ylab("Allele frequency")+
#   scale_y_continuous(breaks = c(1,0.5,0.25,0.1,0.01,0.001, 0.001))+
#   viridis::scale_color_viridis(direction = -1,end = 0.9,option = "magma")+
#   customtheme
# 
# #-------------
# 
# hits <- na.exclude(combine(dirname = "20221219_084028/"))
# hits$FC <- hits$NRG_to_NoLig/hits$NRG
# 
# # Setting for 1, 2 or multiple mutants. Multiple = 3
# hits$N_mutants[hits$N_mutants>2] = 3
# set1$N_mutants[set1$N_mutants>2] = 3
# 
