rm(list=ls())
source("/Users/deepankar/OneDrive - O365 Turun yliopisto/Git/GitHub/KE-Group/ERBB3_iSCREAM/bootstrap/0. setup environ for search.R")

# -----------> RUNNING <-----------
# -----------> Running with Marika's selected hits <-----------
# ERBB3 iSCREAM hits

# Setting the hits from Marika's figure to search the data.
hits <- readRDS("/Users/deepankar/Documents/Seafile/NGS-Data-DC/BaseSpace/20201125 Marika ERBB3 iSCREAM/Analysis/Figures/20220506.plotMuts.RDS")
hits <- hits[hits$FC > 25]

# Finding reference altered nucleotides for the selected mutations.
ref_start <- stringi::stri_locate_first_regex(str = hits$MutID,pattern = "chr[0-9]{1,2}:[0-9]+")[,"end"]
ref_end <- stringi::stri_locate_first_fixed(str = hits$MutID,pattern = ">")[,"end"]

hits[,Ref_nt := substring(hits$MutID,ref_start+1,ref_end-1)]
hits[,Alt_nt := substring(hits$MutID,ref_end+1)]

rm(ref_start,ref_end)

# Find mutIDs to to explore in the PacBio data
chr_mutation_vector <- find_mutID_to_test(hits)

data_analysis(chr_mutation_vector)


# -----------> 18 random hits FC≈1 (1±10%) <-----------

# Setting the hits from Marika's figure to search the data.
hits <- readRDS("/Users/deepankar/Documents/Seafile/NGS-Data-DC/BaseSpace/20201125 Marika ERBB3 iSCREAM/Analysis/Figures/20220506.plotMuts.RDS")

# Taking 18 samples with FC of 1 ± 10% (i.e. between 0.9 and 1)
idx <- sample(
  x = which(hits$FC >= 0.2 & hits$FC <= 5 & hits$VF > 0.05),
  size = 18,
  replace = F
)
hits <- hits[idx,]

# Finding reference altered nucleotides for the selected mutations.
ref_start <- stringi::stri_locate_first_regex(str = hits$MutID,pattern = "chr[0-9]{1,2}:[0-9]+")[,"end"]
ref_end <- stringi::stri_locate_first_fixed(str = hits$MutID,pattern = ">")[,"end"]

hits[,Ref_nt := substring(hits$MutID,ref_start+1,ref_end-1)]
hits[,Alt_nt := substring(hits$MutID,ref_end+1)]

rm(ref_start,ref_end)

# Find mutIDs to to explore in the PacBio data
chr_mutation_vector <- find_mutID_to_test(hits)

data_analysis(chr_mutation_vector)

