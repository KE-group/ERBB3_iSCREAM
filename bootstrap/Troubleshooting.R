source('/Users/deepankar/OneDrive - O365 Turun yliopisto/Git/GitHub/KE-Group/ERBB3_iSCREAM/bootstrap/0. setup environ for search.R')
Sys.setlocale("LC_MESSAGES", 'en_GB.UTF-8')
Sys.setenv(LANG = "en_US.UTF-8")

# Setting the hits from Marika's figure to search the data.
hits <- readRDS("/Users/deepankar/Documents/Seafile/NGS-Data-DC/Illumina/iSCREAM/ERBB3/20201125 Marika ERBB3 iSCREAM/Analysis/Figures/20220506.plotMuts.RDS")
# hits <- hits[hits$FC > 25,]
# hits <- hits[hits$FC > 25 & hits$VF>0.01, ]
 
# # Selecting 18 random mutations with FC < 25 and FC > 1 and VF > 1%
# hits <- hits[sample(x = which(hits$FC < 25 & FC > 1 & hits$VF>0.01),size = 18,replace = F), ]

# Selecting all mutations between FC of 1 - 25 (in Illumina data!!)
hits <- hits[which(hits$FC < 25 & hits$FC > 1 & hits$VF>0.01),]
hits <- hits[-grep(pattern = "*",x =  hits$AAchange, fixed = T),] # removing rows with nonsense mutations

# Finding reference altered nucleotides for the selected mutations.
ref_start <- stringi::stri_locate_first_regex(str = hits$MutID,pattern = "chr[0-9]{1,2}:[0-9]+")[,"end"]
ref_end <- stringi::stri_locate_first_fixed(str = hits$MutID,pattern = ">")[,"end"]

hits[,Ref_nt := substring(hits$MutID,ref_start+1,ref_end-1)]
hits[,Alt_nt := substring(hits$MutID,ref_end+1)]

rm(ref_start,ref_end)

# Find mutIDs to to explore in the PacBio data
chr_mutation_vector <- find_mutID_to_test(hits)

data_analysis(chr_mutation_vector)
