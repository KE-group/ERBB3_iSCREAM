## Description: Script for making Excel tables for Marika's ERBB3 screen
# ---- Sat, Nov 19, 2022 @ 08:33 ----

rm(list = ls())
gc()
library(data.table)

V1 <- readRDS("/Volumes/DATA/Seafile/NGS-Data-DC/PacBio/20210909 iSCREAM Run1/analysis/filtered/R.Binary/ERBB3_V1_iSCREAM.RDS")

V2 <- readRDS("/Volumes/DATA/Seafile/NGS-Data-DC/PacBio/20210909 iSCREAM Run1/analysis/filtered/R.Binary/ERBB3_V2_iSCREAM.RDS")


# creating unique list of mutations in the screen
excel <- unique.data.frame(rbind(V1[,c(1:3)], V2[,c(1:3)]))

# Finding reference altered nucleotides for the selected mutations.
ref_start <- stringi::stri_locate_first_regex(str = excel$MutIDs, pattern = "c\\.")[,"end"]
ref_end <- stringi::stri_locate_first_regex(str = excel$MutIDs,pattern = "[0-9]+>")[,"start"]

alt_start <- stringi::stri_locate_first_fixed(str = excel$MutIDs, pattern = ">")[,"end"]

excel[, Ref_nt := substring(MutIDs, ref_start + 1, ref_end - 1)]
excel[, Alt_nt := substring(MutIDs, alt_start+1)]

rm(ref_start,ref_end,alt_start)

excel[ , AA_change := gsub("p.","",AA_change, fixed = T)]
excel <- excel[,.(gene, MutIDs, Ref_nt,Alt_nt, AA_change)]

excel$Plasmid_PCR.AD <- V1$Plasmid_pcr_V1_V2_AD[match(x = excel$MutID, table = V1$MutID)]
excel$Plasmid_PCR.DP <- V1$Plasmid_pcr_V1_V2_DP[match(x = excel$MutID, table = V1$MutID)]
excel$Plasmid_PCR.VF <- V1$Plasmid_pcr_V1_V2_VF[match(x = excel$MutID, table = V1$MutID)]

excel$V1_WEHI.AD <- V1$WEHI_V1_AD[match(x = excel$MutID, table = V1$MutID)]
excel$V1_WEHI.DP <- V1$WEHI_V1_DP[match(x = excel$MutID, table = V1$MutID)]
excel$V1_WEHI.VF <- V1$WEHI_V1_VF[match(x = excel$MutID, table = V1$MutID)]

excel$V1_NRG.AD <- V1$NRG_V1_AD[match(x = excel$MutID, table = V1$MutID)]
excel$V1_NRG.DP <- V1$NRG_V1_DP[match(x = excel$MutID, table = V1$MutID)]
excel$V1_NRG.VF <- V1$NRG_V1_VF[match(x = excel$MutID, table = V1$MutID)]

excel$V1_NRG_to_noIL3.AD <- V1$NRG_to_NoLig_V1_AD[match(x = excel$MutID, table = V1$MutID)]
excel$V1_NRG_to_noIL3.DP <- V1$NRG_to_NoLig_V1_DP[match(x = excel$MutID, table = V1$MutID)]
excel$V1_NRG_to_noIL3.VF <- V1$NRG_to_NoLig_V1_VF[match(x = excel$MutID, table = V1$MutID)]

# excel$V1_FC <- excel$V1_NRG_to_noIL3.VF/excel$V1_NRG.VF

excel$V2_WEHI.AD <- V2$WEHI_V2_AD[match(x = excel$MutID, table = V2$MutID)]
excel$V2_WEHI.DP <- V2$WEHI_V2_DP[match(x = excel$MutID, table = V2$MutID)]
excel$V2_WEHI.VF <- V2$WEHI_V2_VF[match(x = excel$MutID, table = V2$MutID)]

excel$V2_NRG.AD <- V2$NRG_V2_AD[match(x = excel$MutID, table = V2$MutID)]
excel$V2_NRG.DP <- V2$NRG_V2_DP[match(x = excel$MutID, table = V2$MutID)]
excel$V2_NRG.VF <- V2$NRG_V2_VF[match(x = excel$MutID, table = V2$MutID)]

excel$V2_NRG_to_noIL3.AD <- V2$NRG_to_NoLig_V2_AD[match(x = excel$MutID, table = V2$MutID)]
excel$V2_NRG_to_noIL3.DP <- V2$NRG_to_NoLig_V2_DP[match(x = excel$MutID, table = V2$MutID)]
excel$V2_NRG_to_noIL3.VF <- V2$NRG_to_NoLig_V2_VF[match(x = excel$MutID, table = V2$MutID)]

excel$V2_noIL3.AD <- V2$NoLigand_V2_AD[match(x = excel$MutID, table = V2$MutID)]
excel$V2_noIL3.DP <- V2$NoLigand_V2_DP[match(x = excel$MutID, table = V2$MutID)]
excel$V2_noIL3.VF <- V2$NoLigand_V2_VF[match(x = excel$MutID, table = V2$MutID)]

# excel$V2_FC <- excel$V2_NRG_to_noIL3.VF/excel$V2_NRG.VF

write.table(
  x = excel,
  file = "~/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Marika ERBB3 screen/Figures/Processed PacBio data.tsv",
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = T
)   
