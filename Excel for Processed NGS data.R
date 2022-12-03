## Description: Script for making Excel tables for Marika's ERBB3 screen
# ---- Sat, Nov 19, 2022 @ 08:33 ----

rm(list = ls())
gc()
library(data.table)

V1 <- setDT(readRDS("/Volumes/DATA/Seafile/NGS-Data-DC/BaseSpace/20201125 Marika ERBB3 iSCREAM/Analysis/Figures/V1/20210108.ERBB3_iSCREAM_V1_Mutations.RDS"))

V2 <- setDT(readRDS("/Volumes/DATA/Seafile/NGS-Data-DC/BaseSpace/20201125 Marika ERBB3 iSCREAM/Analysis/Figures/V2/20201216.ERBB3_iSCREAM_V2_Mutations.RDS"))

excel <- unique.data.frame(rbind(V1[,c(10,1:5,26,9)], V2[,c(10,1:5,29,9)]))

excel$Plasmid_PCR.AD <- V1$Plasmid_PCR.AD[match(x = excel$MutID, table = V1$MutID)]
excel$Plasmid_PCR.DP <- V1$Plasmid_PCR.DP[match(x = excel$MutID, table = V1$MutID)]
excel$Plasmid_PCR.VF <- V1$Plasmid_PCR.VF[match(x = excel$MutID, table = V1$MutID)]

excel$V1_WEHI.AD <- V1$V1_WEHI.AD[match(x = excel$MutID, table = V1$MutID)]
excel$V1_WEHI.DP <- V1$V1_WEHI.DP[match(x = excel$MutID, table = V1$MutID)]
excel$V1_WEHI.VF <- V1$V1_WEHI.VF[match(x = excel$MutID, table = V1$MutID)]

excel$V1_NRG.AD <- V1$V1_NRG.AD[match(x = excel$MutID, table = V1$MutID)]
excel$V1_NRG.DP <- V1$V1_NRG.DP[match(x = excel$MutID, table = V1$MutID)]
excel$V1_NRG.VF <- V1$V1_NRG.VF[match(x = excel$MutID, table = V1$MutID)]

excel$V1_NRG_to_noIL3.AD <- V1$V1_NRG2noIL3.AD[match(x = excel$MutID, table = V1$MutID)]
excel$V1_NRG_to_noIL3.DP <- V1$V1_NRG2noIL3.DP[match(x = excel$MutID, table = V1$MutID)]
excel$V1_NRG_to_noIL3.VF <- V1$V1_NRG2noIL3.VF[match(x = excel$MutID, table = V1$MutID)]

# excel$V1_FC <- excel$V1_NRG_to_noIL3.VF/excel$V1_NRG.VF

excel$V2_WEHI.AD <- V2$V2_WEHI.AD[match(x = excel$MutID, table = V2$MutID)]
excel$V2_WEHI.DP <- V2$V2_WEHI.DP[match(x = excel$MutID, table = V2$MutID)]
excel$V2_WEHI.VF <- V2$V2_WEHI.VF[match(x = excel$MutID, table = V2$MutID)]

excel$V2_NRG.AD <- V2$V2_NRG.AD[match(x = excel$MutID, table = V2$MutID)]
excel$V2_NRG.DP <- V2$V2_NRG.DP[match(x = excel$MutID, table = V2$MutID)]
excel$V2_NRG.VF <- V2$V2_NRG.VF[match(x = excel$MutID, table = V2$MutID)]

excel$V2_NRG_to_noIL3.AD <- V2$V2_NRG2noIL3.AD[match(x = excel$MutID, table = V2$MutID)]
excel$V2_NRG_to_noIL3.DP <- V2$V2_NRG2noIL3.DP[match(x = excel$MutID, table = V2$MutID)]
excel$V2_NRG_to_noIL3.VF <- V2$V2_NRG2noIL3.VF[match(x = excel$MutID, table = V2$MutID)]

excel$V2_noIL3.AD <- V2$V2_noIL3.AD[match(x = excel$MutID, table = V2$MutID)]
excel$V2_noIL3.DP <- V2$V2_noIL3.DP[match(x = excel$MutID, table = V2$MutID)]
excel$V2_noIL3.VF <- V2$V2_noIL3.VF[match(x = excel$MutID, table = V2$MutID)]

# excel$V2_FC <- excel$V2_NRG_to_noIL3.VF/excel$V2_NRG.VF

write.table(
  x = excel,
  file = "OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Marika ERBB3 screen/Figures/Excel for Processed NGS data.tsv",
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = T
)   
