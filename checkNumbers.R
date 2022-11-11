## Description: Check the number of mutations in various categories
# ---- Fri, Nov 11, 2022 @ 12:25 ----

library(data.table)

# Theoretical
load("/Users/deepankar/OneDrive - O365 Turun yliopisto/Git/GitHub/KE-Group/ERBB.Miscellaneous/Posssible mutants/ERBB3/ERBB3_Variants.bin")

print(paste0("Theoretical nt SNV = ", length(unique(Variants$mutants))))
print(paste0("Theoretical AA SNVs = ", length(unique(Variants$mutantprotein))))

rm(list=ls()); gc()

# Illumina
# V1
illumina_v1 <- as.data.table(readRDS("/Users/deepankar/Documents/Seafile/NGS-Data-DC/BaseSpace/20201125 Marika ERBB3 iSCREAM/Analysis/Figures/V1/20210108.ERBB3_iSCREAM_V1_Mutations.RDS"))

sub <- illumina_v1[!is.na(Plasmid.VF),]
print(paste0("Plasmid nt SNV = ", length(unique(sub$MutID))))
print(paste0("Plasmid AA SNVs = ", length(unique(sub$AAchange))))

sub <- illumina_v1[!is.na(Plasmid_PCR.VF),]
print(paste0("PCR from Plasmid nt SNV = ", length(unique(sub$MutID))))
print(paste0("PCR from Plasmid AA SNVs = ", length(unique(sub$AAchange))))

sub <- illumina_v1[!is.na(V1_WEHI.VF),]
print(paste0("WEHI nt SNV = ", length(unique(sub$MutID))))
print(paste0("WEHI AA SNVs = ", length(unique(sub$AAchange))))

sub <- illumina_v1[!is.na(V1_NRG.VF),]
print(paste0("NRG nt SNV = ", length(unique(sub$MutID))))
print(paste0("NRG AA SNVs = ", length(unique(sub$AAchange))))

sub <- illumina_v1[!is.na(V1_NRG2noIL3.VF),]
print(paste0("NRG -> NO IL3 nt SNV = ", length(unique(sub$MutID))))
print(paste0("NRG -> NO IL3 AA SNVs = ", length(unique(sub$AAchange))))

illumina_v1[,FC:= V1_NRG2noIL3.VF/V1_NRG.VF]
sub <- illumina_v1[!is.na(FC),]
print(paste0("FC v1 nt SNV = ", length(unique(sub$MutID))))
print(paste0("FC v1 AA SNVs = ", length(unique(sub$AAchange))))

# V2
rm(list=ls());gc()
illumina_v2 <- as.data.table(readRDS("/Users/deepankar/Documents/Seafile/NGS-Data-DC/BaseSpace/20201125 Marika ERBB3 iSCREAM/Analysis/Figures/V2/20201216.ERBB3_iSCREAM_V2_Mutations.RDS"))

sub <- illumina_v2[!is.na(Plasmid.VF),]
print(paste0("Plasmid nt SNV = ", length(unique(sub$MutID))))
print(paste0("Plasmid AA SNVs = ", length(unique(sub$AAchange))))

sub <- illumina_v2[!is.na(Plasmid_PCR.VF),]
print(paste0("PCR from Plasmid nt SNV = ", length(unique(sub$MutID))))
print(paste0("PCR from Plasmid AA SNVs = ", length(unique(sub$AAchange))))

sub <- illumina_v2[!is.na(V2_WEHI.VF),]
print(paste0("WEHI nt SNV = ", length(unique(sub$MutID))))
print(paste0("WEHI AA SNVs = ", length(unique(sub$AAchange))))

sub <- illumina_v2[!is.na(V2_NRG.VF),]
print(paste0("NRG nt SNV = ", length(unique(sub$MutID))))
print(paste0("NRG AA SNVs = ", length(unique(sub$AAchange))))

sub <- illumina_v2[!is.na(V2_NRG2noIL3.VF),]
print(paste0("NRG -> NO IL3 nt SNV = ", length(unique(sub$MutID))))
print(paste0("NRG -> NO IL3 AA SNVs = ", length(unique(sub$AAchange))))

sub <- illumina_v2[!is.na(V2_noIL3.VF),]
print(paste0("WEHI -> NO IL3 nt SNV = ", length(unique(sub$MutID))))
print(paste0("WEHI -> NO IL3 AA SNVs = ", length(unique(sub$AAchange))))

illumina_v2[,FC:= V2_noIL3.VF/V2_WEHI.VF]
sub <- illumina_v2[!is.na(FC),]
print(paste0("FC (WEHI) v2 nt SNV = ", length(unique(sub$MutID))))
print(paste0("FC (WEHI) v2 AA SNVs = ", length(unique(sub$AAchange))))

illumina_v2[,FC:= V2_NRG2noIL3.VF/V2_NRG.VF]
sub <- illumina_v2[!is.na(FC),]
print(paste0("FC (NRG) v2 nt SNV = ", length(unique(sub$MutID))))
print(paste0("FC (NRG) v2 AA SNVs = ", length(unique(sub$AAchange))))


# Combined "all in one"
rm(list=ls());gc()
illumina_combined <- readRDS("/Users/deepankar/Documents/Seafile/NGS-Data-DC/BaseSpace/20201125 Marika ERBB3 iSCREAM/Analysis/Figures/20220506.plotMuts.RDS")
print(paste0("FC (combined) v2 nt SNV = ", length(unique(illumina_combined$MutID))))
print(paste0("FC (combined) v2 AA SNVs = ", length(unique(illumina_combined$AAchange))))
