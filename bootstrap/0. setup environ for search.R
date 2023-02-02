## Description: The code needs list of mutations and finds them in 
##              the PacBio data and calculates their fold changes.
# ---- Wed, Nov 09, 2022 @ 11:32 ----

rm(list=ls())

# Set up code
# source("/Users/deepankar/OneDrive - O365 Turun yliopisto/Git/GitHub/KE-Group/PacBio iSCREAM run1/annotate_cDNAchange_vectorized.R")
source("/Users/chakrobd/Library/CloudStorage/GoogleDrive-robocopalpha404@gmail.com/My Drive/Git/GitHub/KE-Group/PacBio iSCREAM run1/annotate_cDNAchange_vectorized.R")
source("https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/parse_IUPAC_AAchange.R")
library(data.table)
# setwd("/Volumes/DATA/Seafile/NGS-Data-DC/PacBio/20210909 iSCREAM Library and samples/analysis")
setwd("/Volumes/UTU_KE_2TB/Seafile/NGS-Data-DC/PacBio/20210909 iSCREAM Library and samples/analysis")
dir.create("./Figures/Bootstrap/Data",recursive = T, showWarnings = F)
setwd("./Figures/Bootstrap/Data/")

all_samples <- c("1.EGFR Library",
                 "2.EGFR_WEHI",
                 "3.EGFR_NoLigand",
                 "4.ERBB2_HerceptinResistant",
                 "5.ERBB4_Plasmid_V2",
                 "6.ERBB4_NRG_D10_V2",
                 "7.ERBB3_Plasmid_pcr",
                 "8.ERBB3_WEHI_V1",
                 "9.ERBB3_NRG_V1",
                 "10.ERBB3_NRG_to_NoLig_V1",
                 "1.ERBB2 Library",
                 "2.ERBB2 WEHI",
                 "3.ERBB2 NoLigand",
                 "4.ERBB4 Lib Plasmid V1",
                 "5.ERBB4 WEHI V1",
                 "6.ERBB4 NRG D8 V1",
                 "7.ERBB3 WEHI V2",
                 "8.ERBB3 NRG V2",
                 "9.ERBB3_NRG_to_NoLig_V2",
                 "10.ERBB3 NoLigand_V2")
all_samples <- gsub(" ","_",all_samples)
searchSamples <- all_samples[c(7, 8, 9, 10,17,18,19,20)]
rm(all_samples)

# Setting up WT and codonDF for doing translations of detected nt changes (from PacBio data)
referenceSeq <- "/Volumes/UTU_KE_2TB/Seafile/NGS-Data-DC/PacBio/20210909 iSCREAM Library and samples/reference/pBG-ERBB3/ERBB3.fa"
wt_dna <- as.character(Biostrings::readDNAStringSet(filepath = referenceSeq, format = "fasta"))
codonDF <- as.data.table(DNA2codons(inputText = wt_dna)[[2]])
rm(referenceSeq)

# -----------> Define Functions <-------------
data_analysis <- function(chr_mutation_vector){
  ## Description: Runs the analysis for a given list of mutations sample by sample
  
  # Define a function within function for resolving variable scope issues
  plot_reads <- function(Mutation){
    ## Description: For a Given mutation, it finds the "genotypes" of all co-occuring mutations
    cDNAchange <- IUPAC_cDNA_Parse(Mutation)
    Mutation_AA <-
      gsub("p.",
           "",
           annotate_cDNAchange(Mutation),
           fixed = T)
    
    # mutant read names with the exact AA change
    mutant_reads <- unique(Mut[ pos==cDNAchange$Pos & ref == cDNAchange$Ref & alt == cDNAchange$Alt, .(qname)])
    
    if(dim(mutant_reads)[1] > 0 ){
      # Find other mutations on same reads
      reads <- Mut[ qname %in% mutant_reads$qname, ]
      reads[, IUPAC_cDNA := paste0("c.",ref,pos,">",alt)]
      reads[, c("ref","pos","alt") := NULL]
      # reads <- reads[1:min(nrow(reads),1000)] # for troubleshooting. Reduces runtime
      
      Mutant <- unique(reads[,.(IUPAC_cDNA)])
      
      # Translate cDNAchange to AAchange
      
      Mutant[, AAchange:=annotate_cDNAchange(.SD), by=1:nrow(Mutant)]
      idx <- match(reads$IUPAC_cDNA, Mutant$IUPAC_cDNA)
      reads[,AAchange:=Mutant$AAchange[idx]]
      
      reads <- cbind(reads, parse_IUPAC_AAchange(MutationColumn = reads$AAchange))
      reads <- reads[ REF_AA != ALT_AA,]
      
      reads[, AAchange := gsub("p.","",AAchange,fixed = T)]
      reads[, `:=` (mutFreq = .N), by = qname]
      INFO <- reads[, .(muts = paste(AAchange,collapse = "/"),mut_count = .N), by = qname]
      FREQ <- INFO[,.N,.(muts)]
      setorder(FREQ,-N)
      FREQ[,Sample := gsub("[0-9]+\\.", "", SAMPLE)]
    } else {
      INFO <- data.frame(qname = NA,
                         muts = NA,
                         mut_count = NA)
      FREQ <-
        data.frame(
          muts = NA,
          N = NA,
          Sample = gsub("[0-9]+\\.", "", SAMPLE)
        )
    }
    return(FREQ)
    
  }
  dirname <- format(Sys.time(),"%Y%m%d_%H%M%S")
  dir.create(path = paste0("./",dirname))
  setwd(paste0("./",dirname))
  Mutation <- chr_mutation_vector
  write(x = Mutation,file = "searching.txt",ncolumns = 1)
  for(SAMPLE in searchSamples){
    # SAMPLE <- searchSamples[4]
    # Mutation <- Mutation[1]
    FILE <-
      paste0(
        "/Volumes/UTU_KE_2TB/Seafile/NGS-Data-DC/PacBio/20210909 iSCREAM Library and samples/analysis/sam2tsv/results/",
        SAMPLE,
        "_sam.tsv.gz"
      )
    Mut <- fread(file = FILE)
    colnames(Mut) <- c("qname","ref","pos","alt")
    setkey(Mut, qname)
    
    # Parallel processing.
    mergedData <- parallel::mclapply(
      X = as.list(Mutation),
      FUN = function(X)
        plot_reads(Mutation = X),
      mc.cores = parallel::detectCores()
    )
    # print(head(mergedData))
    mergedData <- unique(rbindlist(mergedData))
    setorder(mergedData,-N)
    
    write.table(
      x = mergedData,
      file = paste0(gsub("[0-9]+\\.", "", SAMPLE),
                    ".csv"),
      sep = ",",
      quote = F,
      row.names = F,
      col.names = T
    )
    print(paste0("Processed ",SAMPLE))
  }
  setwd("..")
}


selectTrueHits <- function(index){
  ## Description: finds the cDNA change for the hits (given as AAchange)
  ## Input: index of the hits you want to get cDNA changes for
  ## Output: returns index where the AAchange in the hit has the correct corresponding cDNA change 
  ## Utility: Due to codon degenracy there can be multiple muts for a basic AA change e.g. A123L
  ##          This function finds the correct nt change correspinding to a given AA change.
  
  idx <- which(hits$AAchange == iSCREAM$AA_change[index])
  
  if(iSCREAM$AA_change[index] == hits$AAchange[idx] & iSCREAM$Ref_nt[index] == hits$Ref_nt[idx] & iSCREAM$Alt_nt[index] == hits$Alt_nt[idx] ){
    return(index)
  } else {
    return(NA)
  }
}


find_mutID_to_test <- function(hits=NULL){
  ## Description: Return matching mutation IDs 
  ##              by matching AAchange to exact nucleotide change
  toCheck <- which(iSCREAM$AA_change %in% hits$AAchange)
  trueMatches <- na.exclude(sapply(toCheck, FUN = function(X) selectTrueHits(X)))
  return(iSCREAM$MutIDs[trueMatches])
}

# -----------> Setting up iSCREAM data DF from PacBio <-----------
# reading file
iSCREAM <- readRDS("/Volumes/UTU_KE_2TB/Seafile/NGS-Data-DC/PacBio/20210909 iSCREAM Library and samples/analysis/filtered/R.Binary/ERBB3_V1_iSCREAM.RDS")

# calculating NRG -> No lig Fold change and sorting the DT based on that
iSCREAM[, FC := NRG_to_NoLig_V1_VF / NRG_V1_VF]
setorder(iSCREAM, -FC, na.last = T)
# Finding Ref and Alt nucleutides from the IUPAC nt change
ref_start <- stringi::stri_locate_first_fixed(str = iSCREAM$MutIDs,pattern = "c.")[,"end"]
ref_end <- stringi::stri_locate_first_regex(str = iSCREAM$MutIDs,pattern = "[0-9]+")[,"start"]
iSCREAM[,Ref_nt := substring(iSCREAM$MutIDs,ref_start+1,ref_end-1)]

alt_start <- stringi::stri_locate_first_fixed(str = iSCREAM$MutIDs,pattern = ">")[,"end"]
iSCREAM[,Alt_nt := substring(iSCREAM $MutID,alt_start+1)]
rm(ref_start,alt_start,ref_end)

# Removing "p." from the text in AA_change column
iSCREAM[,AA_change := gsub("p.","",AA_change, fixed = T)]
