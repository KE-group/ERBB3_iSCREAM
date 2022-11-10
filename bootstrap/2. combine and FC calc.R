


#--------> Combining files <------------
rm(list=ls())
setwd("/Volumes/DATA/Seafile/NGS-Data-DC/PacBio/20210909 iSCREAM Run1/analysis/Figures/Bootstrap/Data/")
dir.create("./combined_hits", showWarnings = F)

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

# Pooling NRG_to_NoLig data
for(FILE in data.files){
  print(FILE)
}

tmp <- fread(data.files[1])
idx <- match(x = dat$genotypes,table = tmp$muts)
dat$NRG_to_NoLig <- tmp$N[idx]

tmp <- fread(data.files[2])
idx <- match(x = dat$genotypes,table = tmp$muts)
dat$NRG <- tmp$N[idx]

tmp <- fread(data.files[3])
idx <- match(x = dat$genotypes,table = tmp$muts)
dat$WEHI <- tmp$N[idx]
rm(tmp)
setorder(dat, -NRG_to_NoLig,na.last = T)

write.table(
  x = dat,
  file = "combined_hits/read_counts_ERBB3_V1_NRG_to_NoLig.tsv",
  quote = F,
  row.names = F,
  col.names = T,
  sep = "\t"
)