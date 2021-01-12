rm(list=ls())
setwd("/Users/deepankar/Documents/Seafile/NGS-Data-DC/BaseSpace/20201125 Marika ERBB3 iSCREAM/Analysis/variants/R.Binary/")
allMuts <- readRDS("../../Annovar/ERBB3.iSCREAM.hg19_multianno.RDS")
toKeep <- c("Chr","Start","End","Ref","Alt","Func.refGene","Gene.refGene","ExonicFunc.refGene","AAChange.refGene","MutID")
Mutation.Table <- allMuts[,toKeep]
rm(allMuts,toKeep)

files <- list.files(".",pattern = ".RDS")
files <- c(files[grep(pattern = "Plasmid",
                      x = files)], 
           files[grep(pattern = "V2",
                      x = files)])
for(f in files){
  # f <- files[1]
  print(paste("Processing",f))
  tab=readRDS(f)
  tab$VF=(tab$AD*100)/tab$DP
  # tab$Rel.Ab=(tab$AD*100)/sum(tab$AD)
  tab$MutID=paste(tab$chr,":",tab$pos,tab$ref,">",tab$alt,sep="")
  idx=match(Mutation.Table$MutID,tab$MutID)
  colnames(tab)=paste(gsub("RDS","",f),colnames(tab),sep="")
  Mutation.Table=cbind.data.frame(Mutation.Table,tab[idx,c(5:7)])
}
rm(tab,f,files,idx)

source("https://github.com/dchakro/shared_Rscripts/raw/master/Annovar/AnnovarMutCodeFind.R")
Mutation.Table$AAchange=annovarMutCodeFind(Mutation.Table$AAChange.refGene,isoform = "NM_001982") 
rm(annovarMutCodeFind)

# Filtering

Mutation.Table=Mutation.Table[Mutation.Table$Func.refGene=="exonic",]

Mutation.Table=Mutation.Table[Mutation.Table$ExonicFunc.refGene!="synonymous SNV",]


# writing file
saveRDS(object = Mutation.Table,file = "../../Figures/20201216.ERBB3_iSCREAM_Mutations.RDS")
