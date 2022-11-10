#-------------
# cDNA change
rm(list=ls())

setwd("/Volumes/DATA/Seafile/NGS-Data-DC/BaseSpace/20201125 Marika ERBB3 iSCREAM/Analysis/Figures/")

dat=readRDS("20220506.plotMuts.RDS")

bases=c("A","C","G","T")
mut.mat=matrix(data=NA, nrow = 4, ncol = 4)
colnames(mut.mat)=row.names(mut.mat)=bases
searchI=gregexpr(">",dat$MutID,fixed = T)
start.posI=unlist(lapply(searchI, `[[`, 1),use.names = F)
dat$ref=substr(dat$MutID,start.posI-1,start.posI-1)
dat$alt=substr(dat$MutID,start.posI+1,start.posI+1)
rm(searchI,start.posI)

subDF=dat[,c("AAchange","ref","alt")]


for(i in bases){
  temp.i=subset(subDF,subDF$ref==i)
  for(j in unique(temp.i$alt)){
    temp.j=subset(temp.i,temp.i$alt==j)
    nobs=dim(temp.j)[1]
    mut.mat[i,j]=nobs
    print(paste(i,j,nobs))
  }
}
rm(i,j,temp.i,temp.j)
breakList=seq(1, max(mut.mat,na.rm = T), by = 5)

# colours=colorRampPalette(c("navyblue","maroon4","red1"))(length(breakList))
colours=viridis::plasma(length(breakList), direction = -1, end = 0.9)

library(pheatmap)
pheatmap(t(mut.mat),color=colours,cluster_rows = F,cluster_cols = F,border_color = "black",cellwidth = 25,cellheight = 25,display_numbers = T,number_color = "white",number_format="%.0f",main="EGFR iSCREAM DNA change Distribution",filename = "DNAchange distribution.pdf")

mut.mat[is.na(mut.mat)]=0

longDF=reshape2::melt(mut.mat)
longDF$Var1=as.character(longDF$Var1);longDF$Var2=as.character(longDF$Var2)
longDF$Mutation=paste(longDF$Var1,">",longDF$Var2,sep="")
longDF=longDF[longDF$value!=0,]


transitions=c("A>G","G>A","C>T","T>C")
transversions=c("A>C","C>A","A>T","T>A","C>G","G>C","G>T","T>G")

longDF$Type=longDF$Mutation
longDF$Type[longDF$Type%in%transitions]="Transitions"
longDF$Type[longDF$Type%in%transversions]="Transversions"

a1=longDF$value[longDF$Type=="Transitions"]
a2=longDF$value[longDF$Type=="Transversions"]
longDF <- as.data.table(longDF)
t.test(a1,a2)
wilcox.test(a1,a2)

#-----------

#-------------
# AA change
rm(list=ls())

setwd("/Volumes/DATA/Seafile/NGS-Data-DC/BaseSpace/20201125 Marika ERBB3 iSCREAM/Analysis/Figures/")

dat=readRDS("20220506.plotMuts.RDS")
subDF=dat[-grep("?",AAchange,fixed=T),c("AAchange","AAPos")] # removing start codon muts

Muts=subDF$AAchange;temp=nchar(Muts)
subDF$refAA=substr(x = Muts,start = 0,stop = 1)
subDF$altAA=substr(x = Muts,start = temp,stop = temp)
rm(Muts,temp)

AAs=sort(unique(c(subDF$refAA,subDF$altAA)))

mut.mat=matrix(data=NA, nrow = length(AAs), ncol = length(AAs))
colnames(mut.mat)=row.names(mut.mat)=AAs


notPointMuts=readRDS("~/OneDrive - O365 Turun yliopisto/Git/GitLab.UTU/EleniusGroup/ERBB.Miscellaneous/Posssible mutants/ERBB3/NotPointMuts.ERBB3.RDS")


AAs=sort(unique(c(notPointMuts$refAA,notPointMuts$altAA)))
for(i in AAs){
  temp.i=subset(notPointMuts,notPointMuts$refAA==i)
  if(length(temp.i$ref)>0){
    for(j in unique(temp.i$altAA)){
      temp.j=subset(temp.i,temp.i$altAA==j)
      mut.mat[i,j]=0
      # print(paste(i,j,nobs))
    } 
  }
}
rm(AAs,temp.i,temp.j)

AAs=sort(unique(c(subDF$refAA,subDF$altAA)))

for(i in AAs){
  temp.i=subset(subDF,subDF$refAA==i)
  for(j in unique(temp.i$altAA)){
    temp.j=subset(temp.i,temp.i$altAA==j)
    nobs=dim(temp.j)[1]
    mut.mat[i,j]=nobs
    # print(paste(i,j,nobs))
  }
}
rm(AAs,temp.i,temp.j,nobs)


mut.mat=t(mut.mat) # Rows as refaa and cols as altaa

breakList=seq(1, max(mut.mat,na.rm = T), by = 1)
breakList=c(-1,breakList)

# colours=colorRampPalette(c("navyblue","maroon4","red1"))(length(breakList)-1)
colours=viridis::inferno(length(breakList)-1, direction = -1, end = 0.9)
colours=c("grey",colours)

library(pheatmap)

# pheatmap(mut.mat,color=colours,cluster_rows = F,cluster_cols = F,border_color = "black",cellwidth = 18,cellheight = 18,display_numbers = T,number_color = "white",number_format="%.0f",main="EGFR iSCREAM AA mutation distribution")

pheatmap(
  mut.mat,
  color = colours,
  cluster_rows = F,
  cluster_cols = F,
  border_color = "black",
  cellwidth = 18,
  cellheight = 18,
  display_numbers = T,
  number_color = "white",
  number_format = "%.0f",
  main = "EGFR iSCREAM AA mutation distribution",
  filename = "AAchange distribution.pdf"
)
