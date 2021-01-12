setwd("/Users/deepankar/Documents/Seafile/NGS-Data-DC/BaseSpace/20201125 Marika ERBB3 iSCREAM/Analysis/variants/R.Binary/")

rm(list=ls())

files <- list.files(".",pattern = ".RDS")
files <- c(files[grep(pattern = "Plasmid",
                      x = files)], 
           files[grep(pattern = "V2",
                      x = files)])
# files=gtools::mixedsort(files) # Optional step. This gets the column names in the file as "Sample#1, Sample#2,..." instead of #1,#10,#11,...
allMuts=c()
for(f in files){
  # f <- files[1]
  tab=readRDS(f)
  MutID=paste(tab$chr,":",tab$pos,tab$ref,">",tab$alt,sep="")
  allMuts=c(allMuts,MutID)
}
rm(tab)
allMuts=unique(allMuts)

searchI=gregexpr(":",allMuts,fixed = T)
start.posI=unlist(lapply(searchI, `[[`, 1),use.names = F)
chr=substr(allMuts,0,start.posI-1)

searchI=gregexpr(">",allMuts,fixed = T)
start.posI=unlist(lapply(searchI, `[[`, 1),use.names = F)
mut=substr(allMuts,start.posI-1,nchar(allMuts))

searchI=gregexpr(">",mut,fixed = T)
start.posI=unlist(lapply(searchI, `[[`, 1),use.names = F)
ref=substr(mut,start.posI-1,start.posI-1)
alt=substr(mut,start.posI+1,start.posI+1);rm(mut,searchI,start.posI)

searchI=gregexpr(":",allMuts,fixed = T)
start.pos=unlist(lapply(searchI, `[[`, 1),use.names = F)
searchI=gregexpr(">",allMuts,fixed = T)
end.pos=unlist(lapply(searchI, `[[`, 1),use.names = F)
pos=substr(allMuts,start.pos+1,end.pos-2);rm(searchI,start.pos,end.pos)

Mutation.Table=data.frame(allMuts,chr,pos,ref,alt,stringsAsFactors = F)

rm(list=ls()[!ls() %in% c("Mutation.Table")])

# ERBB3 location chr12:56,473,641-56,497,291
Mutation.Table=subset(Mutation.Table,Mutation.Table$chr=="chr12" & Mutation.Table$pos>=56473641 & Mutation.Table$pos <= 56497291)

AnnoDF=Mutation.Table[,c(2,3,3,4,5,1)]

setwd("/Users/deepankar/Documents/Seafile/NGS-Data-DC/BaseSpace/20201125 Marika ERBB3 iSCREAM/Analysis/Annovar/")
# write.table(AnnoDF,file = "20201215.InputAnnovar.txt",col.names = F,quote = F,row.names = F)

#- run following lines to setup for annovar
# cd '/Users/deepankar/Documents/Seafile/NGS_Seq_Tools/annovar2'
# ln -s "/Users/deepankar/Documents/Seafile/NGS-Data-DC/BaseSpace/20201125 Marika ERBB3 iSCREAM/Analysis/Annovar/20201215.InputAnnovar.txt" "20201215.InputAnnovar.txt"

# Run annovar
# script 20201215.ERBB3.iSCREAM.log command ./table_annovar.pl 20201215.InputAnnovar.txt humandb -buildver hg19 -out ERBB3.iSCREAM -remove -protocol refgene,gnomad211_exome,gnomad211_genome,avsnp150,clinvar_20190305 -operation g,f,f,f,f -nastring . -csvout -polish

# minimal annotations
# annotate_variation.pl -out ERBB2.txt --build hg19 20190103.InputAnnovar.txt humandb

rm(AnnoDF,Mutation.Table)
var=read.table("ERBB3.iSCREAM.hg19_multianno.csv",sep=",",stringsAsFactors = F,header = T)
var$MutID=paste(var$Chr,":",var$Start,var$Ref,">",var$Alt,sep="")
saveRDS(object = var,file = "ERBB3.iSCREAM.hg19_multianno.RDS")
#-----------

