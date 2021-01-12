# To parse VCF files generated from samtools mpileup & bcftools
rm(list=ls())
setwd("/Users/deepankar/Documents/Seafile/NGS-Data-DC/BaseSpace/20201125 Marika ERBB3 iSCREAM/Analysis")

dir.create(file.path(".", "parsedVCF"), showWarnings = FALSE)
dir.create(file.path(".", "parsedVCF/TXT"), showWarnings = FALSE)
dir.create(file.path(".", "parsedVCF/R.Binary/"), showWarnings = FALSE)

setwd("/Users/deepankar/Documents/Seafile/NGS-Data-DC/BaseSpace/20201125 Marika ERBB3 iSCREAM/Analysis/vcf/")

fnames=list.files(pattern = ".vcf")

headers=readLines("FILE.HEADERS.txt")
FORMAT.head.idx=grep("^##FORMAT=<ID=",headers)
FORMAT.head=headers[FORMAT.head.idx];rm(FORMAT.head.idx)
searchI=gregexpr(",",FORMAT.head)
start.posI=unlist(lapply(searchI, `[[`, 1))
temp=substr(FORMAT.head,0,start.posI-1)
FORMAT.head=gsub("##FORMAT=<ID=","",temp);rm(temp,searchI,start.posI)

INFO.head.idx=grep("^##INFO=<ID=",headers)
INFO.head=headers[INFO.head.idx];rm(INFO.head.idx)
searchI=gregexpr(",",INFO.head)
start.posI=unlist(lapply(searchI, `[[`, 1))
temp=substr(INFO.head,0,start.posI-1)
INFO.head=gsub("##INFO=<ID=","",temp);rm(temp,searchI,start.posI)

rm(headers)
# require(progress)

for (f in fnames){
  # f=fnames[1]
  print (paste("Processing",f))
  tempvcf=readr::read_tsv(file = f,comment = "##")
  class(tempvcf)="data.frame"
  tempvcf$SampleName=rep(gsub("../BAM/","",colnames(tempvcf)[10],fixed = T),length(tempvcf[,1]))
  colnames(tempvcf)=c("Chromosome","Position","rsID","Ref_Base","Alt_Base","Quality","Filter","INFO","FORMAT","DATA","SampleID")
  parsed=matrix("",dim(tempvcf)[1],length(FORMAT.head))
  colnames(parsed)=FORMAT.head
  message("Parsing DATA column")
  # pb <- progress_bar$new(format = "  Processing [:bar] :current sites processed in - :elapsed", total = nrow(tempvcf), clear = FALSE, width= 65)
  # pb$tick(0)
  for(i in seq(1,dim(tempvcf)[1])){
    #i=1
    # print(i)
    nm=unlist(strsplit(x = tempvcf[i,"FORMAT"],split = ":",fixed = T),use.names = F)
    d=unlist(strsplit(x =tempvcf[i,"DATA"],split = ":",fixed = T),use.names = F)
    parsed[i,nm]=d
    # if(i%%10==0){
    #   pb$tick(10)
    # }
  }
  parsed=data.frame(parsed,stringsAsFactors = F)
  tempvcf=tempvcf[,-which(colnames(tempvcf)%in%c("FORMAT","DATA"))]
  tempvcf=cbind(tempvcf,parsed)
  #----------------------------------
  parsed=matrix("",dim(tempvcf)[1],length(INFO.head))
  colnames(parsed)=INFO.head
  message("Parsing INFO column")
  # pb <- progress_bar$new(format = "  Processing [:bar] :current sites processed in - :elapsed", total = nrow(tempvcf), clear = FALSE, width= 65)
  # pb$tick(0)
  # for(i in seq(1,10000)){
  for(i in seq(1,dim(tempvcf)[1])){
    # print(i)
    # i=4905
    dat=unlist(strsplit(tempvcf[i,"INFO"],";",fixed = T),use.names = F)
    searchI=gregexpr("=",dat)
    pos=unlist(lapply(searchI, `[[`, 1));rm(searchI)
    keepsame=which(pos==(-1))
    if(length(keepsame)!=0){
      # print(i)
      nm=substr(x = dat,start = 0,stop = pos-1)
      d=substr(dat,pos+1,nchar(dat))
      nm[keepsame]=d[keepsame]=dat[keepsame]
    }
    else{
      nm=substr(x = dat,start = 0,stop = pos-1)
      d=substr(dat,pos+1,nchar(dat))
    }
    parsed[i,nm]=d  
    # if(i%%1000==0){
    #   pb$tick(1000)
    # }
  }
  parsed=data.frame(parsed,stringsAsFactors = F)
  tempvcf=tempvcf[,-which(colnames(tempvcf)%in%c("INFO"))]
  tempvcf=cbind(tempvcf,parsed)
  tempvcf <- tempvcf[,-which(colnames(tempvcf) %in% c("INDEL","IDV","IMF"))]
  rm(list=ls()[!ls() %in% c("tempvcf","f","fnames","FORMAT.head","INFO.head")])
  SampleNm=gsub(".vcf","",f)
  assign(SampleNm,tempvcf)
  rm(tempvcf)
  write.table(get(SampleNm),file = paste("../parsedVCF/TXT/",f,".tsv",sep=""),sep="\t",col.names = T,row.names = F,quote = F,na = "NA")
  saveRDS(object =  get(SampleNm),file = paste("../parsedVCF/R.Binary/",f,".RDS",sep=""))
}
