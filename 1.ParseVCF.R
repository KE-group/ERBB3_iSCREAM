# To parse VCF files generated from samtools mpileup & bcftools
rm(list=ls())
setwd("~/BaseSpace/20180112 EGFR Trimmed HiDP/vcf")
# file="20180112 EGFR_lib.EGFR.locus.HiDP.vcf"
# file="20180112 EGFR.No.Lig.EGFR.locus.HiDP.vcf"
# file="20180808.EGFRwEGF_P1.EGFR.locus.HiDP.vcf"
# file="20180808.EGFRwEGF_P5.EGFR.locus.HiDP.vcf"
file="20180808.EGFRwIL3P1.EGFR.locus.HiDP.vcf"

temptext=readLines(paste(file,sep=""))
temptext=temptext[-grep("^##",x = temptext)]
require(progress)
pb <- progress_bar$new( format = "  Processing [:bar] :percent in :elapsed", total = length(temptext), clear = FALSE, width= 60)
temptext=gsub(";","\t",temptext)
temptext=temptext[-seq(1:49)]
parsing=matrix(NA,length(temptext),19)
colnames(parsing)=c("chr","pos","dot1","Ref","Alt","zero","dot2","DP","I16","QS","VDB","SGB","RPB","MQB","MQSB","BQB","MQ0F","Format","DATA")
for( i in seq(1,length(temptext))){
#for( i in seq(69801,69810)){
  line=temptext[i]
  words=unlist(strsplit(line,"\t"))
  parsing[i,1:10]=words[1:10]
  for (j in seq(11,19)){
    if(is.na(words[j])|words[j]==""){
      break
    }
    #print(regexpr("VDB",words[j]))
    if(regexpr("VDB",words[j])==1){
      parsing[i,"VDB"]=words[j]
    }
    else if(regexpr("SGB",words[j])==1){
      parsing[i,"SGB"]=words[j]
    }
    else if(regexpr("RPB",words[j])==1){
      parsing[i,"RPB"]=words[j]
    }
    else if(regexpr("MQB",words[j])==1){
      parsing[i,"MQB"]=words[j]
    }
    else if(regexpr("MQSB",words[j])==1){
      parsing[i,"MQSB"]=words[j]
    }
    else if(regexpr("BQB",words[j])==1){
      parsing[i,"BQB"]=words[j]
    }
    else if(regexpr("MQ0F",words[j])==1){
      parsing[i,"MQ0F"]=words[j]
    }
    else if(regexpr("PL:AD",words[j])==1){
      parsing[i,"Format"]=words[j]
      parsing[i,"DATA"]=words[j+1]
    }
  }
  if(i%%1000==0){
    pb$tick(1000)
  }
}
temp=data.frame(parsing,stringsAsFactors = F)

pb <- progress_bar$new( format = "  Processing [:bar] :percent in :elapsed", total = length(temp$DATA), clear = FALSE, width= 60)
PL=c();AD=c()
for(i in seq(1,length(temp$DATA))){
#for(i in seq(1,10)){
  #print(temp$DATA[i])
  var=unlist(strsplit(temp$DATA[i],":"))
  PL=c(PL,var[1])
  AD=c(AD,var[2])
  if(i%%1000==0){
    pb$tick(1000)
  }
}
temp=cbind(temp,PL,AD)
temp=temp[,c(-18,-19)]
colnames(temp)[c(19)]="Allelic.Depth"
parsing=temp;rm(temp)

assign(file,parsing);rm(parsing)
file.s=paste("Parsed/P-",file,sep="")
saveRDS(object =get(file),file = paste(file.s,".RDS",sep=""))
write.table(get(file),file=paste(file.s,".tsv",sep=""),sep="\t",row.names = F,col.names = T,quote = F,na = "NA")


