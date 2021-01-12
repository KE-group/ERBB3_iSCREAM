rm(list=ls())
setwd("/Users/deepankar/Documents/Seafile/NGS-Data-DC/BaseSpace/20201125 Marika ERBB3 iSCREAM/Analysis")

dir.create(file.path(".", "variants"), showWarnings = FALSE)
dir.create(file.path(".", "variants/TXT"), showWarnings = FALSE)
dir.create(file.path(".", "variants/R.Binary/"), showWarnings = FALSE)

rm(list=ls())
setwd("/Users/deepankar/Documents/Seafile/NGS-Data-DC/BaseSpace/20201125 Marika ERBB3 iSCREAM/Analysis/parsedVCF/R.Binary/")

fnames=list.files(pattern = ".RDS")

for(file in fnames){
  # file= fnames[1]
  print (paste("Processing",file))
  rm(list=ls()[!ls() %in% c("file","fnames")])
  tab=readRDS(file = file)
  # str(tab)
  tab=tab[,-which(colnames(tab)%in%c("rsID","Filter","Quality"))] # Removing zero and dot (.) columns
  #Change Depth to numeric
  tab$DP=as.numeric(tab$DP)
  subtab <- tab
  subtab$Alt=gsub("<*>","XX",subtab$Alt,fixed = T)
  
  variant_count <- length(unlist(strsplit(subtab$Alt,",",fixed = T),use.names = F))
  
  # require(progress);pb <- progress_bar$new( format = "  Processing [:bar] :current SNVs processed in - :elapsed", total = variant_count, clear = FALSE, width= 65)
  # pb$tick(0)
  chr <- pos <- ref <- alt <- DP <- AD <- rep(NA,variant_count); count <- 0
  for(i in seq(1,length(subtab$DP))){
  #for(i in seq(1,100)){
    datarow=subtab[i,]
    Alts=unlist(strsplit(datarow[,"Alt"],",",fixed = T),use.names = F)
    ADs=unlist(strsplit(datarow[,"AD"],",",fixed = T),use.names = F)
    for(j in seq(1,length(Alts))){
      count=count+1
      #chr=unlist(c(chr,datarow[1]))
      chr[count] <- unlist(datarow["Chromosome"],use.names = F)
      pos[count] <- unlist(datarow["Position"],use.names = F)
      ref[count] <- unlist(datarow["Ref_Base"],use.names = F)
      alt[count] <- Alts[j]
      DP[count] <- unlist(datarow["DP"],use.names = F)
      AD[count] <- ADs[j+1]
      # if(count%%10==0){
      #   pb$tick(10)
      # }
    }
  }
  AD=as.numeric(AD)
  SNV.matrix=data.frame(chr,pos,ref,alt,DP,AD,stringsAsFactors = F)
  rm(chr,pos,ref,alt,DP,AD,Alts,ADs,datarow,i,j)
  
  # Keeping Mutations with at least 10 000 reads 
  # This is an exceptionally high Depth threshold as 
  # the sequencing depth is really astronomical!! mean(SNV.matrix$DP) is 2 million +
  SNV.matrix <- SNV.matrix[SNV.matrix$DP>=10000,]
  
  # Removing all <*> variants
  SNV.matrix <- SNV.matrix[SNV.matrix$alt!="XX",]
  
  fnm=gsub(".vcf.RDS","",file)
  saveRDS(object = SNV.matrix,file = paste("../../variants/R.Binary/",fnm,".RDS",sep="")) # use var=readRDS("FileName")
  write.table(SNV.matrix,file = paste("../../variants/TXT/",fnm,".tsv",sep=""),sep="\t",row.names = F,col.names = T,quote = F,na = "NA")
}
