#---->  Lollipop <------
setwd("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Marika ERBB3 screen/Figures/mutation landscape")
#---------> cBioportal <------
rm(list=ls())
cBio <- data.table::fread("ERBB3_cbioPortal.tsv")
colnames(cBio) <- gsub(" ",".",colnames(cBio))
cBio <- cBio[Mutation.Type!="Fusion",]
cBio <- cBio[!cBio$Mutation.Type %in% c("Splice_Region","Splice_Site"),]
cBio$Sample.ID <- toupper(cBio$Sample.ID)
cBio$MutID <- paste(cBio$Sample.ID,"=",cBio$Protein.Change,sep="")

unifiedDF <- data.frame(MutID = cBio$MutID, source = rep("cBioportal"))

#---------> Analysis begins <------
table(unifiedDF$source)

tmp <- data.frame(stringi::stri_split_fixed(str = unifiedDF$MutID,pattern = "=",simplify = T))
unifiedDF$sampleID <- tmp$X1
unifiedDF$AAchange <- tmp$X2
rm(tmp)

# unifiedDF <- unifiedDF[-which(duplicated(unifiedDF$MutID)),]
# These six muts were removed as they were duplicate entries
# [1] "DS-UTTCC-033-P=M91I"     "DS-UTTCC-060-P=V119I"   
# [3] "DS-BLA-032=Y148H"        "DS-UTTCC-060-P=R1136C"  
# [5] "P-0007096-T01-IM5=P306R" "P-0004688-T01-IM5=R669H"

table(unifiedDF$source)

# Per Mutation Figure ####
rm(list=ls()[!ls() %in% c("unifiedDF")])
uniqueMuts <- table(unifiedDF$AAchange)
df <- data.frame(AAchange=names(uniqueMuts),count=as.numeric(uniqueMuts),stringsAsFactors = F)
source('https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/MutSiteFind.R')
df$AAPos <- as.numeric(MutSiteFind(df$AAchange))
df <- df[order(df$count,decreasing = T),]
rm(MutSiteFind)

source("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/ACCEPTED/2022 ERBB4 BaF3 screen/Supplementary/Code/MutationClassifier.R")
df$Type <- MutationClassifier(df$AAchange)[,2]
rm(MutationClassifier)
DF2 <- df
rm(df)

rang <-  c(
  missense = "#198c19",
  ambiguous = "#a3a3a3",
  nonsense = "#00D0E5",
  complex = "black",
  InDel = "#800080",
  deletion = "#800080",
  oncogenic = "#e10000",
  likelyOncogenic = "#ff7373"
)

likelyOncogenic <- readLines("LikelyOncogenic.txt")
oncogenic <- readLines("oncogenic.txt")

DF2$Type[DF2$AAchange %in% likelyOncogenic] <- "likelyOncogenic"
DF2$Type[DF2$AAchange %in% oncogenic] <- "oncogenic"

library(ggplot2)
source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type='L',legend = T)
# Generate domain markings
tmp <- read.table(file = "~/OneDrive - O365 Turun yliopisto/Git/GitHub/KE-Group/MultipleMutations/ERBBs Multiple Mutations/ProteinDomains.txt",sep = "\t",header = T)
proteins <-  c("ERBB3")
gene <- "ERBB3"
source("https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/roundUp.R")
value_range <- c(0,roundUp(max(DF2$count,na.rm = T)))
domains <- tmp[tmp$Protein %in% proteins,]
tmp <- domains[domains$Symbol=="FULL",c("Protein","End")]
proteinLength <- as.vector(tmp$End)
names(proteinLength) <- tmp$Protein
domain_markings <- domains[domains$Protein==gene,]

# Generate axis

source("https://github.com/dchakro/shared_Rscripts/raw/master/ggplotBreaks.R")
x_brk <- ggplotBreaks(range = c(0, proteinLength[gene]), tick = 100, skip.steps = 0)
x_brk$breaks <- c(x_brk$breaks,unname(proteinLength[gene]))
x_brk$labels <- c(x_brk$labels,as.character(unname(proteinLength[gene])))
rm(tmp)

y_brk <- ggplotBreaks(range = value_range,tick = 25,skip.steps = 0)
source("https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/roundUp.R")

ggplot(DF2,aes(x = AAPos,y=count))+
  geom_segment(aes(yend=-1,xend=AAPos))+
  geom_point(size = 2,
             shape = 21,
             alpha = 1,
             aes(fill=Type), 
             color = "#000000")+
  scale_x_continuous(
    breaks = x_brk$breaks,
    labels = x_brk$labels,
    limits = c(0, proteinLength[gene])
  )+
  scale_y_continuous(
    breaks = y_brk$breaks,
    labels = y_brk$labels,
    limits = c(-1, roundUp(x = max(DF2$count), to = 10))
  )+
  xlab(paste(gene,"Amino Acid Position"))+
  ylab("No. of mutations")+
  scale_fill_manual(values=rang)+
  ggtitle(paste(gene,"mutations cBioPortal, N=",length(unifiedDF$MutID)))+
  customtheme + 
  geom_rect(
    data = domain_markings[domain_markings$Symbol != "FULL", ],
    mapping = aes(
      xmin = Start,
      xmax = End,
      ymin = 0,
      ymax = (max(value_range) * 0.1)
    ),
    color = "black",
    fill = "black",
    inherit.aes = FALSE
  ) + 
  geom_text(
    data = domain_markings[domain_markings$Symbol != "FULL", ],
    mapping = aes(
      x = Start + (End - Start) / 2,
      y = 0 + ((max(value_range) * 0.1) - 0) / 2,
      label = Symbol
    ),
    fontface = "bold",
    color = "#FFFFFF",
    size = 4,
    inherit.aes = FALSE
  ) +
  geom_text(data = DF2[order(DF2$count, decreasing = T)[1:12], ], aes(x = AAPos, y =
                                                                       count, label = AAchange))
# ggsave(filename = "Lollipop_AA.pdf",height = 3,width = 14)

DF2 <- DF2[,c("AAPos","AAchange","count","reducedType","Type","breakdown")]
DF2 <- DF2[order(DF2$count,decreasing = T),]



# Per residue Figure ####

uniqueMuts <- table(unifiedDF$AAchange)
df <- data.frame(AAchange=names(uniqueMuts),count=as.numeric(uniqueMuts),stringsAsFactors = F)
source('https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/MutSiteFind.R')
df$AAPos <- as.numeric(MutSiteFind(df$AAchange))
df <- df[order(df$count,decreasing = T),]
rm(MutSiteFind)

source("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/ACCEPTED/2022 ERBB4 BaF3 screen/Supplementary/Code/MutationClassifier.R")
df$Type <- MutationClassifier(df$AAchange)[,2]
rm(MutationClassifier)

rang=c(missense="#198c19",
       ambiguous="#ff7373",
       nonsense="#800080",
       complex="black",
       InDel="brown",
       deletion="brown",
       insertion="brown")

DF2 <- aggregate(AAchange ~ AAPos,data=df,FUN = function(X) paste(unique(X), collapse="/"))
DF2$count <- aggregate(count ~ AAPos,data=df,FUN = function(X) sum(X))[,2]
DF2$breakdown <- aggregate(count ~ AAPos,data=df,FUN = function(X) paste(X, collapse="/"))[,2]
DF2$Type <- aggregate(Type ~ AAPos,data=df,FUN = function(X) paste(unique(X), collapse="/"))[,2]
DF2 <- DF2[,c("AAchange", "AAPos", "count", "breakdown","Type")]
rm(df,uniqueMuts)

for(i in grep(pattern = "/",x = DF2$AAchange,fixed = T)){
  words <- unlist(strsplit(x = DF2$AAchange[i],split = "/",fixed = T),use.names = F)
  # print(words)
  first <- stringi::stri_extract_first_regex(pattern = "[ACDEFGHIKLMNPQRSTVWYX*]?[0-9]+(_[ACDEFGHIKLMNPQRSTVWYX*]?[0-9]+)?", str = words[1])
  words[2:length(words)] <- gsub(first,"",words[2:length(words)])
  DF2$AAchange[i] <- paste(words,collapse = "/")
}

# Assigning color of "most popular mutation"
DF2$reducedType <- unlist(lapply(strsplit(x = DF2$Type,split = "/",fixed = T),`[[`, 1))


library(ggplot2)
source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type='L',legend = F)
# Generate domain markings
tmp <- read.table(file = "~/OneDrive - O365 Turun yliopisto/Git/GitHub/KE-Group/MultipleMutations/ERBBs Multiple Mutations/ProteinDomains.txt",sep = "\t",header = T)
proteins <-  c("ERBB3")
gene <- "ERBB3"
value_range <- c(0,max(DF2$count,na.rm = T))
domains <- tmp[tmp$Protein %in% proteins,]
tmp <- domains[domains$Symbol=="FULL",c("Protein","End")]
proteinLength <- as.vector(tmp$End)
names(proteinLength) <- tmp$Protein
domain_markings <- domains[domains$Protein==gene,]

# Generate axis

source("https://github.com/dchakro/shared_Rscripts/raw/master/ggplotBreaks.R")
x_brk <- ggplotBreaks(range = c(0, proteinLength[gene]), tick = 50, skip.steps = 1)
x_brk$breaks <- c(x_brk$breaks,unname(proteinLength[gene]))
x_brk$labels <- c(x_brk$labels,as.character(unname(proteinLength[gene])))
rm(tmp)

y_brk <- ggplotBreaks(range = value_range,tick = 25,skip.steps = 0)

ggplot(DF2,aes(x = AAPos,y=count))+
  geom_segment(aes(yend=-1,xend=AAPos))+
  geom_point(size = 2,
             shape = 21,
             alpha = 1,
             aes(fill=reducedType), 
             color = "#000000")+
  scale_x_continuous(breaks = x_brk$breaks, 
                     labels = x_brk$labels, 
                     limits = c(0, proteinLength[gene]))+
  scale_y_continuous(breaks = y_brk$breaks, 
                     labels = y_brk$labels,
                     limits = c(-1,max(DF2$count)))+
  xlab(paste(gene,"Amino Acid Position"))+ylab("No. of mutations")+
  scale_fill_manual(values=rang)+
  ggtitle(paste(gene,"mutations cBioPortal, N=",
                length(unifiedDF$MutID)))+
  customtheme + 
  geom_rect(data=domain_markings[domain_markings$Symbol!="FULL",], 
            mapping = aes(xmin=Start, xmax=End, ymin=0, ymax=(max(value_range) *0.1)),
            color="black",
            fill="black",
            inherit.aes = FALSE) + 
  geom_text(data=domain_markings[domain_markings$Symbol!="FULL",], 
            mapping = aes(x=Start+(End-Start)/2, 
                          y=0+((max(value_range) *0.1)-0)/2, 
                          label=Symbol),
            fontface="bold",
            color="#FFFFFF", 
            size=4, 
            inherit.aes = FALSE) +
  geom_text(data=DF2[order(DF2$count,decreasing = T)[1:5],],
            aes(x=AAPos,
                y=count,
                label=AAchange))

# ggsave(filename = "Lollipop_Residue.pdf",height = 3,width = 12)

DF2 <- DF2[,c("AAPos","AAchange","count","reducedType","Type","breakdown")]
DF2 <- DF2[order(DF2$count,decreasing = T),]
