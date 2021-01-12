options(scipen = 100, digits = 4)
rm(list=ls())

require(ggplot2)
source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type='L',ticks = "in",fontsize.cex = 1.2)

setwd("/Users/deepankar/Documents/Seafile/NGS-Data-DC/BaseSpace/20201125 Marika ERBB3 iSCREAM/Analysis/Figures/")
allMuts <- readRDS("20201216.ERBB3_iSCREAM_Mutations.RDS")


#--------->
# Compare concordance in VFs
value <- lm(allMuts$Plasmid.VF~allMuts$Plasmid_PCR.VF)
anno_text <- paste("Slope:",signif(value$coefficients[[2]],2),"\nIntercept:",signif(value$coefficients[[1]],1));rm(value)

p <- ggplot(data = allMuts,aes(x=Plasmid.VF, y=Plasmid_PCR.VF)) + geom_point(alpha=0.15) + geom_smooth(method = "lm", se = F, color="#f6546a") + xlab("Plasmid Library") + ylab("PCR from the Plasmid Library") + ggtitle("VF (1 dot = 1 mutation) compared between samples") +  annotate("text",x=0.05,y=0.35,label=anno_text) + customtheme
ggsave(filename = "Cor_Lib vs PCR.svg",plot = p,width = 5,height = 5);rm(p,anno_text)

value <- lm(allMuts$Plasmid.VF~allMuts$V2_WEHI.VF)
anno_text <- paste("Slope:",signif(value$coefficients[[2]],2),"\nIntercept:",signif(value$coefficients[[1]],1));rm(value)
p <- ggplot(data = allMuts,aes(x=Plasmid.VF, y=V2_WEHI.VF)) + geom_point(alpha=0.15) + geom_smooth(method = "lm", se = F, color="#f6546a") + xlab("Plasmid Library") + ylab("V2/+ WEHI sample") + ggtitle("VF (1 dot = 1 mutation) compared between samples")+ annotate("text",x=0.3,y=1.5,label=anno_text)  + customtheme

ggsave(filename = "Cor_Lib vs WEHI.svg",plot = p,width = 5,height = 5);rm(p,anno_text)
#--------->
# Depth
subdf=allMuts[,c("MutID","Plasmid_PCR.DP","Plasmid.DP","V2_WEHI.DP","V2_NRG2noIL3.DP","V2_NRG.DP")]
subdf.long=reshape2::melt(data = subdf)
y.breaks=seq(0,ceiling(max(subdf.long$value,na.rm = T)),by = 500000)

p <- ggplot(data = subdf.long,aes(x=variable,y=value))+geom_jitter(height = 0,alpha=0.3,size=1,color="grey70")+geom_violin(color="black",fill=NA,trim = F,adjust=2,na.rm = T,scale = "width",width=0.8)+geom_boxplot(outlier.color = NA,color="black",fill=NA,width=0.1)+scale_y_continuous(breaks = y.breaks,limits = c(0,max(subdf.long$value,na.rm = T)),labels = scales::comma)+ylab("Total Depth for each variant")+xlab("Sample")+customtheme
ggsave(filename = "ReadDepth_distrib.svg",plot = p,width = 10,height = 5);rm(y.breaks,subdf,subdf.long)


#----------------->
# FC for +NRG
rm(list=ls())

require(ggplot2)
source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type='L',ticks = "in",fontsize.cex = 1)

setwd("/Users/deepankar/Documents/Seafile/NGS-Data-DC/BaseSpace/20201125 Marika ERBB3 iSCREAM/Analysis/Figures/")
allMuts <- readRDS("20201216.ERBB3_iSCREAM_Mutations.RDS")

source('https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/MutSiteFind.R')
allMuts$AAPos=as.numeric(MutSiteFind(MutationColumn = allMuts$AAchange))
allMuts$FC_L <- allMuts$V2_NRG.VF/allMuts$Plasmid.VF
allMuts$FC_LPCR <- allMuts$V2_NRG.VF/allMuts$Plasmid_PCR.VF
allMuts$FC_WEHI <- allMuts$V2_NRG.VF/allMuts$V2_WEHI.VF

p <- ggplot(allMuts,aes(y=FC_L,x=AAPos,size=V2_NRG.VF,color=V2_NRG.VF))+geom_point(alpha=0.8)+scale_colour_gradient(low="#003366",high="firebrick1")+customtheme+ggtitle("Fold Change: +NRG (V2) vs Plasmid")+xlab("Amino Acid Position")+ylab("Fold Change")+scale_x_continuous(breaks = c(1,seq(50,1340,by = 50),1342))+geom_text(data=allMuts[allMuts$FC_L>10,],aes(label=allMuts[allMuts$FC_L>10,"AAchange"]),hjust=1.3,size=3)+geom_hline(yintercept = 0,color="black")+geom_hline(yintercept = 1,linetype="dashed",color="red",alpha=0.5);p
# ggsave(filename =  "FC_L.plot.svg",plot = p,width = 10,height = 5)
ggsave(filename =  "NRG/FC_L.pdf",plot = p,width = 10,height = 5)

p <- ggplot(allMuts,aes(y=FC_LPCR,x=AAPos,size=V2_NRG.VF,color=V2_NRG.VF))+geom_point(alpha=0.8)+scale_colour_gradient(low="#003366",high="firebrick1")+customtheme+ggtitle("Fold Change: +NRG (V2) vs Plasmid (PCR)")+xlab("Amino Acid Position")+ylab("Fold Change")+scale_x_continuous(breaks = c(1,seq(50,1340,by = 50),1342))+geom_text(data=allMuts[allMuts$FC_LPCR>10,],aes(label=allMuts[allMuts$FC_LPCR>10,"AAchange"]),hjust=1.3,size=3)+geom_hline(yintercept = 0,color="black")+geom_hline(yintercept = 1,linetype="dashed",color="red",alpha=0.5);p
# ggsave(filename =  "FC_LPCR.plot.svg",plot = p,width = 10,height = 5)
ggsave(filename =  "NRG/FC_LPCR.pdf",plot = p,width = 10,height = 5)

p <- ggplot(allMuts,aes(y=FC_WEHI,x=AAPos,size=V2_NRG.VF,color=V2_NRG.VF))+geom_point(alpha=0.8)+scale_colour_gradient(low="#003366",high="firebrick1")+customtheme+ggtitle("Fold Change: +NRG (V2) vs WEHI (V2)")+xlab("Amino Acid Position")+ylab("Fold Change")+scale_x_continuous(breaks = c(1,seq(50,1340,by = 50),1342))+geom_text(data=allMuts[allMuts$FC_WEHI>2.5,],aes(label=allMuts[allMuts$FC_WEHI>2.5,"AAchange"]),hjust=1.3,size=3)+geom_hline(yintercept = 0,color="black")+geom_hline(yintercept = 1,linetype="dashed",color="red",alpha=0.5);p
# ggsave(filename =  "FC_LPCR.plot.svg",plot = p,width = 10,height = 5)
ggsave(filename =  "NRG/FC_WEHI.pdf",plot = p,width = 10,height = 5)

## QC plot - FC calculations
value <- lm(allMuts$FC_L~allMuts$FC_LPCR)
anno_text <- paste("Slope:",signif(value$coefficients[[2]],2),"\nIntercept:",signif(value$coefficients[[1]],1));rm(value)

p <- ggplot(data = allMuts,aes(x=FC_L, y=FC_LPCR)) + geom_point(alpha=0.15) + geom_smooth(method = "lm", se = F, color="#f6546a") + xlab("Plasmid Library") + ylab("PCR from the Plasmid Library") + ggtitle("FoldChange (1 dot = 1 mutation) compared between samples") +  annotate("text",x=15,y=115,label=anno_text) + customtheme
ggsave(filename = "NRG/Cor_FC_Lib vs PCR.svg",plot = p,width = 5,height = 5);rm(p,anno_text)

value <- lm(allMuts$FC_L~allMuts$FC_WEHI)
anno_text <- paste("Slope:",signif(value$coefficients[[2]],2),"\nIntercept:",signif(value$coefficients[[1]],1));rm(value)

p <- ggplot(data = allMuts,aes(x=FC_L, y=FC_WEHI)) + geom_point(alpha=0.15) + geom_smooth(method = "lm", se = F, color="#f6546a") + xlab("Plasmid Library") + ylab("PCR from the Plasmid Library") + ggtitle("FoldChange (1 dot = 1 mutation) compared between samples") +  annotate("text",x=15,y=6.5,label=anno_text) + customtheme
ggsave(filename = "NRG/Cor_FC_Lib vs WEHI.svg",plot = p,width = 5,height = 5);rm(p,anno_text)

#---
# FC for +NRG -> -IL3
rm(list=ls())

require(ggplot2)
source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type='L',ticks = "in",fontsize.cex = 1)

setwd("/Users/deepankar/Documents/Seafile/NGS-Data-DC/BaseSpace/20201125 Marika ERBB3 iSCREAM/Analysis/Figures/")
allMuts <- readRDS("20201216.ERBB3_iSCREAM_Mutations.RDS")

source('https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/MutSiteFind.R')
allMuts$AAPos=as.numeric(MutSiteFind(MutationColumn = allMuts$AAchange))

allMuts$FC_L <- allMuts$V2_NRG2noIL3.VF/allMuts$Plasmid.VF
allMuts$FC_LPCR <- allMuts$V2_NRG2noIL3.VF/allMuts$Plasmid_PCR.VF
allMuts$FC_WEHI <- allMuts$V2_NRG2noIL3.VF/allMuts$V2_WEHI.VF
allMuts$FC_NRG <- allMuts$V2_NRG2noIL3.VF/allMuts$V2_NRG.VF

View(allMuts[,c("AAchange","FC_NRG")])

p <- ggplot(allMuts,aes(y=FC_L,x=AAPos,size=V2_NRG2noIL3.VF,color=V2_NRG2noIL3.VF))+geom_point(alpha=0.8)+scale_colour_gradient(low="#003366",high="firebrick1")+customtheme+ggtitle("Fold Change: +NRG -> -IL3 (V2) vs Plasmid")+xlab("Amino Acid Position")+ylab("Fold Change")+scale_x_continuous(breaks = c(1,seq(50,1340,by = 50),1342))+geom_text(data=allMuts[allMuts$FC_L>50,],aes(label=allMuts[allMuts$FC_L>50,"AAchange"]),hjust=1.3,size=3)+geom_hline(yintercept = 0,color="black")+geom_hline(yintercept = 1,linetype="dashed",color="red",alpha=0.5);p
ggsave(filename =  "NRG2noIL3/FC_noIL3_L.pdf",plot = p,width = 10,height = 5)

p <- ggplot(allMuts,aes(y=FC_LPCR,x=AAPos,size=V2_NRG2noIL3.VF,color=V2_NRG2noIL3.VF))+geom_point(alpha=0.8)+scale_colour_gradient(low="#003366",high="firebrick1")+customtheme+ggtitle("Fold Change: +NRG -> -IL3 (V2) vs Plasmid (PCR)")+xlab("Amino Acid Position")+ylab("Fold Change")+scale_x_continuous(breaks = c(1,seq(50,1340,by = 50),1342))+geom_text(data=allMuts[allMuts$FC_LPCR>50,],aes(label=allMuts[allMuts$FC_LPCR>50,"AAchange"]),hjust=1.3,size=3)+geom_hline(yintercept = 0,color="black")+geom_hline(yintercept = 1,linetype="dashed",color="red",alpha=0.5);p
ggsave(filename =  "NRG2noIL3/FC_noIL3_LPCR.pdf",plot = p,width = 10,height = 5)

p <- ggplot(allMuts,aes(y=FC_WEHI,x=AAPos,size=V2_NRG2noIL3.VF,color=V2_NRG2noIL3.VF))+geom_point(alpha=0.8)+scale_colour_gradient(low="#003366",high="firebrick1")+customtheme+ggtitle("Fold Change: +NRG -> -IL3 (V2) vs WEHI")+xlab("Amino Acid Position")+ylab("Fold Change")+scale_x_continuous(breaks = c(1,seq(50,1340,by = 50),1342))+geom_text(data=allMuts[allMuts$FC_WEHI>50,],aes(label=allMuts[allMuts$FC_WEHI>50,"AAchange"]),hjust=1.3,size=3)+geom_hline(yintercept = 0,color="black")+geom_hline(yintercept = 1,linetype="dashed",color="red",alpha=0.5);p
ggsave(filename =  "NRG2noIL3/FC_noIL3_WEHI.pdf",plot = p,width = 10,height = 5)

p <- ggplot(allMuts,aes(y=FC_NRG,x=AAPos,size=V2_NRG2noIL3.VF,color=V2_NRG2noIL3.VF))+geom_point(alpha=0.8)+scale_colour_gradient(low="#003366",high="firebrick1")+customtheme+ggtitle("Fold Change: +NRG -> -IL3 (V2) vs +NRG (V2)")+xlab("Amino Acid Position")+ylab("Fold Change")+scale_x_continuous(breaks = c(1,seq(50,1340,by = 50),1342))+geom_text(data=allMuts[allMuts$FC_NRG>50,],aes(label=allMuts[allMuts$FC_NRG>50,"AAchange"]),hjust=1.3,size=3)+geom_hline(yintercept = 0,color="black")+geom_hline(yintercept = 1,linetype="dashed",color="red",alpha=0.5);p
ggsave(filename =  "NRG2noIL3/FC_noIL3_NRG.pdf",plot = p,width = 10,height = 5)

#------->
# FC for - IL3
rm(list=ls())

require(ggplot2)
source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type='L',ticks = "in",fontsize.cex = 1)

setwd("/Users/deepankar/Documents/Seafile/NGS-Data-DC/BaseSpace/20201125 Marika ERBB3 iSCREAM/Analysis/Figures/")
allMuts <- readRDS("20201216.ERBB3_iSCREAM_Mutations.RDS")

source('https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/MutSiteFind.R')
allMuts$AAPos=as.numeric(MutSiteFind(MutationColumn = allMuts$AAchange))

allMuts$FC_L <- allMuts$V2_noIL3.VF/allMuts$Plasmid.VF
allMuts$FC_LPCR <- allMuts$V2_noIL3.VF/allMuts$Plasmid_PCR.VF
allMuts$FC_WEHI <- allMuts$V2_noIL3.VF/allMuts$V2_WEHI.VF

p <- ggplot(allMuts,aes(y=FC_L,x=AAPos,size=V2_noIL3.VF,color=V2_noIL3.VF))+geom_point(alpha=0.8)+scale_colour_gradient(low="#003366",high="firebrick1")+customtheme+ggtitle("Fold Change: -IL3 (V2) vs Plasmid")+xlab("Amino Acid Position")+ylab("Fold Change")+scale_x_continuous(breaks = c(1,seq(50,1340,by = 50),1342))+geom_text(data=allMuts[allMuts$FC_L>50,],aes(label=allMuts[allMuts$FC_L>50,"AAchange"]),hjust=1.3,size=3)+geom_hline(yintercept = 0,color="black")+geom_hline(yintercept = 1,linetype="dashed",color="red",alpha=0.5);p
ggsave(filename =  "noIL3/FC_noIL3_L.pdf",plot = p,width = 10,height = 5)

p <- ggplot(allMuts,aes(y=FC_LPCR,x=AAPos,size=V2_noIL3.VF,color=V2_noIL3.VF))+geom_point(alpha=0.8)+scale_colour_gradient(low="#003366",high="firebrick1")+customtheme+ggtitle("Fold Change: -IL3 (V2) vs Plasmid (PCR)")+xlab("Amino Acid Position")+ylab("Fold Change")+scale_x_continuous(breaks = c(1,seq(50,1340,by = 50),1342))+geom_text(data=allMuts[allMuts$FC_LPCR>50,],aes(label=allMuts[allMuts$FC_LPCR>50,"AAchange"]),hjust=1.3,size=3)+geom_hline(yintercept = 0,color="black")+geom_hline(yintercept = 1,linetype="dashed",color="red",alpha=0.5);p
ggsave(filename =  "noIL3/FC_noIL3_LPCR.pdf",plot = p,width = 10,height = 5)

p <- ggplot(allMuts,aes(y=FC_WEHI,x=AAPos,size=V2_noIL3.VF,color=V2_noIL3.VF))+geom_point(alpha=0.8)+scale_colour_gradient(low="#003366",high="firebrick1")+customtheme+ggtitle("Fold Change: -IL3 (V2) vs WEHI (V2)")+xlab("Amino Acid Position")+ylab("Fold Change")+scale_x_continuous(breaks = c(1,seq(50,1340,by = 50),1342))+geom_text(data=allMuts[allMuts$FC_WEHI>50,],aes(label=allMuts[allMuts$FC_WEHI>50,"AAchange"]),hjust=1.3,size=3)+geom_hline(yintercept = 0,color="black")+geom_hline(yintercept = 1,linetype="dashed",color="red",alpha=0.5);p
ggsave(filename =  "noIL3/FC_noIL3_WEHI.pdf",plot = p,width = 10,height = 5)

#--------->
# cDNA change Distribition Transitions v Transversions
rm(list=ls())

require(ggplot2)
source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type='L',ticks = "in",fontsize.cex = 1.2)

setwd("/Users/deepankar/Documents/Seafile/NGS-Data-DC/BaseSpace/20201125 Marika ERBB3 iSCREAM/Analysis/Figures/")
allMuts <- readRDS("20201216.ERBB3_iSCREAM_Mutations.RDS")

allMuts=na.exclude(allMuts)

bases=c("A","C","G","T")
mut.mat=matrix(data=NA, nrow = 4, ncol = 4)
colnames(mut.mat)=row.names(mut.mat)=bases
rm(bases)
subDF=allMuts[,c("AAchange","Ref","Alt")]
colnames(subDF)=c("AAchange","ref","alt")

for(i in bases){
  temp.i=subset(subDF,subDF$ref==i)
  for(j in unique(temp.i$alt)){
    temp.j=subset(temp.i,temp.i$alt==j)
    nobs=dim(temp.j)[1]
    mut.mat[i,j]=nobs
    print(paste(i,j,nobs))
  }
}
rm(i,j,nobs,temp.i,temp.j)

breakList=seq(1, max(mut.mat,na.rm = T), by = 5)
colours=colorRampPalette(c("navyblue","maroon4","red1"))(length(breakList))

library(pheatmap)
pheatmap(t(mut.mat),color=colours,cluster_rows = F,cluster_cols = F,border_color = "black",cellwidth = 25,cellheight = 25,display_numbers = T,number_color = "white",number_format="%.0f",main="ERBB4 iSCREAM DNA change Distribution",filename = "DNAchange distrib.pdf")
rm(breakList,colours)

#-----------
# Table export KE KK MK
# FC for +NRG -> -IL3 vs NRG
# FC for -IL3 vs WEHI
rm(list=ls())

require(ggplot2)
source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type='L',ticks = "in",fontsize.cex = 1)

setwd("/Users/deepankar/Documents/Seafile/NGS-Data-DC/BaseSpace/20201125 Marika ERBB3 iSCREAM/Analysis/Figures/")
allMuts <- readRDS("20201216.ERBB3_iSCREAM_Mutations.RDS")

source('https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/MutSiteFind.R')
allMuts$AAPos=as.numeric(MutSiteFind(MutationColumn = allMuts$AAchange))


allMuts$FC_NRG2noIL3_vs_NRG <- allMuts$V2_NRG2noIL3.VF/allMuts$V2_NRG.VF
allMuts$FC_noIL3_vs_WEHI <- allMuts$V2_noIL3.VF/allMuts$V2_WEHI.VF
allMuts[is.na(allMuts)] <- 0
write.table(x = allMuts,
            file = "~/Desktop/20201216.ERBB3_iSCREAM_Mutations.tsv",
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = T)

