#----<---->----

###---- Marika asked on 2022.01.19
# Could you prepare an alternative figure for the Illumina FC results 
# so that the +NRG to no ligand results from vol1 + vol2 would be merged together? 
# We could possibly go with this and explain in the figure text that 
# the results are from 2 different experiments. 

#----<---->----
library(data.table)
options(scipen = 100, digits = 4)
#<--------
# FC for +NRG
rm(list=ls())

require(ggplot2)
source('https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/MutSiteFind.R')
source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type='L',
                                  ticks = "in",
                                  fontsize.cex = 1)

setwd("/Users/deepankar/Documents/Seafile/NGS-Data-DC/BaseSpace/20201125 Marika ERBB3 iSCREAM/Analysis/Figures/")
V1 <- as.data.table(readRDS("V1/20210108.ERBB3_iSCREAM_V1_Mutations.RDS"))
V2 <- as.data.table(readRDS("V2/20201216.ERBB3_iSCREAM_V2_Mutations.RDS"))

V1[, FC_L := V1_NRG.VF / Plasmid.VF]
V1[, FC_LPCR := V1_NRG.VF / Plasmid_PCR.VF]
V1[, FC_WEHI := V1_NRG.VF / V1_WEHI.VF]

V2[, FC_L := V2_NRG.VF / Plasmid.VF]
V2[, FC_LPCR := V2_NRG.VF / Plasmid_PCR.VF]
V2[, FC_WEHI := V2_NRG.VF / V2_WEHI.VF]

objects_to_retain <- c("objects_to_retain",ls())

# Plasmid Library
rm(list = ls()[!ls() %in% objects_to_retain])
plotMuts <- data.table(MutID = unique(c(V1$MutID, V2$MutID)))
plotMuts[, V1 := V1$FC_L[match(plotMuts$MutID, V1$MutID)]]
plotMuts[, V2 := V2$FC_L[match(plotMuts$MutID, V2$MutID)]]
plotMuts[, FC_L := max(V1, V2, na.rm = T), by = 1:nrow(plotMuts)]

plotMuts$V1_NRG <- V1$V1_NRG.VF[match(plotMuts$MutID, V1$MutID)]
plotMuts$V2_NRG <-  V2$V2_NRG.VF[match(plotMuts$MutID, V2$MutID)]
plotMuts[, NRG.VF := max(V1_NRG, V2_NRG, na.rm = T), by = 1:nrow(plotMuts)]

plotMuts$AAchange <- V2$AAchange[match(plotMuts$MutID, V2$MutID)]
missing_index <- which(is.na(plotMuts$AAchange))
plotMuts$AAchange[missing_index] <-
  V1$AAchange[match(plotMuts$MutID[missing_index], V1$MutID)]

plotMuts <- data.table:::na.omit.data.table(plotMuts)
plotMuts[, colnames(plotMuts)[grep(pattern = "V[12].*", x = colnames(plotMuts))] := NULL]

plotMuts[, AAPos := as.numeric(MutSiteFind(plotMuts$AAchange))]

label_threshold <- 20
p <- ggplot(plotMuts,aes(y=FC_L,x=AAPos,size=NRG.VF,color=NRG.VF))+
  geom_point(alpha=0.8)+
  scale_colour_gradient(low="#003366",high="firebrick1")+
  customtheme+
  ggtitle(paste0("Fold Change: +NRG (V1 & V2) vs Plasmid (n = ",nrow(plotMuts),")"))+
  xlab("Amino Acid Position")+
  ylab("Fold Change")+
  scale_x_continuous(breaks = c(1,seq(50,1340,by = 50),1342))+
  geom_text(data= plotMuts[plotMuts$FC_L > label_threshold, ],
            aes(x = AAPos, 
                y = FC_L,
                label = plotMuts$AAchange[plotMuts$FC_L > label_threshold]),
                hjust=1.3,
                size=3)+
  geom_hline(yintercept = 0,color="black")+
  geom_hline(yintercept = 1,linetype="dashed",color="red",alpha=0.5);p
# ggsave(filename =  "FC_L.plot.svg",plot = p,width = 10,height = 5)
ggsave(filename =  "NRG_FC_Plasmid.pdf",plot = p,width = 10,height = 5)

# Plasmid_PCR Library
rm(list = ls()[!ls() %in% objects_to_retain])

plotMuts <- data.table(MutID = unique(c(V1$MutID, V2$MutID)))
plotMuts[, V1 := V1$FC_LPCR[match(plotMuts$MutID, V1$MutID)]]
plotMuts[, V2 := V2$FC_LPCR[match(plotMuts$MutID, V2$MutID)]]
plotMuts[, FC_LPCR := max(V1, V2, na.rm = T), by = 1:nrow(plotMuts)]

plotMuts$V1_NRG <- V1$V1_NRG.VF[match(plotMuts$MutID, V1$MutID)]
plotMuts$V2_NRG <-  V2$V2_NRG.VF[match(plotMuts$MutID, V2$MutID)]
plotMuts[, NRG.VF := max(V1_NRG, V2_NRG, na.rm = T), by = 1:nrow(plotMuts)]

plotMuts$AAchange <- V2$AAchange[match(plotMuts$MutID, V2$MutID)]
missing_index <- which(is.na(plotMuts$AAchange))
plotMuts$AAchange[missing_index] <-
  V1$AAchange[match(plotMuts$MutID[missing_index], V1$MutID)]

plotMuts <- data.table:::na.omit.data.table(plotMuts)
plotMuts[, colnames(plotMuts)[grep(pattern = "V[12].*", x = colnames(plotMuts))] := NULL]

plotMuts[, AAPos := as.numeric(MutSiteFind(plotMuts$AAchange))]

label_threshold <- 20
p <- ggplot(plotMuts,aes(y=FC_LPCR,x=AAPos,size=NRG.VF,color=NRG.VF))+
  geom_point(alpha=0.8)+
  scale_colour_gradient(low="#003366",high="firebrick1")+
  customtheme+
  ggtitle(paste0("Fold Change: +NRG (V1 & V2) vs PCR from Plasmid (n = ",nrow(plotMuts),")"))+
  xlab("Amino Acid Position")+
  ylab("Fold Change")+
  scale_x_continuous(breaks = c(1,seq(50,1340,by = 50),1342))+
  geom_text(data= plotMuts[plotMuts$FC_LPCR > label_threshold, ],
            aes(x = AAPos, 
                y = FC_LPCR,
                label = plotMuts$AAchange[plotMuts$FC_LPCR > label_threshold]),
            hjust=1.3,
            size=3)+
  geom_hline(yintercept = 0,color="black")+
  geom_hline(yintercept = 1,linetype="dashed",color="red",alpha=0.5);p
# ggsave(filename =  "FC_LPCR.plot.svg",plot = p,width = 10,height = 5)
ggsave(filename =  "NRG_FC_Plasmid_PCR.pdf",plot = p,width = 10,height = 5)


# WEHI
rm(list = ls()[!ls() %in% objects_to_retain])

plotMuts <- data.table(MutID = unique(c(V1$MutID, V2$MutID)))
plotMuts[, V1 := V1$FC_WEHI[match(plotMuts$MutID, V1$MutID)]]
plotMuts[, V2 := V2$FC_WEHI[match(plotMuts$MutID, V2$MutID)]]
plotMuts[, FC_WEHI := max(V1, V2, na.rm = T), by = 1:nrow(plotMuts)]

plotMuts$V1_NRG <- V1$V1_NRG.VF[match(plotMuts$MutID, V1$MutID)]
plotMuts$V2_NRG <-  V2$V2_NRG.VF[match(plotMuts$MutID, V2$MutID)]
plotMuts[, NRG.VF := max(V1_NRG, V2_NRG, na.rm = T), by = 1:nrow(plotMuts)]

plotMuts$AAchange <- V2$AAchange[match(plotMuts$MutID, V2$MutID)]
missing_index <- which(is.na(plotMuts$AAchange))
plotMuts$AAchange[missing_index] <-
  V1$AAchange[match(plotMuts$MutID[missing_index], V1$MutID)]

plotMuts <- data.table:::na.omit.data.table(plotMuts)
plotMuts[, colnames(plotMuts)[grep(pattern = "V[12].*", x = colnames(plotMuts))] := NULL]

plotMuts[, AAPos := as.numeric(MutSiteFind(plotMuts$AAchange))]

label_threshold <- 3
p <- ggplot(plotMuts,aes(y=FC_WEHI,x=AAPos,size=NRG.VF,color=NRG.VF))+
  geom_point(alpha=0.8)+
  scale_colour_gradient(low="#003366",high="firebrick1")+
  customtheme+
  ggtitle(paste0("Fold Change: +NRG (V1 & V2) vs WEHI (n = ",nrow(plotMuts),")"))+
  xlab("Amino Acid Position")+
  ylab("Fold Change")+
  scale_x_continuous(breaks = c(1,seq(50,1340,by = 50),1342))+
  geom_text(data= plotMuts[plotMuts$FC_WEHI > label_threshold, ],
            aes(x = AAPos, 
                y = FC_WEHI,
                label = plotMuts$AAchange[plotMuts$FC_WEHI > label_threshold]),
            hjust=1.3,
            size=3)+
  geom_hline(yintercept = 0,color="black")+
  geom_hline(yintercept = 1,linetype="dashed",color="red",alpha=0.5);p
# ggsave(filename =  "FC_WEHI.plot.svg",plot = p,width = 10,height = 5)
ggsave(filename =  "NRG_FC_WEHI.pdf",plot = p,width = 10,height = 5)
