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
V1 <- fread("")

V1 <- as.data.table(readRDS("V1/20210108.ERBB3_iSCREAM_V1_Mutations.RDS"))
V2 <- as.data.table(readRDS("V2/20201216.ERBB3_iSCREAM_V2_Mutations.RDS"))

V1[, FC :=  V1_NRG2noIL3.VF / V1_NRG.VF]
V2[, FC :=  V2_NRG2noIL3.VF / V2_NRG.VF]

objects_to_retain <- c("objects_to_retain",ls())

# Plasmid Library
rm(list = ls()[!ls() %in% objects_to_retain])
plotMuts <- data.table(MutID = unique(c(V1$MutID, V2$MutID)))
plotMuts[, V1 := V1$FC[match(plotMuts$MutID, V1$MutID)]]
plotMuts[, V2 := V2$FC[match(plotMuts$MutID, V2$MutID)]]
plotMuts[, FC := max(V1, V2, na.rm = T), by = 1:nrow(plotMuts)]

plotMuts$V1_VF <- V1$V1_NRG2noIL3.VF[match(plotMuts$MutID, V1$MutID)]
plotMuts$V2_VF <-  V2$V2_NRG2noIL3.VF[match(plotMuts$MutID, V2$MutID)]
plotMuts[, VF := max(V1_VF, V2_VF, na.rm = T), by = 1:nrow(plotMuts)]

plotMuts$AAchange <- V2$AAchange[match(plotMuts$MutID, V2$MutID)]
missing_index <- which(is.na(plotMuts$AAchange))
plotMuts$AAchange[missing_index] <-
  V1$AAchange[match(plotMuts$MutID[missing_index], V1$MutID)]

plotMuts <- data.table:::na.omit.data.table(plotMuts)
plotMuts[, colnames(plotMuts)[grep(pattern = "V[12].*", x = colnames(plotMuts))] := NULL]

plotMuts[, AAPos := as.numeric(MutSiteFind(plotMuts$AAchange))]

label_threshold <- 20
p <- ggplot(plotMuts,aes(y=FC,x=AAPos,size=VF,color=VF))+
  geom_point(alpha=0.8)+
  scale_colour_gradient(low="#003366",high="firebrick1")+
  customtheme+
  ggtitle(paste0("Fold Change: +NRG_to_NoIL3 (V1 & V2) vs NRG (n = ",nrow(plotMuts),")"))+
  xlab("Amino Acid Position")+
  ylab("Fold Change")+
  scale_x_continuous(breaks = c(1,seq(50,1340,by = 50),1342))+
  geom_text(data= plotMuts[plotMuts$FC > label_threshold, ],
            aes(x = AAPos, 
                y = FC,
                label = plotMuts$AAchange[plotMuts$FC > label_threshold]),
                hjust=1.3,
                size=3)+
  geom_hline(yintercept = 0,color="black")+
  geom_hline(yintercept = 1,linetype="dashed",color="red",alpha=0.5);p
# ggsave(filename =  "FC.plot.svg",plot = p,width = 10,height = 5)
ggsave(filename =  "NRG2noIL3_FC_NRG.pdf",plot = p,width = 10,height = 5)

