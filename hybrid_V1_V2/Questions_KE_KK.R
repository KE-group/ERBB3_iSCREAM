#---- <----> ----
library(data.table)
setDTthreads(1)
options(scipen = 100, digits = 4)

# Read Data <-------
rm(list=ls())

require(ggplot2)
source('https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/MutSiteFind.R')
source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')

setwd("/Users/deepankar/Documents/Seafile/NGS-Data-DC/BaseSpace/20201125 Marika ERBB3 iSCREAM/Analysis/Figures/")

V1.bak <- as.data.table(readRDS("V1/20210108.ERBB3_iSCREAM_V1_Mutations.RDS"))
V2.bak <- as.data.table(readRDS("V2/20201216.ERBB3_iSCREAM_V2_Mutations.RDS"))
objects_to_retain <- c("objects_to_retain",ls())

#---- <----> ----
# Check inverse fold Change and enriched muts
rm(list = ls()[!ls() %in% objects_to_retain]);gc()
V1 <- V1.bak; V2 <- V2.bak

V1[, FC :=  V1_NRG2noIL3.VF / V1_WEHI.VF]
V2[, FC :=  V2_NRG2noIL3.VF / V2_WEHI.VF]
V2[, FC_noIL3 :=  V2_noIL3.VF / V2_WEHI.VF]

plotMuts <- data.table(MutID = unique(c(V1$MutID, V2$MutID)))
plotMuts[, V1 := V1$FC[match(plotMuts$MutID, V1$MutID)]]
plotMuts[, V2 := V2$FC[match(plotMuts$MutID, V2$MutID)]]
# plotMuts[, V2_noIL3 := V2$FC_noIL3[match(plotMuts$MutID, V2$MutID)]] # gives error
plotMuts$V2_noIL3 <-  V2$FC_noIL3[match(plotMuts$MutID, V2$MutID)]

plotMuts[, FC := min(V1, V2, V2_noIL3, na.rm = T), by = 1:nrow(plotMuts)]
plotMuts[, FC := (-1*(1/FC))]

plotMuts$V1_VF <- V1$V1_NRG2noIL3.VF[match(plotMuts$MutID, V1$MutID)]
plotMuts$V2_VF <-  V2$V2_NRG2noIL3.VF[match(plotMuts$MutID, V2$MutID)]
plotMuts$V2_noIL3_VF <-  V2$V2_noIL3.VF[match(plotMuts$MutID, V2$MutID)]
plotMuts[, VF := min(V1_VF, V2_VF, V2_noIL3_VF, na.rm = T), by = 1:nrow(plotMuts)]

plotMuts$AAchange <- V2$AAchange[match(plotMuts$MutID, V2$MutID)]
missing_index <- which(is.na(plotMuts$AAchange))
plotMuts$AAchange[missing_index] <-
  V1$AAchange[match(plotMuts$MutID[missing_index], V1$MutID)]

plotMuts <- data.table:::na.omit.data.table(plotMuts)
plotMuts[, colnames(plotMuts)[grep(pattern = "V[12].*", x = colnames(plotMuts))] := NULL]

plotMuts[, AAPos := as.numeric(MutSiteFind(plotMuts$AAchange))]
customtheme <- DC_theme_generator(type='L',
                                  ticks = "in",
                                  fontsize.cex = 1,
                                  title.fontsize.cex = 1)

#---- <----> ----
# Drawing Inverse Fold change scatter plot
label_threshold <- -150
p <- ggplot(plotMuts,aes(y=FC,x=AAPos,size=VF,color=VF))+
  geom_point(alpha=0.8)+
  scale_colour_gradient(low="#003366",high="firebrick1")+
  customtheme+
  ggtitle(paste0("Fold Change: +NRG_to_NoIL3 (V1 & V2) vs NRG (n = ",nrow(plotMuts),")"))+
  xlab("ERBB3 primary sequence")+
  ylab("Fold Change\n(most depleted mutations)")+
  scale_x_continuous(breaks = c(1,seq(50,1340,by = 50),1342))+
  geom_text(data= plotMuts[plotMuts$FC < label_threshold, ],
            aes(x = AAPos, 
                y = FC,
                label = plotMuts$AAchange[plotMuts$FC < label_threshold]),
            hjust=1.3,
            size=3)+
  geom_hline(yintercept = 0,color="black");p
# ggsave(filename =  "FC.plot.svg",plot = p,width = 10,height = 5)
ggsave(filename =  "Inverse-fold-change.pdf",plot = p,width = 10,height = 5)

label_threshold <- -100
p <- ggplot(plotMuts,aes(y=FC,x=AAPos,size=VF,color=VF))+
  geom_point(alpha=0.8)+
  scale_colour_gradient(low="#003366",high="firebrick1")+
  customtheme+
  ggtitle(paste0("Fold Change: +NRG_to_NoIL3 (V1 & V2) vs NRG (n = ",nrow(plotMuts),")"))+
  xlab("ERBB3 primary sequence")+
  ylab("Fold Change\n(most depleted mutations)")+
  scale_x_continuous(breaks = c(1,seq(50,1340,by = 50),1342))+
  geom_text(data= plotMuts[plotMuts$FC < label_threshold, ],
            aes(x = AAPos, 
                y = FC,
                label = plotMuts$AAchange[plotMuts$FC < label_threshold]),
            hjust=1.3,
            size=3)+
  geom_hline(yintercept = 0,color="black")+
  scale_y_continuous(limits = c(-500,0));p
# ggsave(filename =  "FC.plot.svg",plot = p,width = 10,height = 5)
ggsave(filename =  "Inverse-fold-change_zoom.pdf",plot = p,width = 10,height = 5)


#---- <----> ----
# Check PacBio Data
#---- <----> ----
library(data.table)
setDTthreads(1)
options(scipen = 100, digits = 4)

# Read Data <-------
rm(list=ls())

require(ggplot2)
source('https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/MutSiteFind.R')
source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')

setwd("/Users/deepankar/Documents/Seafile/NGS-Data-DC/PacBio/20210909 iSCREAM Run1/analysis/")

V1.bak <- as.data.table(readRDS("filtered/R.Binary/ERBB3_V1_iSCREAM.RDS"))
V2.bak <- as.data.table(readRDS("filtered/R.Binary/ERBB3_V2_iSCREAM.RDS"))
objects_to_retain <- c("objects_to_retain",ls())

#---- <----> ----
# Check inverse fold Change and enriched muts
rm(list = ls()[!ls() %in% objects_to_retain]);gc()
V1 <- V1.bak; V2 <- V2.bak

V1[, FC :=  V1_NRG2noIL3.VF / V1_WEHI.VF]
V2[, FC :=  V2_NRG2noIL3.VF / V2_WEHI.VF]
V2[, FC_noIL3 :=  V2_noIL3.VF / V2_WEHI.VF]

plotMuts <- data.table(MutID = unique(c(V1$MutID, V2$MutID)))
plotMuts[, V1 := V1$FC[match(plotMuts$MutID, V1$MutID)]]
plotMuts[, V2 := V2$FC[match(plotMuts$MutID, V2$MutID)]]
# plotMuts[, V2_noIL3 := V2$FC_noIL3[match(plotMuts$MutID, V2$MutID)]] # gives error
plotMuts$V2_noIL3 <-  V2$FC_noIL3[match(plotMuts$MutID, V2$MutID)]

plotMuts[, FC := min(V1, V2, V2_noIL3, na.rm = T), by = 1:nrow(plotMuts)]
plotMuts[, FC := (-1*(1/FC))]

plotMuts$V1_VF <- V1$V1_NRG2noIL3.VF[match(plotMuts$MutID, V1$MutID)]
plotMuts$V2_VF <-  V2$V2_NRG2noIL3.VF[match(plotMuts$MutID, V2$MutID)]
plotMuts$V2_noIL3_VF <-  V2$V2_noIL3.VF[match(plotMuts$MutID, V2$MutID)]
plotMuts[, VF := min(V1_VF, V2_VF, V2_noIL3_VF, na.rm = T), by = 1:nrow(plotMuts)]

plotMuts$AAchange <- V2$AAchange[match(plotMuts$MutID, V2$MutID)]
missing_index <- which(is.na(plotMuts$AAchange))
plotMuts$AAchange[missing_index] <-
  V1$AAchange[match(plotMuts$MutID[missing_index], V1$MutID)]

plotMuts <- data.table:::na.omit.data.table(plotMuts)
plotMuts[, colnames(plotMuts)[grep(pattern = "V[12].*", x = colnames(plotMuts))] := NULL]

plotMuts[, AAPos := as.numeric(MutSiteFind(plotMuts$AAchange))]
customtheme <- DC_theme_generator(type='L',
                                  ticks = "in",
                                  fontsize.cex = 1,
                                  title.fontsize.cex = 1)

#---- <----> ----
# Drawing Inverse Fold change scatter plot
label_threshold <- -150
p <- ggplot(plotMuts,aes(y=FC,x=AAPos,size=VF,color=VF))+
  geom_point(alpha=0.8)+
  scale_colour_gradient(low="#003366",high="firebrick1")+
  customtheme+
  ggtitle(paste0("Fold Change: +NRG_to_NoIL3 (V1 & V2) vs NRG (n = ",nrow(plotMuts),")"))+
  xlab("ERBB3 primary sequence")+
  ylab("Fold Change\n(most depleted mutations)")+
  scale_x_continuous(breaks = c(1,seq(50,1340,by = 50),1342))+
  geom_text(data= plotMuts[plotMuts$FC < label_threshold, ],
            aes(x = AAPos, 
                y = FC,
                label = plotMuts$AAchange[plotMuts$FC < label_threshold]),
            hjust=1.3,
            size=3)+
  geom_hline(yintercept = 0,color="black");p
# ggsave(filename =  "FC.plot.svg",plot = p,width = 10,height = 5)
ggsave(filename =  "Inverse-fold-change_PACBIO.pdf",plot = p,width = 10,height = 5)

label_threshold <- -100
p <- ggplot(plotMuts,aes(y=FC,x=AAPos,size=VF,color=VF))+
  geom_point(alpha=0.8)+
  scale_colour_gradient(low="#003366",high="firebrick1")+
  customtheme+
  ggtitle(paste0("Fold Change: +NRG_to_NoIL3 (V1 & V2) vs NRG (n = ",nrow(plotMuts),")"))+
  xlab("ERBB3 primary sequence")+
  ylab("Fold Change\n(most depleted mutations)")+
  scale_x_continuous(breaks = c(1,seq(50,1340,by = 50),1342))+
  geom_text(data= plotMuts[plotMuts$FC < label_threshold, ],
            aes(x = AAPos, 
                y = FC,
                label = plotMuts$AAchange[plotMuts$FC < label_threshold]),
            hjust=1.3,
            size=3)+
  geom_hline(yintercept = 0,color="black")+
  scale_y_continuous(limits = c(-500,0));p
# ggsave(filename =  "FC.plot.svg",plot = p,width = 10,height = 5)
ggsave(filename =  "Inverse-fold-change_PACBIO_zoom.pdf",plot = p,width = 10,height = 5)

