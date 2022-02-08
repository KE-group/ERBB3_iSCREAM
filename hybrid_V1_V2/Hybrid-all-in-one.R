#----<---->----

###---- Marika asked on 2022.02.07
# Prepare a figure with WEHI -> No IL3 (V2) and add in data from 
# WEHI -> NRG2NoIL3 (V1) + WEHI -> NRG2NoIL3 (V2)
# Do correlations of FC between WEHI -> NRG (V1) and (V2)

#---- <----> ----
# Differences in Depth between samples creates a doubt in my mind 
# about merging VFs of mutations from different samples into a single scale.

setwd("/Volumes/DATA/Seafile/NGS-Data-DC/BaseSpace/20201125 Marika ERBB3 iSCREAM/Analysis/stats/depth")
averageDepth <- function(FILE, FUN){
  dat <- data.table::fread(FILE)
  return(FUN(dat$V3,na.rm = T))
}
list_of_files <-
  c("V2_noIL3.dep.txt",
    "V1_NRG2noIL3.dep.txt",
    "V2_NRG2noIL3.dep.txt")

depth <-
  parallel::mclapply(
    X = list_of_files,
    FUN = function(X)
      averageDepth(X, median), # change mean or median here
    mc.cores = 2
  )
names(depth) <- list_of_files
print(format(depth,nsmall = 0, big.mark = ","))
rm(list=ls())

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

# Correlation <-------
rm(list = ls()[!ls() %in% objects_to_retain]); gc()

V1 <- V1.bak; V2 <- V2.bak
V1[, FC :=  V1_NRG.VF / V1_WEHI.VF]
V2[, FC :=  V2_NRG.VF / V2_WEHI.VF]
plotMuts <- data.table(MutID = unique(c(V1$MutID, V2$MutID)))
plotMuts[, V1 := V1$FC[match(plotMuts$MutID, V1$MutID)]]
plotMuts[, V2 := V2$FC[match(plotMuts$MutID, V2$MutID)]]
plotMuts <- na.exclude(plotMuts)

customtheme <- DC_theme_generator(type='L',
                                  ticks = "out",
                                  fontsize.cex = 1,
                                  title.fontsize.cex = 1)

# r <- cor.test(plotMuts$V1, plotMuts$V2, method = "spearman")
r <- cor.test(plotMuts$V1, plotMuts$V2, method = "pearson")
correlation_text <-
  paste0("r=", round(unname(r$estimate), digits = 3), ", P < 2e-16")

ggplot(data = plotMuts,aes(x=V1,y=V2))+
  geom_point(alpha=0.25,size=2)+
  geom_smooth(method = "lm", 
              se = F, 
              na.rm = T, 
              color = "red",
              alpha = 0.25)+
  scale_y_continuous(expand=c(0,0),limits = c(0,10))+
  scale_x_continuous(expand=c(0,0),limits = c(0,10))+
  ggtitle("Correlation: fold change(WEHI to NRG) in V1 vs V2")+
  xlab("fold change - V1")+
  ylab("fold change - V2")+
  annotate(geom = "text",
           x = 1.5,
           y = 8.5,
           label = correlation_text)+
  customtheme

# ggsave(
#   "correlation.pdf",
#   device = cairo_pdf,
#   height = 5,
#   width = 5
# )

#---- <----> ----
# Merge V1 & V2 of NRG vs NRG -> No IL3
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

plotMuts[, FC := max(V1, V2, V2_noIL3, na.rm = T), by = 1:nrow(plotMuts)]

plotMuts$V1_VF <- V1$V1_NRG2noIL3.VF[match(plotMuts$MutID, V1$MutID)]
plotMuts$V2_VF <-  V2$V2_NRG2noIL3.VF[match(plotMuts$MutID, V2$MutID)]
plotMuts$V2_noIL3_VF <-  V2$V2_noIL3.VF[match(plotMuts$MutID, V2$MutID)]
plotMuts[, VF := max(V1_VF, V2_VF, V2_noIL3_VF, na.rm = T), by = 1:nrow(plotMuts)]

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
# Fitting distribution and finding critical value

library(fitdistrplus)
set.seed(2022)
plotMuts[,log_FC := log2(FC)]
p <- ggplot(data = plotMuts,aes(x=log_FC))+
  geom_density()+ 
  geom_vline(xintercept = mean(plotMuts$log_FC),
               linetype = "dashed",
               color="black")+
  scale_y_continuous(limits = c(0,1))+
  customtheme

p <- p + geom_density(data = plotMuts, 
                  aes(x = log2(V1)),
                  fill="#0057e7",
                  alpha=0.3,
                  color=NA)+
  geom_density(data = plotMuts, 
               aes(x = log2(V2_noIL3)),
               fill="#008744",
               alpha=0.3,
               color=NA)+
  geom_density(data = plotMuts, 
               aes(x = log2(V2)),
               fill="#d62d20",
               alpha=0.3,
               color=NA)+
  scale_y_continuous(limits = c(0,0.5))

ggsave(filename = "Distribution problems/ERBB3_individual.png",
       plot = p, 
       height = 4,
       width = 5)


Mixture <- data.frame(model1=rep(NA, 50000),model2=NA)
set.seed(2022)
Mixture$model1 <- fGarch::rsnorm(n = 50000, mean = -2.5, sd = 2.5, xi = -1.1 )
Mixture$model2 <- fGarch::rsnorm(n = 50000, mean = -0.1, sd = 0.75 , xi = 1)

p + geom_density(data = Mixture, 
                 aes(x = model1),
                 fill="#0000FF",
                 alpha=0.3,
                 color=NA)+
  geom_density(data = Mixture, 
               aes(x = model2),
               fill="#FF0000",
               alpha=0.3,
               color=NA)+
  geom_vline(xintercept = -0.1,
             linetype = "dotted",
             color="red")+
  geom_vline(xintercept = -2.5,
             linetype = "dotted",
             color="blue")

#-----------

normfit <- fitdist(plotMuts$FC, "norm", method = 'mle') # 'mle' default
gammafit<- fitdist(plotMuts$FC, "gamma")
weibullfit <- fitdist(plotMuts$FC, "weibull")
lnormfit <- fitdist(plotMuts$FC, "lnorm")
p <-  ppcomp(list(gammafit,weibullfit,lnormfit,normfit),plotstyle = "ggplot")
p+ customtheme

normfit <- fitdist(data = plotMuts$log_FC, distr = "norm", method = 'mle') # 'mle' default
summary(normfit)
p <-
  denscomp(
    list(normfit),
    demp = T,
    plotstyle = "ggplot",
    ylim = c(0, 1)
  )
p+customtheme

X=sort(plotMuts$log_FC)
Y=dnorm(x = X,mean = normfit$estimate[['mean']],sd = normfit$estimate[['sd']])
plot(x = X,y = Y,xlab="X",ylab="dnorm(x)",pch=20,col="#00000011")

Y=2*pnorm(q = -abs(X),mean = normfit$estimate[['mean']],sd = normfit$estimate[['sd']])
plot(x = X,y = Y,xlab="X",ylab="2*pnorm(-abs(x))",pch=20,col="#00000011")


# lnormfit<- fitdist(plotMuts$FC, "lnorm")
# p <- denscomp(list(lnormfit),demp = T,plotstyle = "ggplot",ylim=c(0,1))
# p+customtheme

plotMuts$P.value <- 2*pnorm(q= -(abs(plotMuts$log_FC)),mean = normfit$estimate[['mean']],sd = normfit$estimate[['sd']])
plotMuts$FDR <- p.adjust(p = plotMuts$P.value, method = "fdr")

# Calculating threshold and selecting mutations above threshold
(critical.value <-
    qnorm(
      p = 5.613646e-07,
      mean = normfit$estimate[['mean']],
      sd = normfit$estimate[['sd']]))
FC.sub <- subset(plotMuts, plotMuts$FDR < 0.05)
# FC.sub <- subset(plotMuts, plotMuts$FDR < 0.01 & plotMuts$FC > 1)


#---- <----> ----
# Drawing Fold chane scatter plot
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
ggsave(filename =  "Hybrid-all-in-one-fold-change.pdf",plot = p,width = 10,height = 5)

