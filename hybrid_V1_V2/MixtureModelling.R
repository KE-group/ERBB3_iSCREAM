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

source("https://github.com/dchakro/shared_Rscripts/raw/master/Annovar/Annovar_cDNA_Find.R")

V1.bak <- as.data.table(readRDS("V1/20210108.ERBB3_iSCREAM_V1_Mutations.RDS"))
V1.bak$cDNAchange=annovar_cDNA_Find(V1.bak$AAChange.refGene,
                                    isoform = "NM_001982") 
V2.bak <- as.data.table(readRDS("V2/20201216.ERBB3_iSCREAM_V2_Mutations.RDS"))
V2.bak$cDNAchange=annovar_cDNA_Find(V2.bak$AAChange.refGene,
                                    isoform = "NM_001982") 
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

plotMuts$cDNAchange <- V2$cDNAchange[match(plotMuts$MutID, V2$MutID)]
missing_index <- which(is.na(plotMuts$cDNAchange))
plotMuts$cDNAchange[missing_index] <-
  V1$cDNAchange[match(plotMuts$MutID[missing_index], V1$MutID)]

plotMuts <- data.table:::na.omit.data.table(plotMuts)
# plotMuts[, colnames(plotMuts)[grep(pattern = "V[12].*", x = colnames(plotMuts))] := NULL]

plotMuts[, AAPos := as.numeric(MutSiteFind(plotMuts$AAchange))]
plotMuts[, cDNAPos := as.numeric(MutSiteFind(plotMuts$cDNAchange))]

# saveRDS(object = plotMuts,file = "20220506.plotMuts.RDS")

customtheme <- DC_theme_generator(type='L',
                                  ticks = "in",
                                  fontsize.cex = 1,
                                  title.fontsize.cex = 1)

#---- <----> ----
# Fitting distribution and finding critical value

library(fitdistrplus)
set.seed(2022)
plotMuts[,log_FC := log2(FC)]
plotMuts.bak <- plotMuts
library(mixtools)

model <- normalmixEM(plotMuts$log_FC)
mixtools::plot.mixEM(
  model,
  which = 2,
  density = TRUE,
  cex.axis = 1.4,
  cex.lab = 1.4,
  cex.main = 1.8,
  main2 = "Gaussian mixture",
  xlab2 = "log2(fold change)"
)
lines(density(plotMuts$log_FC),lty=2,lwd=2)
summary(model)
abline(v = c(4,6),lty=2) 

# Sampling from fitted distribution
probabilities <- model$lambda
mu <- model$mu # -0.02026 &  -1.96985
sigma <- model$sigma # 0.4074  & 1.7337
# N <- 1e5
# group <- sample(x = length(probabilities), size = N, replace = T , prob = probabilities)
# x <- rnorm(n = N, mean = mu[group], sd = sigma[group])
# hist(x,freq = F,ylim = c(0,0.5))
# lines(density(x),lty=2,lwd=2)


# Distribution 1
X <- sort(plotMuts$log_FC)
Y <- dnorm(x = X, mean = mu[1], sd = sigma[1])
plot(
  x = X,
  y = Y,
  xlab = "X",
  ylab = "dnorm(x)",
  pch = 20,
  col = "#00000011"
)

Y=2*pnorm(q = -abs(X), mean = mu[1], sd = sigma[1])
plot(x = X,y = Y,xlab="X",ylab="2*pnorm(-abs(x))",pch=20,col="#00000011")
plotMuts$P.value <-
  2 * pnorm(q = -(abs(plotMuts$log_FC)),
            mean = mu[1],
            sd = sigma[1])

plotMuts$FDR_1 <- p.adjust(p = plotMuts$P.value, method = "bonferroni")
plotMuts[,P.value:=NULL]
View(plotMuts[FDR_1<0.05 & log_FC>1,])
sub <- plotMuts[FDR_1<0.05 & log_FC>1,]
write.table(
  sub,
  file = "FDR_0.05_dist1_hits.csv",
  sep = ",",
  quote = F,
  col.names = T,
  row.names = F
)

# Distribution 2
X <- sort(plotMuts$log_FC)
Y <- dnorm(x = X, mean = mu[2], sd = sigma[2])
plot(
  x = X,
  y = Y,
  xlab = "X",
  ylab = "dnorm(x)",
  pch = 20,
  col = "#00000011"
)

Y=2*pnorm(q = -abs(X), mean = mu[2], sd = sigma[2])
plot(x = X,y = Y,xlab="X",ylab="2*pnorm(-abs(x))",pch=20,col="#00000011")
plotMuts$P.value <-
  2 * pnorm(q = -(abs(plotMuts$log_FC)),
            mean = mu[2],
            sd = sigma[2])

plotMuts$FDR_2 <- p.adjust(p = plotMuts$P.value, method = "bonferroni")
plotMuts[,P.value:=NULL]
View(plotMuts[FDR_2<0.00001 & log_FC>1,])
sub <- plotMuts[FDR_2<0.00001 & log_FC>1,]
write.table(
  sub,
  file = "FDR_0.00001_dist2_hits.csv",
  sep = ",",
  quote = F,
  col.names = T,
  row.names = F
)


# model <- spEMsymloc(plotMuts$log_FC,mu=c(0,-2.5),bw=1)
# lines(density(plotMuts$log_FC),lty=2,lwd=2)
# plot(
# model,
# which = 2,
# density = TRUE,
# cex.axis = 1.4,
# cex.lab = 1.4,
# cex.main = 1.8,
# main2 = "Semi-parametric mixture ",
# xlab2 = "log2(fold change)"
# )
