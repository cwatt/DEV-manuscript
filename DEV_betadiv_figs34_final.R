
# Cassandra Wattenburger
# 01/24/18

# Cleaned up scripts to generate beta diversity data and figs 2 and 4

############################################

library("phyloseq")
library("vegan")
library("ggplot2")
library("plyr")
library("ggplot2")
library("tidyr")
library("reshape2")
library("cowplot")

setwd("Y:/Cassi/Sequencing Data/Argonne 2017/01 DEV Project/03 Final Scripts and figures")


#########
#16s Data

rm(list=ls())

###################
#Create ordinations

#load phyloseq object with normalized and pre-processed count data
load("DEV16S.physeq.want.RData")
physeq.want

############################################################################################################################
load("DEV16S.ord.nmds.RData")
############################################################################################################################

#calculate distance using Bray-Curtis method
dist <- distance(physeq.want, method="bray")

#PCoA ordination
ord.pcoa <- ordinate(physeq.want, "PCoA", distance=dist)

#graph PCoA
ord.pcoa.plot <- plot_ordination(physeq.want, ord.pcoa, color = "Compartment", shape = "Soil") + 
  facet_wrap(~Timepoint) +
  ggtitle("16S PCoA using Bray-Curtis Dissimilarity") + 
  stat_ellipse(type="t")
ord.pcoa.plot 

#NMDS ordination
set.seed(2)
ord.nmds <- ordinate(physeq.want, "NMDS", distance=dist)
ord.nmds

#graph NMDS
ord.nmds.plot <- plot_ordination(physeq.want, ord.nmds, color = "Compartment", shape = "Soil") +
  facet_wrap(~Timepoint) +
  ggtitle("16S NMDS using Bray-Curtis Dissimilarity") + 
  stat_ellipse(type="t")
ord.nmds.plot
#stress = 0.142
#PCoA and NMDS give similar ordinations
#I'll stick with NMDS because it's non-euclidean


###########
#Statistics

#make a dataframe to run PErMANOVA on
adonis.df = as(sample_data(physeq.want), "data.frame")

#Full model without block effect, since it is a random effect
#constrained to block
set.seed(2)
full <- adonis(dist~Soil*Compartment*Timepoint, data = adonis.df, strata=adonis.df$Block, permutations = 9999)
full

#                             Df  SumsOfSqs MeanSqs F.Model R2    Pr(>F)    
#Soil                         1    1.6901  1.6901  24.248 0.06457 0.0001 ***
#Compartment                  2    7.2639  3.6320  52.108 0.27753 0.0001 ***
#Timepoint                    3    1.4585  0.4862   6.975 0.05573 0.0001 ***
#Soil:Compartment             2    0.3913  0.1957   2.807 0.01495 0.0011 ** 
#Soil:Timepoint               3    0.2740  0.0913   1.310 0.01047 0.0943 .  
#Compartment:Timepoint        6    1.5458  0.2576   3.696 0.05906 0.0001 ***
#Soil:Compartment:Timepoint   6    0.3767  0.0628   0.901 0.01439 0.5692    
#Residuals                  189   13.1734  0.0697         0.50330           
#Total                      212   26.1738                 1.00000      

#another full model to see if significance changes
set.seed(2)
full2 <- adonis(dist~Compartment*Soil*Timepoint, data = adonis.df, strata=adonis.df$Block, permutations = 9999)
full2
#no change in significance

#Mixed model to investigate block effect
#Block is a random effect so exclude from interactions
set.seed(2)
mixed <- adonis(dist~Soil*Compartment*Timepoint+Block, data = adonis.df, strata=adonis.df$Block, permutations = 9999)
mixed

#Call:
#adonis(formula = dist ~ Soil * Compartment * Timepoint + Block,      data = adonis.df, permutations = 9999, strata = adonis.df$Block) 

#Blocks:  strata 
#Permutation: free
#Number of permutations: 9999

#Terms added sequentially (first to last)

#Df                               SumsOfSqs MeanSqs F.Model R2      Pr(>F)    
#Soil                         1    1.6901  1.6901  26.530 0.06457 0.0001 ***
#Compartment                  2    7.2639  3.6320  57.012 0.27753 0.0001 ***
#Timepoint                    3    1.4585  0.4862   7.632 0.05573 0.0001 ***
#Block                        2    1.2601  0.6300   9.890 0.04814 0.0001 ***
#Soil:Compartment             2    0.3930  0.1965   3.084 0.01501 0.0009 ***
#Soil:Timepoint               3    0.2746  0.0915   1.437 0.01049 0.0714 .  
#Compartment:Timepoint        6    1.5450  0.2575   4.042 0.05903 0.0001 ***
#Soil:Compartment:Timepoint   6    0.3758  0.0626   0.983 0.01436 0.4732    
#Residuals                  187   11.9128  0.0637         0.45514           
#Total                      212   26.1738                 1.00000

#graph to try to see interaction of comp*time
ord.nmds.plot2 <- plot_ordination(physeq.want, ord.nmds, color = "Timepoint") + 
  facet_grid(Soil~Compartment) +
  ggtitle("16S NMDS using Bray-Curtis Dissimilarity") + 
  stat_ellipse(type="t")
ord.nmds.plot2 
#rhizoplane seems to change across time differently than bulk and rhizosphere

##########
#Contrasts
#To pick apart significance between communities of interest
#Will use PERMANOVA with Bonferroni correction for multiple testing

#create data frame

#OTU data
otu = as(otu_table(physeq.want), "matrix")
head(otu[,1:5])
otu.trans <- t(otu)
head(otu.trans[,1:5])
dim(otu.trans)

#meta data
meta <- sample_data(physeq.want)
head(meta)
dim(meta)

#Merge otu and metadata
otu.meta <- merge(meta, otu.trans, by=0)
head(otu.meta[,1:10])
dim(otu.meta)
rownames(otu.meta) <- otu.meta[,1]
head(otu.meta[,1:10])
otu.meta <- otu.meta[,-1]
head(otu.meta[,1:10])


#################################################
#Bulk vs rhizosphere for each time point and soil

#TP1 Conv.
head(otu.meta[,1:10])
tp1.conv.brs <- subset(otu.meta, Timepoint=="TP1")
tp1.conv.brs <- subset(tp1.conv.brs, Soil=="Conv.")
tp1.conv.brs <- tp1.conv.brs[!(tp1.conv.brs$Compartment=="rhizoplane"),]
head(tp1.conv.brs[,1:10])
dim(tp1.conv.brs)

dist.tp1.conv.brs <- vegdist(tp1.conv.brs[,-(1:5)], method="bray")
set.seed(3)
test1 <- adonis(dist.tp1.conv.brs~Compartment, data=tp1.conv.brs, strata=tp1.conv.brs$Block, permutation=9999)
tp1.conv.brs.p <- test1$aov.tab$`Pr(>F)`[1]
tp1.conv.brs.p
#P=0.0058

#TP1 Div.
tp1.div.brs <- subset(otu.meta, Timepoint=="TP1")
tp1.div.brs <- subset(tp1.div.brs, Soil=="Div.")
tp1.div.brs <- tp1.div.brs[!(tp1.div.brs$Compartment=="rhizoplane"),]
head(tp1.div.brs[,1:10])
dim(tp1.div.brs)

dist.tp1.div.brs <- vegdist(tp1.div.brs[,-(1:5)], method="bray")
set.seed(3)
test2 <- adonis(dist.tp1.div.brs~Compartment, data=tp1.div.brs, strata=tp1.div.brs$Block, permutation=9999)
tp1.div.brs.p <- test2$aov.tab$`Pr(>F)`[1]
tp1.div.brs.p
#P=0.0551

#TP2 Conv.
tp2.conv.brs <- subset(otu.meta, Timepoint=="TP2")
tp2.conv.brs <- subset(tp2.conv.brs, Soil=="Conv.")
tp2.conv.brs <- tp2.conv.brs[!(tp2.conv.brs$Compartment=="rhizoplane"),]
head(tp2.conv.brs[,1:10])
dim(tp2.conv.brs)

dist.tp2.conv.brs <- vegdist(tp2.conv.brs[,-(1:5)], method="bray")
set.seed(3)
test3 <- adonis(dist.tp2.conv.brs~Compartment, data=tp2.conv.brs, strata=tp2.conv.brs$Block, permutation=9999)
tp2.conv.brs.p <- test3$aov.tab$`Pr(>F)`[1]
tp2.conv.brs.p
#P=0.0007

#TP2 Div.
tp2.div.brs <- subset(otu.meta, Timepoint=="TP2")
tp2.div.brs <- subset(tp2.div.brs, Soil=="Div.")
tp2.div.brs <- tp2.div.brs[!(tp2.div.brs$Compartment=="rhizoplane"),]
head(tp2.div.brs[,1:10])
dim(tp2.div.brs)

dist.tp2.div.brs <- vegdist(tp2.div.brs[,-(1:5)], method="bray")
set.seed(3)
test4 <- adonis(dist.tp2.div.brs~Compartment, data=tp2.div.brs, strata=tp2.div.brs$Block, permutation=9999)
tp2.div.brs.p <- test4$aov.tab$`Pr(>F)`[1]
tp2.div.brs.p
#P=0.0086

#TP3 Conv.
head(otu.meta[,1:10])
tp3.conv.brs <- subset(otu.meta, Timepoint=="TP3")
tp3.conv.brs <- subset(tp3.conv.brs, Soil=="Conv.")
tp3.conv.brs <- tp3.conv.brs[!(tp3.conv.brs$Compartment=="rhizoplane"),]
head(tp3.conv.brs[,1:10])
dim(tp3.conv.brs)

dist.tp3.conv.brs <- vegdist(tp3.conv.brs[,-(1:5)], method="bray")
set.seed(3)
test5 <- adonis(dist.tp3.conv.brs~Compartment, data=tp3.conv.brs, strata=tp3.conv.brs$Block, permutation=9999)
tp3.conv.brs.p <- test5$aov.tab$`Pr(>F)`[1]
tp3.conv.brs.p
#P=0.0002

#TP3 Div.
tp3.div.brs <- subset(otu.meta, Timepoint=="TP3")
tp3.div.brs <- subset(tp3.div.brs, Soil=="Div.")
tp3.div.brs <- tp3.div.brs[!(tp3.div.brs$Compartment=="rhizoplane"),]
head(tp3.div.brs[,1:10])
dim(tp3.div.brs)

dist.tp3.div.brs <- vegdist(tp3.div.brs[,-(1:5)], method="bray")
set.seed(3)
test6 <- adonis(dist.tp3.div.brs~Compartment, data=tp3.div.brs, strata=tp3.div.brs$Block, permutation=9999)
tp3.div.brs.p <- test6$aov.tab$`Pr(>F)`[1]
tp3.div.brs.p
#P=0.1463

#TP4 Conv.
head(otu.meta[,1:10])
tp4.conv.brs <- subset(otu.meta, Timepoint=="TP4")
tp4.conv.brs <- subset(tp4.conv.brs, Soil=="Conv.")
tp4.conv.brs <- tp4.conv.brs[!(tp4.conv.brs$Compartment=="rhizoplane"),]
head(tp4.conv.brs[,1:10])
dim(tp4.conv.brs)

dist.tp4.conv.brs <- vegdist(tp4.conv.brs[,-(1:5)], method="bray")
set.seed(3)
test7 <- adonis(dist.tp4.conv.brs~Compartment, data=tp4.conv.brs, strata=tp4.conv.brs$Block, permutation=9999)
tp4.conv.brs.p <- test7$aov.tab$`Pr(>F)`[1]
tp4.conv.brs.p
#P=0.0037

#TP4 Div.
tp4.div.brs <- subset(otu.meta, Timepoint=="TP4")
tp4.div.brs <- subset(tp4.div.brs, Soil=="Div.")
tp4.div.brs <- tp4.div.brs[!(tp4.div.brs$Compartment=="rhizoplane"),]
head(tp4.div.brs[,1:10])
dim(tp4.div.brs)

dist.tp4.div.brs <- vegdist(tp4.div.brs[,-(1:5)], method="bray")
set.seed(3)
test8 <- adonis(dist.tp4.div.brs~Compartment, data=tp4.div.brs, strata=tp4.div.brs$Block, permutation=9999)
tp4.div.brs.p <- test8$aov.tab$`Pr(>F)`[1]
tp4.div.brs.p
#P=0.0504

#Adjust for multiple comparisons
#54 is the total number of contrasts done in this script
p.adj1 <- p.adjust(c(tp1.conv.brs.p, tp1.div.brs.p, tp2.conv.brs.p, tp2.div.brs.p,
                     tp3.conv.brs.p, tp3.div.brs.p, tp4.conv.brs.p, tp4.div.brs.p), method="bonferroni", n=54)
p.adj1


###############################################
#Bulk vs rhizoplane for each timepoint and soil

#TP1 Conv.
head(otu.meta[,1:10])
tp1.conv.brp <- subset(otu.meta, Timepoint=="TP1")
tp1.conv.brp <- subset(tp1.conv.brp, Soil=="Conv.")
tp1.conv.brp <- tp1.conv.brp[!(tp1.conv.brp$Compartment=="rhizosphere"),]
head(tp1.conv.brp[,1:10])
dim(tp1.conv.brp)

dist.tp1.conv.brp <- vegdist(tp1.conv.brp[,-(1:5)], method="bray")
set.seed(3)
test9 <- adonis(dist.tp1.conv.brp~Compartment, data=tp1.conv.brp, strata=tp1.conv.brp$Block, permutation=9999)
tp1.conv.brp.p <- test9$aov.tab$`Pr(>F)`[1]
tp1.conv.brp.p
#P=0.0004

#TP1 Div.
tp1.div.brp <- subset(otu.meta, Timepoint=="TP1")
tp1.div.brp <- subset(tp1.div.brp, Soil=="Div.")
tp1.div.brp <- tp1.div.brp[!(tp1.div.brp$Compartment=="rhizosphere"),]
head(tp1.div.brp[,1:10])
dim(tp1.div.brp)

dist.tp1.div.brp <- vegdist(tp1.div.brp[,-(1:5)], method="bray")
set.seed(3)
test10 <- adonis(dist.tp1.div.brp~Compartment, data=tp1.div.brp, strata=tp1.div.brp$Block, permutation=9999)
tp1.div.brp.p <- test10$aov.tab$`Pr(>F)`[1]
tp1.div.brp.p
#P=0.0001

#TP2 Conv.
tp2.conv.brp <- subset(otu.meta, Timepoint=="TP2")
tp2.conv.brp <- subset(tp2.conv.brp, Soil=="Conv.")
tp2.conv.brp <- tp2.conv.brp[!(tp2.conv.brp$Compartment=="rhizosphere"),]
head(tp2.conv.brp[,1:10])
dim(tp2.conv.brp)

dist.tp2.conv.brp <- vegdist(tp2.conv.brp[,-(1:5)], method="bray")
set.seed(3)
test11 <- adonis(dist.tp2.conv.brp~Compartment, data=tp2.conv.brp, strata=tp2.conv.brp$Block, permutation=9999)
tp2.conv.brp.p <- test11$aov.tab$`Pr(>F)`[1]
tp2.conv.brp.p
#P=0.0005, significant

#TP2 Div.
tp2.div.brp <- subset(otu.meta, Timepoint=="TP2")
tp2.div.brp <- subset(tp2.div.brp, Soil=="Div.")
tp2.div.brp <- tp2.div.brp[!(tp2.div.brp$Compartment=="rhizosphere"),]
head(tp2.div.brp[,1:10])
dim(tp2.div.brp)

dist.tp2.div.brp <- vegdist(tp2.div.brp[,-(1:5)], method="bray")
set.seed(3)
test12 <- adonis(dist.tp2.div.brp~Compartment, data=tp2.div.brp, strata=tp2.div.brp$Block, permutation=9999)
tp2.div.brp.p <- test12$aov.tab$`Pr(>F)`[1]
tp2.div.brp.p
#P=0.0005

#TP3 Conv.
head(otu.meta[,1:10])
tp3.conv.brp <- subset(otu.meta, Timepoint=="TP3")
tp3.conv.brp <- subset(tp3.conv.brp, Soil=="Conv.")
tp3.conv.brp <- tp3.conv.brp[!(tp3.conv.brp$Compartment=="rhizosphere"),]
head(tp3.conv.brp[,1:10])
dim(tp3.conv.brp)

dist.tp3.conv.brp <- vegdist(tp3.conv.brp[,-(1:5)], method="bray")
set.seed(3)
test13 <- adonis(dist.tp3.conv.brp~Compartment, data=tp3.conv.brp, strata=tp3.conv.brp$Block, permutation=9999)
tp3.conv.brp.p <- test13$aov.tab$`Pr(>F)`[1]
tp3.conv.brp.p
#P=0.0004

#TP3 Div.
tp3.div.brp <- subset(otu.meta, Timepoint=="TP3")
tp3.div.brp <- subset(tp3.div.brp, Soil=="Div.")
tp3.div.brp <- tp3.div.brp[!(tp3.div.brp$Compartment=="rhizosphere"),]
head(tp3.div.brp[,1:10])
dim(tp3.div.brp)

dist.tp3.div.brp <- vegdist(tp3.div.brp[,-(1:5)], method="bray")
set.seed(3)
test14 <- adonis(dist.tp3.div.brp~Compartment, data=tp3.div.brp, strata=tp3.div.brp$Block, permutation=9999)
tp3.div.brp.p <- test14$aov.tab$`Pr(>F)`[1]
tp3.div.brp.p
#P=0.0005

#TP4 Conv.
head(otu.meta[,1:10])
tp4.conv.brp <- subset(otu.meta, Timepoint=="TP4")
tp4.conv.brp <- subset(tp4.conv.brp, Soil=="Conv.")
tp4.conv.brp <- tp4.conv.brp[!(tp4.conv.brp$Compartment=="rhizosphere"),]
head(tp4.conv.brp[,1:10])
dim(tp4.conv.brp)

dist.tp4.conv.brp <- vegdist(tp4.conv.brp[,-(1:5)], method="bray")
set.seed(3)
test15 <- adonis(dist.tp4.conv.brp~Compartment, data=tp4.conv.brp, strata=tp4.conv.brp$Block, permutation=9999)
tp4.conv.brp.p <- test15$aov.tab$`Pr(>F)`[1]
tp4.conv.brp.p
#P=0.0005

#TP4 Div.
tp4.div.brp <- subset(otu.meta, Timepoint=="TP4")
tp4.div.brp <- subset(tp4.div.brp, Soil=="Div.")
tp4.div.brp <- tp4.div.brp[!(tp4.div.brp$Compartment=="rhizosphere"),]
head(tp4.div.brp[,1:10])
dim(tp4.div.brp)

dist.tp4.div.brp <- vegdist(tp4.div.brp[,-(1:5)], method="bray")
set.seed(3)
test16 <- adonis(dist.tp4.div.brp~Compartment, data=tp4.div.brp, strata=tp4.div.brp$Block, permutation=9999)
tp4.div.brp.p <- test16$aov.tab$`Pr(>F)`[1]
tp4.div.brp.p
#P=0.001

#Adjust p-values for multiple comparisons
p.adj2 <- p.adjust(c(tp1.conv.brp.p, tp1.div.brp.p, tp2.conv.brp.p, tp2.div.brp.p,
                     tp3.conv.brp.p, tp3.div.brp.p, tp4.conv.brp.p, tp4.div.brp.p), method="bonferroni", n=54)
p.adj2
#all signficant


###############################
#Bulk vs bulk for each timepoint

#TP1
tp1.b <- subset(otu.meta, Timepoint=="TP1")
tp1.b <- subset(tp1.b, Compartment=="bulk")
head(tp1.b[,1:5])
dim(tp1.b)

dist.tp1.b <- vegdist(tp1.b[,-(1:5)], method="bray")
set.seed(3)
test17 <- adonis(dist.tp1.b~Soil, data=tp1.b, strata=tp1.b$Block, permutation=9999)
tp1.b.p <- test17$aov.tab$`Pr(>F)`[1]
tp1.b.p
#P=0.0005

#TP2
tp2.b <- subset(otu.meta, Timepoint=="TP2")
tp2.b <- subset(tp2.b, Compartment=="bulk")
head(tp2.b[,1:5])
dim(tp2.b)

dist.tp2.b <- vegdist(tp2.b[,-(1:5)], method="bray")
set.seed(3)
test18 <- adonis(dist.tp2.b~Soil, data=tp2.b, strata=tp2.b$Block, permutation=9999)
tp2.b.p <- test18$aov.tab$`Pr(>F)`[1]
tp2.b.p
#P=0.0005

#TP3
tp3.b <- subset(otu.meta, Timepoint=="TP3")
tp3.b <- subset(tp3.b, Compartment=="bulk")
head(tp3.b[,1:5])
dim(tp3.b)

dist.tp3.b <- vegdist(tp3.b[,-(1:5)], method="bray")
set.seed(3)
test19 <- adonis(dist.tp3.b~Soil, data=tp3.b, strata=tp3.b$Block, permutation=9999)
tp3.b.p <- test19$aov.tab$`Pr(>F)`[1]
tp3.b.p
#P=0.0009

#TP4
tp4.b <- subset(otu.meta, Timepoint=="TP4")
tp4.b <- subset(tp4.b, Compartment=="bulk")
head(tp4.b[,1:5])
dim(tp4.b)

dist.tp4.b <- vegdist(tp4.b[,-(1:5)], method="bray")
set.seed(3)
test20 <- adonis(dist.tp4.b~Soil, data=tp4.b, strata=tp4.b$Block, permutation=9999)
tp4.b.p <- test20$aov.tab$`Pr(>F)`[1]
tp4.b.p
#P=0.0005

#Mulitple comparisons p adjustment
p.adj3 <- p.adjust(c(tp1.b.p, tp2.b.p, tp3.b.p, tp4.b.p), method="bonferroni", n=54)
p.adj3
#all significant


############################################
#Rhizosphere vs rhizosphere at each timepoint

#TP1
tp1.rs <- subset(otu.meta, Timepoint=="TP1")
tp1.rs <- subset(tp1.rs, Compartment=="rhizosphere")
head(tp1.rs[,1:5])
dim(tp1.rs)

dist.tp1.rs <- vegdist(tp1.rs[,-(1:5)], method="bray")
set.seed(3)
test21 <- adonis(dist.tp1.rs~Soil, data=tp1.rs, strata=tp1.rs$Block, permutation=9999)
tp1.rs.p <- test21$aov.tab$`Pr(>F)`[1]
tp1.rs.p
#P=0.0005

#TP2
tp2.rs <- subset(otu.meta, Timepoint=="TP2")
tp2.rs <- subset(tp2.rs, Compartment=="rhizosphere")
head(tp2.rs[,1:5])
dim(tp2.rs)

dist.tp2.rs <- vegdist(tp2.rs[,-(1:5)], method="bray")
set.seed(3)
test22 <- adonis(dist.tp2.rs~Soil, data=tp2.rs, strata=tp2.rs$Block, permutation=9999)
tp2.rs.p <- test22$aov.tab$`Pr(>F)`[1]
tp2.rs.p
#P=0.0005

#TP3
tp3.rs <- subset(otu.meta, Timepoint=="TP3")
tp3.rs <- subset(tp3.rs, Compartment=="rhizosphere")
head(tp3.rs[,1:5])
dim(tp3.rs)

dist.tp3.rs <- vegdist(tp3.rs[,-(1:5)], method="bray")
set.seed(3)
test23 <- adonis(dist.tp3.rs~Soil, data=tp3.rs, strata=tp3.rs$Block, permutation=9999)
tp3.rs.p <- test23$aov.tab$`Pr(>F)`[1]
tp3.rs.p
#P=0.0009

#TP4
tp4.rs <- subset(otu.meta, Timepoint=="TP4")
tp4.rs <- subset(tp4.rs, Compartment=="rhizosphere")
head(tp4.rs[,1:5])
dim(tp4.rs)

set.seed(3)
dist.tp4.rs <- vegdist(tp4.rs[,-(1:5)], method="bray")
test24 <- adonis(dist.tp4.rs~Soil, data=tp4.rs, strata=tp4.rs$Block, permutation=9999)
tp4.rs.p <- test24$aov.tab$`Pr(>F)`[1]
tp4.rs.p
#P=0.0005, significant

#Multiple comparisons p-adjustment
p.adj4 <- p.adjust(c(tp1.rs.p, tp2.rs.p, tp3.rs.p, tp4.rs.p), method="bonferroni", n=54)
p.adj4
#all significant

############################################
#Rhizoplane vs rhizoplane at each timepoint

#TP1
tp1.rp <- subset(otu.meta, Timepoint=="TP1")
tp1.rp <- subset(tp1.rp, Compartment=="rhizoplane")
head(tp1.rp[,1:5])
dim(tp1.rp)

dist.tp1.rp <- vegdist(tp1.rp[,-(1:5)], method="bray")
set.seed(3)
test25 <- adonis(dist.tp1.rp~Soil, data=tp1.rp, strata=tp1.rp$Block, permutation=9999)
tp1.rp.p <- test25$aov.tab$`Pr(>F)`[1]
tp1.rp.p
#P=0.0001

#TP2
tp2.rp <- subset(otu.meta, Timepoint=="TP2")
tp2.rp <- subset(tp2.rp, Compartment=="rhizoplane")
head(tp2.rp[,1:5])
dim(tp2.rp)

dist.tp2.rp <- vegdist(tp2.rp[,-(1:5)], method="bray")
set.seed(3)
test26 <- adonis(dist.tp2.rp~Soil, data=tp2.rp, strata=tp2.rp$Block, permutation=9999)
tp2.rp.p <- test26$aov.tab$`Pr(>F)`[1]
tp2.rp.p
#P=0.0005

#TP3
tp3.rp <- subset(otu.meta, Timepoint=="TP3")
tp3.rp <- subset(tp3.rp, Compartment=="rhizoplane")
head(tp3.rp[,1:5])
dim(tp3.rp)

dist.tp3.rp <- vegdist(tp3.rp[,-(1:5)], method="bray")
set.seed(3)
test27 <- adonis(dist.tp3.rp~Soil, data=tp3.rp, strata=tp3.rp$Block, permutation=9999)
tp3.rp.p <- test27$aov.tab$'Pr(>F)'[1]
tp3.rp.p
#P=0.0011

#TP4
tp4.rp <- subset(otu.meta, Timepoint=="TP4")
tp4.rp <- subset(tp4.rp, Compartment=="rhizoplane")
head(tp4.rp[,1:5])
dim(tp4.rp)

dist.tp4.rp <- vegdist(tp4.rp[,-(1:5)], method="bray")
set.seed(3)
test28 <- adonis(dist.tp4.rp~Soil, data=tp4.rp, strata=tp4.rp$Block, permutation=9999)
tp4.rp.p <- test28$aov.tab$'Pr(>F)'[1]
tp4.rp.p
#P=0.0030

#Multiple comparisons p-adjustment
p.adj5 <- p.adjust(c(tp1.rp.p, tp2.rp.p, tp3.rp.p, tp4.rp.p), method="bonferroni", n=54)
p.adj5
#TP1 and TP2 significant, TP3 and TP4 not


#Explore CompxTime interaction because it's the most biologically significant interaction with visual support in ordinations

##########
#Bulk soil

#TP1 vs TP2 Conv.
conv.b.tp12 <- subset(otu.meta, Compartment=="bulk")
conv.b.tp12 <- subset(conv.b.tp12, Soil=="Conv.")
conv.b.tp12 <- conv.b.tp12[!(conv.b.tp12$Timepoint=="TP3"),]
conv.b.tp12 <- conv.b.tp12[!(conv.b.tp12$Timepoint=="TP4"),]  
head(conv.b.tp12[,1:10])

dist.conv.b.tp12 <- vegdist(conv.b.tp12[,-(1:5)], method="bray")
set.seed(3)
test29 <- adonis(dist.conv.b.tp12~Timepoint, data=conv.b.tp12, strata=conv.b.tp12$Block, permutation=9999)
conv.b.tp12.p <- test29$aov.tab$'Pr(>F)'[1]
conv.b.tp12.p
#P=0.003

#TP1 vs TP2 Div.
div.b.tp12 <- subset(otu.meta, Compartment=="bulk")
div.b.tp12 <- subset(div.b.tp12, Soil=="Div.")
div.b.tp12 <- div.b.tp12[!(div.b.tp12$Timepoint=="TP3"),]
div.b.tp12 <- div.b.tp12[!(div.b.tp12$Timepoint=="TP4"),]  
head(div.b.tp12[,1:10])

dist.div.b.tp12 <- vegdist(div.b.tp12[,-(1:5)], method="bray")
set.seed(3)
test30 <- adonis(dist.div.b.tp12~Timepoint, data=div.b.tp12, strata=div.b.tp12$Block, permutation=9999)
div.b.tp12.p <- test30$aov.tab$'Pr(>F)'[1]
div.b.tp12.p
#0.1492

#TP2 vs TP3 Conv.
conv.b.tp23 <- subset(otu.meta, Compartment=="bulk")
conv.b.tp23 <- subset(conv.b.tp23, Soil=="Conv.")
conv.b.tp23 <- conv.b.tp23[!(conv.b.tp23$Timepoint=="TP1"),]
conv.b.tp23 <- conv.b.tp23[!(conv.b.tp23$Timepoint=="TP4"),]  
head(conv.b.tp23[,1:10])

dist.conv.b.tp23 <- vegdist(conv.b.tp23[,-(1:5)], method="bray")
set.seed(3)
test31 <- adonis(dist.conv.b.tp23~Timepoint, data=conv.b.tp23, strata=conv.b.tp23$Block, permutation=9999)
conv.b.tp23.p <- test31$aov.tab$'Pr(>F)'[1]
conv.b.tp23.p
#P=0.0153

#TP2 vs TP3 Div.
div.b.tp23 <- subset(otu.meta, Compartment=="bulk")
div.b.tp23 <- subset(div.b.tp23, Soil=="Div.")
div.b.tp23 <- div.b.tp23[!(div.b.tp23$Timepoint=="TP1"),]
div.b.tp23 <- div.b.tp23[!(div.b.tp23$Timepoint=="TP4"),]  
head(div.b.tp23[,1:10])

dist.div.b.tp23 <- vegdist(div.b.tp23[,-(1:5)], method="bray")
set.seed(3)
test32 <- adonis(dist.div.b.tp23~Timepoint, data=div.b.tp23, strata=div.b.tp23$Block, permutation=9999)
div.b.tp23.p <- test32$aov.tab$'Pr(>F)'[1]
div.b.tp23.p
#P=0.6171

#TP3 vs TP4 Conv.
conv.b.tp34 <- subset(otu.meta, Compartment=="bulk")
conv.b.tp34 <- subset(conv.b.tp34, Soil=="Conv.")
conv.b.tp34 <- conv.b.tp34[!(conv.b.tp34$Timepoint=="TP1"),]
conv.b.tp34 <- conv.b.tp34[!(conv.b.tp34$Timepoint=="TP2"),]  
head(conv.b.tp34[,1:10])

dist.conv.b.tp34 <- vegdist(conv.b.tp34[,-(1:5)], method="bray")
set.seed(3)
test33 <- adonis(dist.conv.b.tp34~Timepoint, data=conv.b.tp34, strata=conv.b.tp34$Block, permutation=9999)
conv.b.tp34.p <- test33$aov.tab$'Pr(>F)'[1]
conv.b.tp34.p
#0.1566

#TP3 vs TP4 Div.
div.b.tp34 <- subset(otu.meta, Compartment=="bulk")
div.b.tp34 <- subset(div.b.tp34, Soil=="Div.")
div.b.tp34 <- div.b.tp34[!(div.b.tp34$Timepoint=="TP1"),]
div.b.tp34 <- div.b.tp34[!(div.b.tp34$Timepoint=="TP2"),]  
head(div.b.tp34[,1:10])

dist.div.b.tp34 <- vegdist(div.b.tp34[,-(1:5)], method="bray")
set.seed(3)
test34 <- adonis(dist.div.b.tp34~Timepoint, data=div.b.tp34, strata=div.b.tp34$Block, permutation=9999)
div.b.tp34.p <- test34$aov.tab$'Pr(>F)'[1]
div.b.tp34.p
#0.6474

#Multiple comparisons p adjustment
padj6 <- p.adjust(c(conv.b.tp12.p, div.b.tp12.p, conv.b.tp23.p, div.b.tp23.p, conv.b.tp34.p, div.b.tp34.p),
                  method="bonferroni", n=54)
padj6
#no significance


#################
#Rhizosphere soil

#TP1 vs TP2 Conv.
conv.rs.tp12 <- subset(otu.meta, Compartment=="rhizosphere")
conv.rs.tp12 <- subset(conv.rs.tp12, Soil=="Conv.")
conv.rs.tp12 <- conv.rs.tp12[!(conv.rs.tp12$Timepoint=="TP3"),]
conv.rs.tp12 <- conv.rs.tp12[!(conv.rs.tp12$Timepoint=="TP4"),]  
head(conv.rs.tp12[,1:10])

dist.conv.rs.tp12 <- vegdist(conv.rs.tp12[,-(1:5)], method="bray")
set.seed(3)
test35 <- adonis(dist.conv.rs.tp12~Timepoint, data=conv.rs.tp12, strata=conv.rs.tp12$Block, permutation=9999)
conv.rs.tp12.p <- test35$aov.tab$'Pr(>F)'[1]
conv.rs.tp12.p
#0.0021

#TP1 vs TP2 Div.
div.rs.tp12 <- subset(otu.meta, Compartment=="rhizosphere")
div.rs.tp12 <- subset(div.rs.tp12, Soil=="Div.")
div.rs.tp12 <- div.rs.tp12[!(div.rs.tp12$Timepoint=="TP3"),]
div.rs.tp12 <- div.rs.tp12[!(div.rs.tp12$Timepoint=="TP4"),]  
head(div.rs.tp12[,1:10])

dist.div.rs.tp12 <- vegdist(div.rs.tp12[,-(1:5)], method="bray")
set.seed(3)
test36 <- adonis(dist.div.rs.tp12~Timepoint, data=div.rs.tp12, strata=div.rs.tp12$Block, permutation=9999)
div.rs.tp12.p <- test36$aov.tab$'Pr(>F)'[1]
div.rs.tp12.p
#0.0021

#TP2 vs TP3 Conv.
conv.rs.tp23 <- subset(otu.meta, Compartment=="rhizosphere")
conv.rs.tp23 <- subset(conv.rs.tp23, Soil=="Conv.")
conv.rs.tp23 <- conv.rs.tp23[!(conv.rs.tp23$Timepoint=="TP1"),]
conv.rs.tp23 <- conv.rs.tp23[!(conv.rs.tp23$Timepoint=="TP4"),]  
head(conv.rs.tp23[,1:10])

dist.conv.rs.tp23 <- vegdist(conv.rs.tp23[,-(1:5)], method="bray")
set.seed(3)
test37 <- adonis(dist.conv.rs.tp23~Timepoint, data=conv.rs.tp23, strata=conv.rs.tp23$Block, permutation=9999)
conv.rs.tp23.p <- test37$aov.tab$'Pr(>F)'[1]
conv.rs.tp23.p
#0.0092

#TP2 vs TP3 Div.
div.rs.tp23 <- subset(otu.meta, Compartment=="rhizosphere")
div.rs.tp23 <- subset(div.rs.tp23, Soil=="Div.")
div.rs.tp23 <- div.rs.tp23[!(div.rs.tp23$Timepoint=="TP1"),]
div.rs.tp23 <- div.rs.tp23[!(div.rs.tp23$Timepoint=="TP4"),]  
head(div.rs.tp23[,1:10])

dist.div.rs.tp23 <- vegdist(div.rs.tp23[,-(1:5)], method="bray")
set.seed(3)
test38 <- adonis(dist.div.rs.tp23~Timepoint, data=div.rs.tp23, strata=div.rs.tp23$Block, permutation=9999)
div.rs.tp23.p <- test38$aov.tab$'Pr(>F)'[1]
div.rs.tp23.p
#0.0177

#TP3 vs TP4 Conv.
conv.rs.tp34 <- subset(otu.meta, Compartment=="rhizosphere")
conv.rs.tp34 <- subset(conv.rs.tp34, Soil=="Conv.")
conv.rs.tp34 <- conv.rs.tp34[!(conv.rs.tp34$Timepoint=="TP1"),]
conv.rs.tp34 <- conv.rs.tp34[!(conv.rs.tp34$Timepoint=="TP2"),]  
head(conv.rs.tp34[,1:10])

dist.conv.rs.tp34 <- vegdist(conv.rs.tp34[,-(1:5)], method="bray")
set.seed(3)
test39 <- adonis(dist.conv.rs.tp34~Timepoint, data=conv.rs.tp34, strata=conv.rs.tp34$Block, permutation=9999)
conv.rs.tp34.p <- test39$aov.tab$'Pr(>F)'[1]
conv.rs.tp34.p
#0.0357

#TP3 vs TP4 Div.
div.rs.tp34 <- subset(otu.meta, Compartment=="rhizosphere")
div.rs.tp34 <- subset(div.rs.tp34, Soil=="Div.")
div.rs.tp34 <- div.rs.tp34[!(div.rs.tp34$Timepoint=="TP1"),]
div.rs.tp34 <- div.rs.tp34[!(div.rs.tp34$Timepoint=="TP2"),]  
head(div.rs.tp34[,1:10])

dist.div.rs.tp34 <- vegdist(div.rs.tp34[,-(1:5)], method="bray")
set.seed(3)
test40 <- adonis(dist.div.rs.tp34~Timepoint, data=div.rs.tp34, strata=div.rs.tp34$Block, permutation=9999)
div.rs.tp34.p <- test40$aov.tab$'Pr(>F)'[1]
div.rs.tp34.p
#0.055

#Multiple comparisons p adjustment
p.adj7 <- p.adjust(c(conv.rs.tp12.p, div.rs.tp12.p, conv.rs.tp23.p, div.rs.tp23.p, conv.rs.tp34.p, div.rs.tp34.p),
                   method="bonferroni", n=54)
p.adj7
#no significance


################
#Rhizoplane soil

#TP1 vs TP2 Conv.
conv.rp.tp12 <- subset(otu.meta, Compartment=="rhizoplane")
conv.rp.tp12 <- subset(conv.rp.tp12, Soil=="Conv.")
conv.rp.tp12 <- conv.rp.tp12[!(conv.rp.tp12$Timepoint=="TP3"),]
conv.rp.tp12 <- conv.rp.tp12[!(conv.rp.tp12$Timepoint=="TP4"),]  
head(conv.rp.tp12[,1:10])

dist.conv.rp.tp12 <- vegdist(conv.rp.tp12[,-(1:5)], method="bray")
set.seed(3)
test41 <- adonis(dist.conv.rp.tp12~Timepoint, data=conv.rp.tp12, strata=conv.rp.tp12$Block, permutation=9999)
conv.rp.tp12.p <- test41$aov.tab$'Pr(>F)'[1]
conv.rp.tp12.p
#P=0.0005

#TP1 vs TP2 Div.
div.rp.tp12 <- subset(otu.meta, Compartment=="rhizoplane")
div.rp.tp12 <- subset(div.rp.tp12, Soil=="Div.")
div.rp.tp12 <- div.rp.tp12[!(div.rp.tp12$Timepoint=="TP3"),]
div.rp.tp12 <- div.rp.tp12[!(div.rp.tp12$Timepoint=="TP4"),]  
head(div.rp.tp12[,1:10])

dist.div.rp.tp12 <- vegdist(div.rp.tp12[,-(1:5)], method="bray")
set.seed(3)
test42 <- adonis(dist.div.rp.tp12~Timepoint, data=div.rp.tp12, strata=div.rp.tp12$Block, permutation=9999)
div.rp.tp12.p <- test42$aov.tab$'Pr(>F)'[1]
div.rp.tp12.p
#P=0.0004

#TP2 vs TP3 Conv.
conv.rp.tp23 <- subset(otu.meta, Compartment=="rhizoplane")
conv.rp.tp23 <- subset(conv.rp.tp23, Soil=="Conv.")
conv.rp.tp23 <- conv.rp.tp23[!(conv.rp.tp23$Timepoint=="TP1"),]
conv.rp.tp23 <- conv.rp.tp23[!(conv.rp.tp23$Timepoint=="TP4"),]  
head(conv.rp.tp23[,1:10])

dist.conv.rp.tp23 <- vegdist(conv.rp.tp23[,-(1:5)], method="bray")
set.seed(3)
test43 <- adonis(dist.conv.rp.tp23~Timepoint, data=conv.rp.tp23, strata=conv.rp.tp23$Block, permutation=9999)
conv.rp.tp23.p <- test43$aov.tab$'Pr(>F)'[1]
conv.rp.tp23.p
#P=0.0004

#TP2 vs TP3 Div.
div.rp.tp23 <- subset(otu.meta, Compartment=="rhizoplane")
div.rp.tp23 <- subset(div.rp.tp23, Soil=="Div.")
div.rp.tp23 <- div.rp.tp23[!(div.rp.tp23$Timepoint=="TP1"),]
div.rp.tp23 <- div.rp.tp23[!(div.rp.tp23$Timepoint=="TP4"),]  
head(div.rp.tp23[,1:10])

dist.div.rp.tp23 <- vegdist(div.rp.tp23[,-(1:5)], method="bray")
set.seed(3)
test44 <- adonis(dist.div.rp.tp23~Timepoint, data=div.rp.tp23, strata=div.rp.tp23$Block, permutation=9999)
div.rp.tp23.p <- test44$aov.tab$'Pr(>F)'[1]
div.rp.tp23.p
#P=0.0005

#TP3 vs TP4 Conv.
conv.rp.tp34 <- subset(otu.meta, Compartment=="rhizoplane")
conv.rp.tp34 <- subset(conv.rp.tp34, Soil=="Conv.")
conv.rp.tp34 <- conv.rp.tp34[!(conv.rp.tp34$Timepoint=="TP1"),]
conv.rp.tp34 <- conv.rp.tp34[!(conv.rp.tp34$Timepoint=="TP2"),]  
head(conv.rp.tp34[,1:10])

dist.conv.rp.tp34 <- vegdist(conv.rp.tp34[,-(1:5)], method="bray")
set.seed(3)
test45 <- adonis(dist.conv.rp.tp34~Timepoint, data=conv.rp.tp34, strata=conv.rp.tp34$Block, permutation=9999)
conv.rp.tp34.p <- test45$aov.tab$'Pr(>F)'[1]
conv.rp.tp34.p
#P=0.0003

#TP3 vs TP4 Div.
div.rp.tp34 <- subset(otu.meta, Compartment=="rhizoplane")
div.rp.tp34 <- subset(div.rp.tp34, Soil=="Div.")
div.rp.tp34 <- div.rp.tp34[!(div.rp.tp34$Timepoint=="TP1"),]
div.rp.tp34 <- div.rp.tp34[!(div.rp.tp34$Timepoint=="TP2"),]  
head(div.rp.tp34[,1:10])

dist.div.rp.tp34 <- vegdist(div.rp.tp34[,-(1:5)], method="bray")
set.seed(3)
test46 <- adonis(dist.div.rp.tp34~Timepoint, data=div.rp.tp34, strata=div.rp.tp34$Block, permutation=9999)
div.rp.tp34.p <- test46$aov.tab$'Pr(>F)'[1]
div.rp.tp34.p
#P=0.0005

#Correct for multiple comparisons
padj8 <- p.adjust(c(conv.rp.tp12.p, div.rp.tp12.p, conv.rp.tp23.p, div.rp.tp23.p, conv.rp.tp34.p, div.rp.tp34.p),
                  method="bonferroni", n=54)
padj8
#all significant



###################################################
#Rhizosphere vs Rhizoplane for each cropping system

#adding this in after all other comparisons, to make analysis more complete

#TP1 Conv.
head(otu.meta[,1:10])
tp1.conv.rsrp <- subset(otu.meta, Timepoint=="TP1")
tp1.conv.rsrp <- subset(tp1.conv.rsrp, Soil=="Conv.")
tp1.conv.rsrp <- tp1.conv.rsrp[!(tp1.conv.rsrp$Compartment=="bulk"),]
head(tp1.conv.rsrp[,1:10])
dim(tp1.conv.rsrp)

dist.tp1.conv.rsrp <- vegdist(tp1.conv.rsrp[,-(1:5)], method="bray")
set.seed(3)
test47 <- adonis(dist.tp1.conv.rsrp~Compartment, data=tp1.conv.rsrp, strata=tp1.conv.rsrp$Block, permutation=9999)
tp1.conv.rsrp.p <- test47$aov.tab$`Pr(>F)`[1]
tp1.conv.rsrp.p
#P=0.0004

#TP1 Div.
tp1.div.rsrp <- subset(otu.meta, Timepoint=="TP1")
tp1.div.rsrp <- subset(tp1.div.rsrp, Soil=="Div.")
tp1.div.rsrp <- tp1.div.rsrp[!(tp1.div.rsrp$Compartment=="bulk"),]
head(tp1.div.rsrp[,1:10])
dim(tp1.div.rsrp)

dist.tp1.div.rsrp <- vegdist(tp1.div.rsrp[,-(1:5)], method="bray")
set.seed(3)
test48 <- adonis(dist.tp1.div.rsrp~Compartment, data=tp1.div.rsrp, strata=tp1.div.rsrp$Block, permutation=9999)
tp1.div.rsrp.p <- test48$aov.tab$`Pr(>F)`[1]
tp1.div.rsrp.p
#P=0.0001

#TP2 Conv.
tp2.conv.rsrp <- subset(otu.meta, Timepoint=="TP2")
tp2.conv.rsrp <- subset(tp2.conv.rsrp, Soil=="Conv.")
tp2.conv.rsrp <- tp2.conv.rsrp[!(tp2.conv.rsrp$Compartment=="bulk"),]
head(tp2.conv.rsrp[,1:10])
dim(tp2.conv.rsrp)

dist.tp2.conv.rsrp <- vegdist(tp2.conv.rsrp[,-(1:5)], method="bray")
set.seed(3)
test49 <- adonis(dist.tp2.conv.rsrp~Compartment, data=tp2.conv.rsrp, strata=tp2.conv.rsrp$Block, permutation=9999)
tp2.conv.rsrp.p <- test49$aov.tab$`Pr(>F)`[1]
tp2.conv.rsrp.p
#P=0.0005, significant

#TP2 Div.
tp2.div.rsrp <- subset(otu.meta, Timepoint=="TP2")
tp2.div.rsrp <- subset(tp2.div.rsrp, Soil=="Div.")
tp2.div.rsrp <- tp2.div.rsrp[!(tp2.div.rsrp$Compartment=="bulk"),]
head(tp2.div.rsrp[,1:10])
dim(tp2.div.rsrp)

dist.tp2.div.rsrp <- vegdist(tp2.div.rsrp[,-(1:5)], method="bray")
set.seed(3)
test50 <- adonis(dist.tp2.div.rsrp~Compartment, data=tp2.div.rsrp, strata=tp2.div.rsrp$Block, permutation=9999)
tp2.div.rsrp.p <- test50$aov.tab$`Pr(>F)`[1]
tp2.div.rsrp.p
#P=0.0005

#TP3 Conv.
head(otu.meta[,1:10])
tp3.conv.rsrp <- subset(otu.meta, Timepoint=="TP3")
tp3.conv.rsrp <- subset(tp3.conv.rsrp, Soil=="Conv.")
tp3.conv.rsrp <- tp3.conv.rsrp[!(tp3.conv.rsrp$Compartment=="bulk"),]
head(tp3.conv.rsrp[,1:10])
dim(tp3.conv.rsrp)

dist.tp3.conv.rsrp <- vegdist(tp3.conv.rsrp[,-(1:5)], method="bray")
set.seed(3)
test51 <- adonis(dist.tp3.conv.rsrp~Compartment, data=tp3.conv.rsrp, strata=tp3.conv.rsrp$Block, permutation=9999)
tp3.conv.rsrp.p <- test51$aov.tab$`Pr(>F)`[1]
tp3.conv.rsrp.p
#P=0.0004

#TP3 Div.
tp3.div.rsrp <- subset(otu.meta, Timepoint=="TP3")
tp3.div.rsrp <- subset(tp3.div.rsrp, Soil=="Div.")
tp3.div.rsrp <- tp3.div.rsrp[!(tp3.div.rsrp$Compartment=="bulk"),]
head(tp3.div.rsrp[,1:10])
dim(tp3.div.rsrp)

dist.tp3.div.rsrp <- vegdist(tp3.div.rsrp[,-(1:5)], method="bray")
set.seed(3)
test52 <- adonis(dist.tp3.div.rsrp~Compartment, data=tp3.div.rsrp, strata=tp3.div.rsrp$Block, permutation=9999)
tp3.div.rsrp.p <- test52$aov.tab$`Pr(>F)`[1]
tp3.div.rsrp.p
#P=0.0005

#TP4 Conv.
head(otu.meta[,1:10])
tp4.conv.rsrp <- subset(otu.meta, Timepoint=="TP4")
tp4.conv.rsrp <- subset(tp4.conv.rsrp, Soil=="Conv.")
tp4.conv.rsrp <- tp4.conv.rsrp[!(tp4.conv.rsrp$Compartment=="bulk"),]
head(tp4.conv.rsrp[,1:10])
dim(tp4.conv.rsrp)

dist.tp4.conv.rsrp <- vegdist(tp4.conv.rsrp[,-(1:5)], method="bray")
set.seed(3)
test53 <- adonis(dist.tp4.conv.rsrp~Compartment, data=tp4.conv.rsrp, strata=tp4.conv.rsrp$Block, permutation=9999)
tp4.conv.rsrp.p <- test53$aov.tab$`Pr(>F)`[1]
tp4.conv.rsrp.p
#P=0.0005

#TP4 Div.
tp4.div.rsrp <- subset(otu.meta, Timepoint=="TP4")
tp4.div.rsrp <- subset(tp4.div.rsrp, Soil=="Div.")
tp4.div.rsrp <- tp4.div.rsrp[!(tp4.div.rsrp$Compartment=="bulk"),]
head(tp4.div.rsrp[,1:10])
dim(tp4.div.rsrp)

dist.tp4.div.rsrp <- vegdist(tp4.div.rsrp[,-(1:5)], method="bray")
set.seed(3)
test54 <- adonis(dist.tp4.div.rsrp~Compartment, data=tp4.div.rsrp, strata=tp4.div.rsrp$Block, permutation=9999)
tp4.div.rsrp.p <- test54$aov.tab$`Pr(>F)`[1]
tp4.div.rsrp.p
#P=0.001

#Adjust p-values for multiple comparisons
p.adj2 <- p.adjust(c(tp1.conv.rsrp.p, tp1.div.rsrp.p, tp2.conv.rsrp.p, tp2.div.rsrp.p,
                     tp3.conv.rsrp.p, tp3.div.rsrp.p, tp4.conv.rsrp.p, tp4.div.rsrp.p), method="bonferroni", n=54)
p.adj2
#all signficant


##############
#Create graphs

# Figure 3a

# create dataframe from nmds ordination
ord.nmds.df <- as.data.frame(ord.nmds$points)
dim(ord.nmds.df)

# add sample metadata
meta <- as.data.frame(sample_data(physeq.want))
dim(meta)
ord.nmds.meta <- merge(meta, ord.nmds.df, by.x=0, by.y=0)
head(ord.nmds.meta)

#Reorder levels of compartment
levels(ord.nmds.meta$Compartment)
ord.nmds.meta$Compartment <- gsub("bulk", "Bulk soil", ord.nmds.meta$Compartment)
ord.nmds.meta$Compartment <- gsub("rhizosphere", "Rhizosphere soil", ord.nmds.meta$Compartment)
ord.nmds.meta$Compartment <- gsub("rhizoplane", "Rhizoplane", ord.nmds.meta$Compartment)
ord.nmds.meta$Compartment <- factor(ord.nmds.meta$Compartment, levels=c("Bulk soil", "Rhizosphere soil", "Rhizoplane"))

#Rename timepoint to dev stage
levels(ord.nmds.meta$Timepoint)
ord.nmds.meta$Timepoint <- revalue(ord.nmds.meta$Timepoint, c("TP1"="V4", "TP2"="V11", "TP3"="R2", "TP4"="R5")) 

#graph
head(ord.nmds.meta)
tp.nmds.plot <- ggplot(ord.nmds.meta, aes(x=MDS1, y=MDS2, color=Compartment, shape=Soil)) +
  geom_point(size=1) +
  facet_wrap(~Timepoint) +
  stat_ellipse(type="t", size=0.5, aes(linetype=Soil)) +
  labs(title="", x="NMDS1", y="NMDS2") +
  theme_bw() +
  theme(axis.title=element_text(size=8),
        axis.text=element_text(size=8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x = element_text(size = 8, margin=margin(3,3,3,3)),
        legend.position="none",
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  geom_hline(yintercept=0, color="grey69", linetype=2) +
  geom_vline(xintercept=0, color="grey69", linetype=2) +
  guides(fill=FALSE)
tp.nmds.plot

save(tp.nmds.plot, file="betadiv.tp.16s.graph.RData")

# Figure 4a

#spell out conv and div completely 
ord.nmds.meta$Soil <- revalue(sample_data(physeq.want)$Soil, c("Conv."="Conventional", "Div."="Diversified"))
levels(sample_data(physeq.want)$Soil)

#graph
prox.nmds.plot <- ggplot(ord.nmds.meta, aes(x=MDS1, y=MDS2, color=Timepoint)) +
  geom_point(size=1) +
  facet_grid(Soil~Compartment) +
  stat_ellipse(type="t", size=0.5) +
  labs(title="", x="NMDS1", y="NMDS2") +
  theme_bw() +
  theme(axis.title=element_text(size=8),
        axis.text=element_text(size=8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x = element_text(size = 8, margin=margin(3,3,3,3)),
        strip.text.y = element_text(size = 8, margin=margin(3,3,3,3)),
        legend.position="none",
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  geom_hline(yintercept=0, color="grey69", linetype=2) +
  geom_vline(xintercept=0, color="grey69", linetype=2) +
  guides(fill=FALSE)
prox.nmds.plot

save(prox.nmds.plot, file="betadiv.prox.16s.graph.RData")


##########
#ITS Data

rm(list=ls())

###################
#Create ordinations

#load phyloseq object
load("DEVITS.physeq.want.RData")
physeq.want

##############################################################################################################################
load("DEVITS.ord.nmds.RData")
#############################################################################################################################

#calculate distance using Bray-Curtis method
dist <- distance(physeq.want, method="bray")

#PCoA
ord.pcoa <- ordinate(physeq.want, "PCoA", distance=dist)

ord.pcoa.plot <- plot_ordination(physeq.want, ord.pcoa, color = "Compartment", shape = "Soil") + 
  facet_wrap(~Timepoint) +
  ggtitle("ITS PCoA using Bray-Curtis Dissimilarity") + 
  stat_ellipse(type="t")
ord.pcoa.plot 

#NMDS
ord.nmds <- ordinate(physeq.want, "NMDS", distance=dist, k=2, sratmax=0.999999)
ord.nmds
#No covergent solutions, stress = 0.200885, try another dimension
ord.nmds2 <- ordinate(physeq.want, "NMDS", distance=dist, k=3)
#No convergence
ord.nmds3 <- ordinate(physeq.want, "NMDS", distance=dist, k=5)
#No convergence
#output: *** No convergence -- monoMDS stopping criteria: 20: stress ratio > sratmax
#Try changing sratmax using a different NMDS ordination function
#plot solutions to see how they differ
meta.nmds <- metaMDS(dist, k=2, plot=TRUE)
meta.nmds
#try changing arguments
meta.nmds2 <- metaMDS(dist, k=2, plot=TRUE, maxit=300, sratmax=0.999999)
meta.nmds2
#Same results, but plots look similar, so this is likely an OK ordination

ord.nmds.plot <- plot_ordination(physeq.want, ord.nmds, color = "Compartment", shape = "Soil") +
  facet_wrap(~Timepoint) +
  ggtitle("ITS NMDS using Bray-Curtis Dissimilarity") + 
  stat_ellipse(type="t")
ord.nmds.plot
#using NMDS ordination to match 16S data


###########
#Statistics

#make a dataframe to run PERMANOVA on
adonis.df = as(sample_data(physeq.want), "data.frame")

#Full model without block effect, since it is a random effect
#constrained to block
set.seed(2)
full <- adonis(dist~Soil*Compartment*Timepoint, data = adonis.df, strata=adonis.df$Block, permutations = 9999)
full

#                             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Soil                         1     4.184  4.1837 28.0677 0.08250 0.0001 ***
#Compartment                  2     8.394  4.1969 28.1561 0.16551 0.0001 ***
#Timepoint                    3     4.809  1.6032 10.7553 0.09484 0.0001 ***
#Soil:Compartment             2     0.918  0.4590  3.0790 0.01810 0.0001 ***
#Soil:Timepoint               3     1.115  0.3717  2.4935 0.02199 0.0001 ***
#Compartment:Timepoint        6     3.329  0.5548  3.7221 0.06564 0.0001 ***
#Soil:Compartment:Timepoint   6     0.985  0.1641  1.1011 0.01942 0.1789    
#Residuals                  181    26.979  0.1491         0.53200           
#Total 


#another full model to see if significance changes
set.seed(2)
full2 <- adonis(dist~Compartment*Soil*Timepoint, data = adonis.df, strata=adonis.df$Block, permutations = 9999)
full2
#no change

#Mixed model to investigate block effect
#Block is a random effect so exclude from interactions
set.seed(2)
mixed <- adonis(dist~Soil*Compartment*Timepoint+Block, data = adonis.df, strata=adonis.df$Block, permutations = 9999)
mixed

#                             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Soil                         1     4.184  4.1837 29.1084 0.08250 0.0001 ***
#Compartment                  2     8.394  4.1969 29.2001 0.16551 0.0001 ***
#Timepoint                    3     4.809  1.6032 11.1541 0.09484 0.0001 ***
#Block                        2     1.255  0.6277  4.3670 0.02475 0.0001 ***
#Soil:Compartment             2     0.908  0.4542  3.1599 0.01791 0.0001 ***
#Soil:Timepoint               3     1.111  0.3705  2.5775 0.02192 0.0001 ***
#Compartment:Timepoint        6     3.340  0.5566  3.8727 0.06585 0.0001 ***
#Soil:Compartment:Timepoint   6     0.984  0.1640  1.1409 0.01940 0.1470    
#Residuals                  179    25.727  0.1437         0.50731           
#Total                      204    50.713                 1.00000 


#graph to try to see interaction of comp*time
ord.nmds.plot2 <- plot_ordination(physeq.want, ord.nmds, color = "Timepoint") + 
  facet_grid(Soil~Compartment) +
  ggtitle("16S NMDS using Bray-Curtis Dissimilarity") + 
  stat_ellipse(type="t")
ord.nmds.plot2 
#rhizoplane seems to change across time differently than bulk and rhizosphere


##########
#Contrasts
#To pick apart significance between communities of interest
#Will use PERMANOVA with Bonferroni correction for multiple testing

#create data frame

#OTU data
otu = as(otu_table(physeq.want), "matrix")
head(otu[,1:5])
otu.trans <- t(otu)
head(otu.trans[,1:5])
dim(otu.trans)

#meta data
meta <- sample_data(physeq.want)
head(meta)
dim(meta)

#Merge otu and metadata
otu.meta <- merge(meta, otu.trans, by=0)
head(otu.meta[,1:10])
dim(otu.meta)
rownames(otu.meta) <- otu.meta[,1]
head(otu.meta[,1:10])
otu.meta <- otu.meta[,-1]
head(otu.meta[,1:10])


#################################################
#Bulk vs rhizosphere for each time point and soil

#TP1 Conv.
head(otu.meta[,1:10])
tp1.conv.brs <- subset(otu.meta, Timepoint=="TP1")
tp1.conv.brs <- subset(tp1.conv.brs, Soil=="Conv.")
tp1.conv.brs <- tp1.conv.brs[!(tp1.conv.brs$Compartment=="rhizoplane"),]
head(tp1.conv.brs[,1:10])
dim(tp1.conv.brs)

dist.tp1.conv.brs <- vegdist(tp1.conv.brs[,-(1:5)], method="bray")
set.seed(3)
test1 <- adonis(dist.tp1.conv.brs~Compartment, data=tp1.conv.brs, strata=tp1.conv.brs$Block, permutation=9999)
tp1.conv.brs.p <- test1$aov.tab$`Pr(>F)`[1]
tp1.conv.brs.p
#P=0.0032

#TP1 Div.
tp1.div.brs <- subset(otu.meta, Timepoint=="TP1")
tp1.div.brs <- subset(tp1.div.brs, Soil=="Div.")
tp1.div.brs <- tp1.div.brs[!(tp1.div.brs$Compartment=="rhizoplane"),]
head(tp1.div.brs[,1:10])
dim(tp1.div.brs)

dist.tp1.div.brs <- vegdist(tp1.div.brs[,-(1:5)], method="bray")
set.seed(3)
test2 <- adonis(dist.tp1.div.brs~Compartment, data=tp1.div.brs, strata=tp1.div.brs$Block, permutation=9999)
tp1.div.brs.p <- test2$aov.tab$`Pr(>F)`[1]
tp1.div.brs.p
#P=0.005

#TP2 Conv.
tp2.conv.brs <- subset(otu.meta, Timepoint=="TP2")
tp2.conv.brs <- subset(tp2.conv.brs, Soil=="Conv.")
tp2.conv.brs <- tp2.conv.brs[!(tp2.conv.brs$Compartment=="rhizoplane"),]
head(tp2.conv.brs[,1:10])
dim(tp2.conv.brs)

dist.tp2.conv.brs <- vegdist(tp2.conv.brs[,-(1:5)], method="bray")
set.seed(3)
test3 <- adonis(dist.tp2.conv.brs~Compartment, data=tp2.conv.brs, strata=tp2.conv.brs$Block, permutation=9999)
tp2.conv.brs.p <- test3$aov.tab$`Pr(>F)`[1]
tp2.conv.brs.p
#P=0.0023

#TP2 Div.
tp2.div.brs <- subset(otu.meta, Timepoint=="TP2")
tp2.div.brs <- subset(tp2.div.brs, Soil=="Div.")
tp2.div.brs <- tp2.div.brs[!(tp2.div.brs$Compartment=="rhizoplane"),]
head(tp2.div.brs[,1:10])
dim(tp2.div.brs)

dist.tp2.div.brs <- vegdist(tp2.div.brs[,-(1:5)], method="bray")
set.seed(3)
test4 <- adonis(dist.tp2.div.brs~Compartment, data=tp2.div.brs, strata=tp2.div.brs$Block, permutation=9999)
tp2.div.brs.p <- test4$aov.tab$`Pr(>F)`[1]
tp2.div.brs.p
#P=0.0977

#TP3 Conv.
head(otu.meta[,1:10])
tp3.conv.brs <- subset(otu.meta, Timepoint=="TP3")
tp3.conv.brs <- subset(tp3.conv.brs, Soil=="Conv.")
tp3.conv.brs <- tp3.conv.brs[!(tp3.conv.brs$Compartment=="rhizoplane"),]
head(tp3.conv.brs[,1:10])
dim(tp3.conv.brs)

dist.tp3.conv.brs <- vegdist(tp3.conv.brs[,-(1:5)], method="bray")
set.seed(3)
test5 <- adonis(dist.tp3.conv.brs~Compartment, data=tp3.conv.brs, strata=tp3.conv.brs$Block, permutation=9999)
tp3.conv.brs.p <- test5$aov.tab$`Pr(>F)`[1]
tp3.conv.brs.p
#P=0.0002

#TP3 Div.
tp3.div.brs <- subset(otu.meta, Timepoint=="TP3")
tp3.div.brs <- subset(tp3.div.brs, Soil=="Div.")
tp3.div.brs <- tp3.div.brs[!(tp3.div.brs$Compartment=="rhizoplane"),]
head(tp3.div.brs[,1:10])
dim(tp3.div.brs)

dist.tp3.div.brs <- vegdist(tp3.div.brs[,-(1:5)], method="bray")
set.seed(3)
test6 <- adonis(dist.tp3.div.brs~Compartment, data=tp3.div.brs, strata=tp3.div.brs$Block, permutation=9999)
tp3.div.brs.p <- test6$aov.tab$`Pr(>F)`[1]
tp3.div.brs.p
#P=0.0062

#TP4 Conv.
head(otu.meta[,1:10])
tp4.conv.brs <- subset(otu.meta, Timepoint=="TP4")
tp4.conv.brs <- subset(tp4.conv.brs, Soil=="Conv.")
tp4.conv.brs <- tp4.conv.brs[!(tp4.conv.brs$Compartment=="rhizoplane"),]
head(tp4.conv.brs[,1:10])
dim(tp4.conv.brs)

dist.tp4.conv.brs <- vegdist(tp4.conv.brs[,-(1:5)], method="bray")
set.seed(3)
test7 <- adonis(dist.tp4.conv.brs~Compartment, data=tp4.conv.brs, strata=tp4.conv.brs$Block, permutation=9999)
tp4.conv.brs.p <- test7$aov.tab$`Pr(>F)`[1]
tp4.conv.brs.p
#P=0.0008

#TP4 Div.
tp4.div.brs <- subset(otu.meta, Timepoint=="TP4")
tp4.div.brs <- subset(tp4.div.brs, Soil=="Div.")
tp4.div.brs <- tp4.div.brs[!(tp4.div.brs$Compartment=="rhizoplane"),]
head(tp4.div.brs[,1:10])
dim(tp4.div.brs)

dist.tp4.div.brs <- vegdist(tp4.div.brs[,-(1:5)], method="bray")
set.seed(3)
test8 <- adonis(dist.tp4.div.brs~Compartment, data=tp4.div.brs, strata=tp4.div.brs$Block, permutation=9999)
tp4.div.brs.p <- test8$aov.tab$`Pr(>F)`[1]
tp4.div.brs.p
#P=0.0005

#Adjust for multiple comparisons
#46 is the total number of contrasts done in this script
p.adj1 <- p.adjust(c(tp1.conv.brs.p, tp1.div.brs.p, tp2.conv.brs.p, tp2.div.brs.p,
                     tp3.conv.brs.p, tp3.div.brs.p, tp4.conv.brs.p, tp4.div.brs.p), method="bonferroni", n=54)
p.adj1
#Only TP3 Conv, TP4 both significant


###############################################
#Bulk vs rhizoplane for each timepoint and soil

#TP1 Conv.
head(otu.meta[,1:10])
tp1.conv.brp <- subset(otu.meta, Timepoint=="TP1")
tp1.conv.brp <- subset(tp1.conv.brp, Soil=="Conv.")
tp1.conv.brp <- tp1.conv.brp[!(tp1.conv.brp$Compartment=="rhizosphere"),]
head(tp1.conv.brp[,1:10])
dim(tp1.conv.brp)

dist.tp1.conv.brp <- vegdist(tp1.conv.brp[,-(1:5)], method="bray")
set.seed(3)
test9 <- adonis(dist.tp1.conv.brp~Compartment, data=tp1.conv.brp, strata=tp1.conv.brp$Block, permutation=9999)
tp1.conv.brp.p <- test9$aov.tab$`Pr(>F)`[1]
tp1.conv.brp.p
#P=0.0001

#TP1 Div.
tp1.div.brp <- subset(otu.meta, Timepoint=="TP1")
tp1.div.brp <- subset(tp1.div.brp, Soil=="Div.")
tp1.div.brp <- tp1.div.brp[!(tp1.div.brp$Compartment=="rhizosphere"),]
head(tp1.div.brp[,1:10])
dim(tp1.div.brp)

dist.tp1.div.brp <- vegdist(tp1.div.brp[,-(1:5)], method="bray")
set.seed(3)
test10 <- adonis(dist.tp1.div.brp~Compartment, data=tp1.div.brp, strata=tp1.div.brp$Block, permutation=9999)
tp1.div.brp.p <- test10$aov.tab$`Pr(>F)`[1]
tp1.div.brp.p
#P=0.0004

#TP2 Conv.
tp2.conv.brp <- subset(otu.meta, Timepoint=="TP2")
tp2.conv.brp <- subset(tp2.conv.brp, Soil=="Conv.")
tp2.conv.brp <- tp2.conv.brp[!(tp2.conv.brp$Compartment=="rhizosphere"),]
head(tp2.conv.brp[,1:10])
dim(tp2.conv.brp)

dist.tp2.conv.brp <- vegdist(tp2.conv.brp[,-(1:5)], method="bray")
set.seed(3)
test11 <- adonis(dist.tp2.conv.brp~Compartment, data=tp2.conv.brp, strata=tp2.conv.brp$Block, permutation=9999)
tp2.conv.brp.p <- test11$aov.tab$`Pr(>F)`[1]
tp2.conv.brp.p
#P=0.0007, significant

#TP2 Div.
tp2.div.brp <- subset(otu.meta, Timepoint=="TP2")
tp2.div.brp <- subset(tp2.div.brp, Soil=="Div.")
tp2.div.brp <- tp2.div.brp[!(tp2.div.brp$Compartment=="rhizosphere"),]
head(tp2.div.brp[,1:10])
dim(tp2.div.brp)

dist.tp2.div.brp <- vegdist(tp2.div.brp[,-(1:5)], method="bray")
set.seed(3)
test12 <- adonis(dist.tp2.div.brp~Compartment, data=tp2.div.brp, strata=tp2.div.brp$Block, permutation=9999)
tp2.div.brp.p <- test12$aov.tab$`Pr(>F)`[1]
tp2.div.brp.p
#P=0.0008

#TP3 Conv.
head(otu.meta[,1:10])
tp3.conv.brp <- subset(otu.meta, Timepoint=="TP3")
tp3.conv.brp <- subset(tp3.conv.brp, Soil=="Conv.")
tp3.conv.brp <- tp3.conv.brp[!(tp3.conv.brp$Compartment=="rhizosphere"),]
head(tp3.conv.brp[,1:10])
dim(tp3.conv.brp)

dist.tp3.conv.brp <- vegdist(tp3.conv.brp[,-(1:5)], method="bray")
set.seed(3)
test13 <- adonis(dist.tp3.conv.brp~Compartment, data=tp3.conv.brp, strata=tp3.conv.brp$Block, permutation=9999)
tp3.conv.brp.p <- test13$aov.tab$`Pr(>F)`[1]
tp3.conv.brp.p
#P=0.0006

#TP3 Div.
tp3.div.brp <- subset(otu.meta, Timepoint=="TP3")
tp3.div.brp <- subset(tp3.div.brp, Soil=="Div.")
tp3.div.brp <- tp3.div.brp[!(tp3.div.brp$Compartment=="rhizosphere"),]
head(tp3.div.brp[,1:10])
dim(tp3.div.brp)

dist.tp3.div.brp <- vegdist(tp3.div.brp[,-(1:5)], method="bray")
set.seed(3)
test14 <- adonis(dist.tp3.div.brp~Compartment, data=tp3.div.brp, strata=tp3.div.brp$Block, permutation=9999)
tp3.div.brp.p <- test14$aov.tab$`Pr(>F)`[1]
tp3.div.brp.p
#P=0.0005

#TP4 Conv.
head(otu.meta[,1:10])
tp4.conv.brp <- subset(otu.meta, Timepoint=="TP4")
tp4.conv.brp <- subset(tp4.conv.brp, Soil=="Conv.")
tp4.conv.brp <- tp4.conv.brp[!(tp4.conv.brp$Compartment=="rhizosphere"),]
head(tp4.conv.brp[,1:10])
dim(tp4.conv.brp)

dist.tp4.conv.brp <- vegdist(tp4.conv.brp[,-(1:5)], method="bray")
set.seed(3)
test15 <- adonis(dist.tp4.conv.brp~Compartment, data=tp4.conv.brp, strata=tp4.conv.brp$Block, permutation=9999)
tp4.conv.brp.p <- test15$aov.tab$`Pr(>F)`[1]
tp4.conv.brp.p
#P=0.0005

#TP4 Div.
tp4.div.brp <- subset(otu.meta, Timepoint=="TP4")
tp4.div.brp <- subset(tp4.div.brp, Soil=="Div.")
tp4.div.brp <- tp4.div.brp[!(tp4.div.brp$Compartment=="rhizosphere"),]
head(tp4.div.brp[,1:10])
dim(tp4.div.brp)

dist.tp4.div.brp <- vegdist(tp4.div.brp[,-(1:5)], method="bray")
set.seed(3)
test16 <- adonis(dist.tp4.div.brp~Compartment, data=tp4.div.brp, strata=tp4.div.brp$Block, permutation=9999)
tp4.div.brp.p <- test16$aov.tab$`Pr(>F)`[1]
tp4.div.brp.p
#P=0.0008

#Adjust p-values for multiple comparisons
p.adj2 <- p.adjust(c(tp1.conv.brp.p, tp1.div.brp.p, tp2.conv.brp.p, tp2.div.brp.p,
                     tp3.conv.brp.p, tp3.div.brp.p, tp4.conv.brp.p, tp4.div.brp.p), method="bonferroni", n=54)
p.adj2
#all signficant


###############################
#Bulk vs bulk for each timepoint

#TP1
tp1.b <- subset(otu.meta, Timepoint=="TP1")
tp1.b <- subset(tp1.b, Compartment=="bulk")
head(tp1.b[,1:5])
dim(tp1.b)

dist.tp1.b <- vegdist(tp1.b[,-(1:5)], method="bray")
set.seed(3)
test17 <- adonis(dist.tp1.b~Soil, data=tp1.b, strata=tp1.b$Block, permutation=9999)
tp1.b.p <- test17$aov.tab$`Pr(>F)`[1]
tp1.b.p
#P=0.0005

#TP2
tp2.b <- subset(otu.meta, Timepoint=="TP2")
tp2.b <- subset(tp2.b, Compartment=="bulk")
head(tp2.b[,1:5])
dim(tp2.b)

dist.tp2.b <- vegdist(tp2.b[,-(1:5)], method="bray")
set.seed(3)
test18 <- adonis(dist.tp2.b~Soil, data=tp2.b, strata=tp2.b$Block, permutation=9999)
tp2.b.p <- test18$aov.tab$`Pr(>F)`[1]
tp2.b.p
#P=0.0003

#TP3
tp3.b <- subset(otu.meta, Timepoint=="TP3")
tp3.b <- subset(tp3.b, Compartment=="bulk")
head(tp3.b[,1:5])
dim(tp3.b)

dist.tp3.b <- vegdist(tp3.b[,-(1:5)], method="bray")
set.seed(3)
test19 <- adonis(dist.tp3.b~Soil, data=tp3.b, strata=tp3.b$Block, permutation=9999)
tp3.b.p <- test19$aov.tab$`Pr(>F)`[1]
tp3.b.p
#P=0.0005

#TP4
tp4.b <- subset(otu.meta, Timepoint=="TP4")
tp4.b <- subset(tp4.b, Compartment=="bulk")
head(tp4.b[,1:5])
dim(tp4.b)

dist.tp4.b <- vegdist(tp4.b[,-(1:5)], method="bray")
set.seed(3)
test20 <- adonis(dist.tp4.b~Soil, data=tp4.b, strata=tp4.b$Block, permutation=9999)
tp4.b.p <- test20$aov.tab$`Pr(>F)`[1]
tp4.b.p
#P=0.0005

#Mulitple comparisons p adjustment
p.adj3 <- p.adjust(c(tp1.b.p, tp2.b.p, tp3.b.p, tp4.b.p), method="bonferroni", n=54)
p.adj3
#all significant


############################################
#Rhizosphere vs rhizosphere at each timepoint

#TP1
tp1.rs <- subset(otu.meta, Timepoint=="TP1")
tp1.rs <- subset(tp1.rs, Compartment=="rhizosphere")
head(tp1.rs[,1:5])
dim(tp1.rs)

dist.tp1.rs <- vegdist(tp1.rs[,-(1:5)], method="bray")
set.seed(3)
test21 <- adonis(dist.tp1.rs~Soil, data=tp1.rs, strata=tp1.rs$Block, permutation=9999)
tp1.rs.p <- test21$aov.tab$`Pr(>F)`[1]
tp1.rs.p
#P=0.0005

#TP2
tp2.rs <- subset(otu.meta, Timepoint=="TP2")
tp2.rs <- subset(tp2.rs, Compartment=="rhizosphere")
head(tp2.rs[,1:5])
dim(tp2.rs)

dist.tp2.rs <- vegdist(tp2.rs[,-(1:5)], method="bray")
set.seed(3)
test22 <- adonis(dist.tp2.rs~Soil, data=tp2.rs, strata=tp2.rs$Block, permutation=9999)
tp2.rs.p <- test22$aov.tab$`Pr(>F)`[1]
tp2.rs.p
#P=0.0005

#TP3
tp3.rs <- subset(otu.meta, Timepoint=="TP3")
tp3.rs <- subset(tp3.rs, Compartment=="rhizosphere")
head(tp3.rs[,1:5])
dim(tp3.rs)

dist.tp3.rs <- vegdist(tp3.rs[,-(1:5)], method="bray")
set.seed(3)
test23 <- adonis(dist.tp3.rs~Soil, data=tp3.rs, strata=tp3.rs$Block, permutation=9999)
tp3.rs.p <- test23$aov.tab$`Pr(>F)`[1]
tp3.rs.p
#P=0.0002

#TP4
tp4.rs <- subset(otu.meta, Timepoint=="TP4")
tp4.rs <- subset(tp4.rs, Compartment=="rhizosphere")
head(tp4.rs[,1:5])
dim(tp4.rs)

set.seed(3)
dist.tp4.rs <- vegdist(tp4.rs[,-(1:5)], method="bray")
test24 <- adonis(dist.tp4.rs~Soil, data=tp4.rs, strata=tp4.rs$Block, permutation=9999)
tp4.rs.p <- test24$aov.tab$`Pr(>F)`[1]
tp4.rs.p
#P=0.0005, significant

#Multiple comparisons p-adjustment
p.adj4 <- p.adjust(c(tp1.rs.p, tp2.rs.p, tp3.rs.p, tp4.rs.p), method="bonferroni", n=54)
p.adj4
#all significant


############################################
#Rhizoplane vs rhizoplane at each timepoint

#TP1
tp1.rp <- subset(otu.meta, Timepoint=="TP1")
tp1.rp <- subset(tp1.rp, Compartment=="rhizoplane")
head(tp1.rp[,1:5])
dim(tp1.rp)

dist.tp1.rp <- vegdist(tp1.rp[,-(1:5)], method="bray")
set.seed(3)
test25 <- adonis(dist.tp1.rp~Soil, data=tp1.rp, strata=tp1.rp$Block, permutation=9999)
tp1.rp.p <- test25$aov.tab$`Pr(>F)`[1]
tp1.rp.p
#P=0.0006

#TP2
tp2.rp <- subset(otu.meta, Timepoint=="TP2")
tp2.rp <- subset(tp2.rp, Compartment=="rhizoplane")
head(tp2.rp[,1:5])
dim(tp2.rp)

dist.tp2.rp <- vegdist(tp2.rp[,-(1:5)], method="bray")
set.seed(3)
test26 <- adonis(dist.tp2.rp~Soil, data=tp2.rp, strata=tp2.rp$Block, permutation=9999)
tp2.rp.p <- test26$aov.tab$`Pr(>F)`[1]
tp2.rp.p
#P=0.0021

#TP3
tp3.rp <- subset(otu.meta, Timepoint=="TP3")
tp3.rp <- subset(tp3.rp, Compartment=="rhizoplane")
head(tp3.rp[,1:5])
dim(tp3.rp)

dist.tp3.rp <- vegdist(tp3.rp[,-(1:5)], method="bray")
set.seed(3)
test27 <- adonis(dist.tp3.rp~Soil, data=tp3.rp, strata=tp3.rp$Block, permutation=9999)
tp3.rp.p <- test27$aov.tab$'Pr(>F)'[1]
tp3.rp.p
#P=0.0017

#TP4
tp4.rp <- subset(otu.meta, Timepoint=="TP4")
tp4.rp <- subset(tp4.rp, Compartment=="rhizoplane")
head(tp4.rp[,1:5])
dim(tp4.rp)

dist.tp4.rp <- vegdist(tp4.rp[,-(1:5)], method="bray")
set.seed(3)
test28 <- adonis(dist.tp4.rp~Soil, data=tp4.rp, strata=tp4.rp$Block, permutation=9999)
tp4.rp.p <- test28$aov.tab$'Pr(>F)'[1]
tp4.rp.p
#P=0.0008

#Multiple comparisons p-adjustment
p.adj5 <- p.adjust(c(tp1.rp.p, tp2.rp.p, tp3.rp.p, tp4.rp.p), method="bonferroni", n=54)
p.adj5
#TP1 and TP4 significant


#####################################
#Compartment x time point interaction

#Explore CompxTime interaction, most biologically significant interaction with visual support in ordinations


##########
#Bulk soil

#TP1 vs TP2 Conv.
conv.b.tp12 <- subset(otu.meta, Compartment=="bulk")
conv.b.tp12 <- subset(conv.b.tp12, Soil=="Conv.")
conv.b.tp12 <- conv.b.tp12[!(conv.b.tp12$Timepoint=="TP3"),]
conv.b.tp12 <- conv.b.tp12[!(conv.b.tp12$Timepoint=="TP4"),]  
head(conv.b.tp12[,1:10])

dist.conv.b.tp12 <- vegdist(conv.b.tp12[,-(1:5)], method="bray")
set.seed(3)
test29 <- adonis(dist.conv.b.tp12~Timepoint, data=conv.b.tp12, strata=conv.b.tp12$Block, permutation=9999)
conv.b.tp12.p <- test29$aov.tab$'Pr(>F)'[1]
conv.b.tp12.p
#P=0.0219

#TP1 vs TP2 Div.
div.b.tp12 <- subset(otu.meta, Compartment=="bulk")
div.b.tp12 <- subset(div.b.tp12, Soil=="Div.")
div.b.tp12 <- div.b.tp12[!(div.b.tp12$Timepoint=="TP3"),]
div.b.tp12 <- div.b.tp12[!(div.b.tp12$Timepoint=="TP4"),]  
head(div.b.tp12[,1:10])

dist.div.b.tp12 <- vegdist(div.b.tp12[,-(1:5)], method="bray")
set.seed(3)
test30 <- adonis(dist.div.b.tp12~Timepoint, data=div.b.tp12, strata=div.b.tp12$Block, permutation=9999)
div.b.tp12.p <- test30$aov.tab$'Pr(>F)'[1]
div.b.tp12.p
#0.0059

#TP2 vs TP3 Conv.
conv.b.tp23 <- subset(otu.meta, Compartment=="bulk")
conv.b.tp23 <- subset(conv.b.tp23, Soil=="Conv.")
conv.b.tp23 <- conv.b.tp23[!(conv.b.tp23$Timepoint=="TP1"),]
conv.b.tp23 <- conv.b.tp23[!(conv.b.tp23$Timepoint=="TP4"),]  
head(conv.b.tp23[,1:10])

dist.conv.b.tp23 <- vegdist(conv.b.tp23[,-(1:5)], method="bray")
set.seed(3)
test31 <- adonis(dist.conv.b.tp23~Timepoint, data=conv.b.tp23, strata=conv.b.tp23$Block, permutation=9999)
conv.b.tp23.p <- test31$aov.tab$'Pr(>F)'[1]
conv.b.tp23.p
#P=0.0061

#TP2 vs TP3 Div.
div.b.tp23 <- subset(otu.meta, Compartment=="bulk")
div.b.tp23 <- subset(div.b.tp23, Soil=="Div.")
div.b.tp23 <- div.b.tp23[!(div.b.tp23$Timepoint=="TP1"),]
div.b.tp23 <- div.b.tp23[!(div.b.tp23$Timepoint=="TP4"),]  
head(div.b.tp23[,1:10])

dist.div.b.tp23 <- vegdist(div.b.tp23[,-(1:5)], method="bray")
set.seed(3)
test32 <- adonis(dist.div.b.tp23~Timepoint, data=div.b.tp23, strata=div.b.tp23$Block, permutation=9999)
div.b.tp23.p <- test32$aov.tab$'Pr(>F)'[1]
div.b.tp23.p
#P=0.5224

#TP3 vs TP4 Conv.
conv.b.tp34 <- subset(otu.meta, Compartment=="bulk")
conv.b.tp34 <- subset(conv.b.tp34, Soil=="Conv.")
conv.b.tp34 <- conv.b.tp34[!(conv.b.tp34$Timepoint=="TP1"),]
conv.b.tp34 <- conv.b.tp34[!(conv.b.tp34$Timepoint=="TP2"),]  
head(conv.b.tp34[,1:10])

dist.conv.b.tp34 <- vegdist(conv.b.tp34[,-(1:5)], method="bray")
set.seed(3)
test33 <- adonis(dist.conv.b.tp34~Timepoint, data=conv.b.tp34, strata=conv.b.tp34$Block, permutation=9999)
conv.b.tp34.p <- test33$aov.tab$'Pr(>F)'[1]
conv.b.tp34.p
#0.0287

#TP3 vs TP4 Div.
div.b.tp34 <- subset(otu.meta, Compartment=="bulk")
div.b.tp34 <- subset(div.b.tp34, Soil=="Div.")
div.b.tp34 <- div.b.tp34[!(div.b.tp34$Timepoint=="TP1"),]
div.b.tp34 <- div.b.tp34[!(div.b.tp34$Timepoint=="TP2"),]  
head(div.b.tp34[,1:10])

dist.div.b.tp34 <- vegdist(div.b.tp34[,-(1:5)], method="bray")
set.seed(3)
test34 <- adonis(dist.div.b.tp34~Timepoint, data=div.b.tp34, strata=div.b.tp34$Block, permutation=9999)
div.b.tp34.p <- test34$aov.tab$'Pr(>F)'[1]
div.b.tp34.p
#0.6787

#Multiple comparisons p adjustment
padj6 <- p.adjust(c(conv.b.tp12.p, div.b.tp12.p, conv.b.tp23.p, div.b.tp23.p, conv.b.tp34.p, div.b.tp34.p),
                  method="bonferroni", n=54)
padj6
#no significance


#################
#Rhizosphere soil

#TP1 vs TP2 Conv.
conv.rs.tp12 <- subset(otu.meta, Compartment=="rhizosphere")
conv.rs.tp12 <- subset(conv.rs.tp12, Soil=="Conv.")
conv.rs.tp12 <- conv.rs.tp12[!(conv.rs.tp12$Timepoint=="TP3"),]
conv.rs.tp12 <- conv.rs.tp12[!(conv.rs.tp12$Timepoint=="TP4"),]  
head(conv.rs.tp12[,1:10])

dist.conv.rs.tp12 <- vegdist(conv.rs.tp12[,-(1:5)], method="bray")
set.seed(3)
test35 <- adonis(dist.conv.rs.tp12~Timepoint, data=conv.rs.tp12, strata=conv.rs.tp12$Block, permutation=9999)
conv.rs.tp12.p <- test35$aov.tab$'Pr(>F)'[1]
conv.rs.tp12.p
#0.0008

#TP1 vs TP2 Div.
div.rs.tp12 <- subset(otu.meta, Compartment=="rhizosphere")
div.rs.tp12 <- subset(div.rs.tp12, Soil=="Div.")
div.rs.tp12 <- div.rs.tp12[!(div.rs.tp12$Timepoint=="TP3"),]
div.rs.tp12 <- div.rs.tp12[!(div.rs.tp12$Timepoint=="TP4"),]  
head(div.rs.tp12[,1:10])

dist.div.rs.tp12 <- vegdist(div.rs.tp12[,-(1:5)], method="bray")
set.seed(3)
test36 <- adonis(dist.div.rs.tp12~Timepoint, data=div.rs.tp12, strata=div.rs.tp12$Block, permutation=9999)
div.rs.tp12.p <- test36$aov.tab$'Pr(>F)'[1]
div.rs.tp12.p
#0.0024

#TP2 vs TP3 Conv.
conv.rs.tp23 <- subset(otu.meta, Compartment=="rhizosphere")
conv.rs.tp23 <- subset(conv.rs.tp23, Soil=="Conv.")
conv.rs.tp23 <- conv.rs.tp23[!(conv.rs.tp23$Timepoint=="TP1"),]
conv.rs.tp23 <- conv.rs.tp23[!(conv.rs.tp23$Timepoint=="TP4"),]  
head(conv.rs.tp23[,1:10])

dist.conv.rs.tp23 <- vegdist(conv.rs.tp23[,-(1:5)], method="bray")
set.seed(3)
test37 <- adonis(dist.conv.rs.tp23~Timepoint, data=conv.rs.tp23, strata=conv.rs.tp23$Block, permutation=9999)
conv.rs.tp23.p <- test37$aov.tab$'Pr(>F)'[1]
conv.rs.tp23.p
#0.0005

#TP2 vs TP3 Div.
div.rs.tp23 <- subset(otu.meta, Compartment=="rhizosphere")
div.rs.tp23 <- subset(div.rs.tp23, Soil=="Div.")
div.rs.tp23 <- div.rs.tp23[!(div.rs.tp23$Timepoint=="TP1"),]
div.rs.tp23 <- div.rs.tp23[!(div.rs.tp23$Timepoint=="TP4"),]  
head(div.rs.tp23[,1:10])

dist.div.rs.tp23 <- vegdist(div.rs.tp23[,-(1:5)], method="bray")
set.seed(3)
test38 <- adonis(dist.div.rs.tp23~Timepoint, data=div.rs.tp23, strata=div.rs.tp23$Block, permutation=9999)
div.rs.tp23.p <- test38$aov.tab$'Pr(>F)'[1]
div.rs.tp23.p
#0.0001

#TP3 vs TP4 Conv.
conv.rs.tp34 <- subset(otu.meta, Compartment=="rhizosphere")
conv.rs.tp34 <- subset(conv.rs.tp34, Soil=="Conv.")
conv.rs.tp34 <- conv.rs.tp34[!(conv.rs.tp34$Timepoint=="TP1"),]
conv.rs.tp34 <- conv.rs.tp34[!(conv.rs.tp34$Timepoint=="TP2"),]  
head(conv.rs.tp34[,1:10])

dist.conv.rs.tp34 <- vegdist(conv.rs.tp34[,-(1:5)], method="bray")
set.seed(3)
test39 <- adonis(dist.conv.rs.tp34~Timepoint, data=conv.rs.tp34, strata=conv.rs.tp34$Block, permutation=9999)
conv.rs.tp34.p <- test39$aov.tab$'Pr(>F)'[1]
conv.rs.tp34.p
#0.0048

#TP3 vs TP4 Div.
div.rs.tp34 <- subset(otu.meta, Compartment=="rhizosphere")
div.rs.tp34 <- subset(div.rs.tp34, Soil=="Div.")
div.rs.tp34 <- div.rs.tp34[!(div.rs.tp34$Timepoint=="TP1"),]
div.rs.tp34 <- div.rs.tp34[!(div.rs.tp34$Timepoint=="TP2"),]  
head(div.rs.tp34[,1:10])

dist.div.rs.tp34 <- vegdist(div.rs.tp34[,-(1:5)], method="bray")
set.seed(3)
test40 <- adonis(dist.div.rs.tp34~Timepoint, data=div.rs.tp34, strata=div.rs.tp34$Block, permutation=9999)
div.rs.tp34.p <- test40$aov.tab$'Pr(>F)'[1]
div.rs.tp34.p
#0.0005

#Multiple comparisons p adjustment
p.adj7 <- p.adjust(c(conv.rs.tp12.p, div.rs.tp12.p, conv.rs.tp23.p, div.rs.tp23.p, conv.rs.tp34.p, div.rs.tp34.p),
                   method="bonferroni", n=54)
p.adj7
#TP1 Div, TP4 Conv not significant


################
#Rhizoplane soil

#TP1 vs TP2 Conv.
conv.rp.tp12 <- subset(otu.meta, Compartment=="rhizoplane")
conv.rp.tp12 <- subset(conv.rp.tp12, Soil=="Conv.")
conv.rp.tp12 <- conv.rp.tp12[!(conv.rp.tp12$Timepoint=="TP3"),]
conv.rp.tp12 <- conv.rp.tp12[!(conv.rp.tp12$Timepoint=="TP4"),]  
head(conv.rp.tp12[,1:10])

dist.conv.rp.tp12 <- vegdist(conv.rp.tp12[,-(1:5)], method="bray")
set.seed(3)
test41 <- adonis(dist.conv.rp.tp12~Timepoint, data=conv.rp.tp12, strata=conv.rp.tp12$Block, permutation=9999)
conv.rp.tp12.p <- test41$aov.tab$'Pr(>F)'[1]
conv.rp.tp12.p
#P=0.0003

#TP1 vs TP2 Div.
div.rp.tp12 <- subset(otu.meta, Compartment=="rhizoplane")
div.rp.tp12 <- subset(div.rp.tp12, Soil=="Div.")
div.rp.tp12 <- div.rp.tp12[!(div.rp.tp12$Timepoint=="TP3"),]
div.rp.tp12 <- div.rp.tp12[!(div.rp.tp12$Timepoint=="TP4"),]  
head(div.rp.tp12[,1:10])

dist.div.rp.tp12 <- vegdist(div.rp.tp12[,-(1:5)], method="bray")
set.seed(3)
test42 <- adonis(dist.div.rp.tp12~Timepoint, data=div.rp.tp12, strata=div.rp.tp12$Block, permutation=9999)
div.rp.tp12.p <- test42$aov.tab$'Pr(>F)'[1]
div.rp.tp12.p
#P=0.0012

#TP2 vs TP3 Conv.
conv.rp.tp23 <- subset(otu.meta, Compartment=="rhizoplane")
conv.rp.tp23 <- subset(conv.rp.tp23, Soil=="Conv.")
conv.rp.tp23 <- conv.rp.tp23[!(conv.rp.tp23$Timepoint=="TP1"),]
conv.rp.tp23 <- conv.rp.tp23[!(conv.rp.tp23$Timepoint=="TP4"),]  
head(conv.rp.tp23[,1:10])

dist.conv.rp.tp23 <- vegdist(conv.rp.tp23[,-(1:5)], method="bray")
set.seed(3)
test43 <- adonis(dist.conv.rp.tp23~Timepoint, data=conv.rp.tp23, strata=conv.rp.tp23$Block, permutation=9999)
conv.rp.tp23.p <- test43$aov.tab$'Pr(>F)'[1]
conv.rp.tp23.p
#P=0.0015

#TP2 vs TP3 Div.
div.rp.tp23 <- subset(otu.meta, Compartment=="rhizoplane")
div.rp.tp23 <- subset(div.rp.tp23, Soil=="Div.")
div.rp.tp23 <- div.rp.tp23[!(div.rp.tp23$Timepoint=="TP1"),]
div.rp.tp23 <- div.rp.tp23[!(div.rp.tp23$Timepoint=="TP4"),]  
head(div.rp.tp23[,1:10])

dist.div.rp.tp23 <- vegdist(div.rp.tp23[,-(1:5)], method="bray")
set.seed(3)
test44 <- adonis(dist.div.rp.tp23~Timepoint, data=div.rp.tp23, strata=div.rp.tp23$Block, permutation=9999)
div.rp.tp23.p <- test44$aov.tab$'Pr(>F)'[1]
div.rp.tp23.p
#P=0.0005

#TP3 vs TP4 Conv.
conv.rp.tp34 <- subset(otu.meta, Compartment=="rhizoplane")
conv.rp.tp34 <- subset(conv.rp.tp34, Soil=="Conv.")
conv.rp.tp34 <- conv.rp.tp34[!(conv.rp.tp34$Timepoint=="TP1"),]
conv.rp.tp34 <- conv.rp.tp34[!(conv.rp.tp34$Timepoint=="TP2"),]  
head(conv.rp.tp34[,1:10])

dist.conv.rp.tp34 <- vegdist(conv.rp.tp34[,-(1:5)], method="bray")
set.seed(3)
test45 <- adonis(dist.conv.rp.tp34~Timepoint, data=conv.rp.tp34, strata=conv.rp.tp34$Block, permutation=9999)
conv.rp.tp34.p <- test45$aov.tab$'Pr(>F)'[1]
conv.rp.tp34.p
#P=0.0001

#TP3 vs TP4 Div.
div.rp.tp34 <- subset(otu.meta, Compartment=="rhizoplane")
div.rp.tp34 <- subset(div.rp.tp34, Soil=="Div.")
div.rp.tp34 <- div.rp.tp34[!(div.rp.tp34$Timepoint=="TP1"),]
div.rp.tp34 <- div.rp.tp34[!(div.rp.tp34$Timepoint=="TP2"),]  
head(div.rp.tp34[,1:10])

dist.div.rp.tp34 <- vegdist(div.rp.tp34[,-(1:5)], method="bray")
set.seed(3)
test46 <- adonis(dist.div.rp.tp34~Timepoint, data=div.rp.tp34, strata=div.rp.tp34$Block, permutation=9999)
div.rp.tp34.p <- test46$aov.tab$'Pr(>F)'[1]
div.rp.tp34.p
#P=0.0008

#Correct for multiple comparisons
padj7 <- p.adjust(c(conv.rp.tp12.p, div.rp.tp12.p, conv.rp.tp23.p, div.rp.tp23.p, conv.rp.tp34.p, div.rp.tp34.p),
                  method="bonferroni", n=54)
padj7
#All sig except TP12 Div and TP23 Conv


###################################################
#Rhizosphere vs Rhizoplane for each cropping system

#adding this in after all other comparisons, to make analysis more complete

#TP1 Conv.
head(otu.meta[,1:10])
tp1.conv.rsrp <- subset(otu.meta, Timepoint=="TP1")
tp1.conv.rsrp <- subset(tp1.conv.rsrp, Soil=="Conv.")
tp1.conv.rsrp <- tp1.conv.rsrp[!(tp1.conv.rsrp$Compartment=="bulk"),]
head(tp1.conv.rsrp[,1:10])
dim(tp1.conv.rsrp)

dist.tp1.conv.rsrp <- vegdist(tp1.conv.rsrp[,-(1:5)], method="bray")
set.seed(3)
test47 <- adonis(dist.tp1.conv.rsrp~Compartment, data=tp1.conv.rsrp, strata=tp1.conv.rsrp$Block, permutation=9999)
tp1.conv.rsrp.p <- test47$aov.tab$`Pr(>F)`[1]
tp1.conv.rsrp.p
#P=0.0004

#TP1 Div.
tp1.div.rsrp <- subset(otu.meta, Timepoint=="TP1")
tp1.div.rsrp <- subset(tp1.div.rsrp, Soil=="Div.")
tp1.div.rsrp <- tp1.div.rsrp[!(tp1.div.rsrp$Compartment=="bulk"),]
head(tp1.div.rsrp[,1:10])
dim(tp1.div.rsrp)

dist.tp1.div.rsrp <- vegdist(tp1.div.rsrp[,-(1:5)], method="bray")
set.seed(3)
test48 <- adonis(dist.tp1.div.rsrp~Compartment, data=tp1.div.rsrp, strata=tp1.div.rsrp$Block, permutation=9999)
tp1.div.rsrp.p <- test48$aov.tab$`Pr(>F)`[1]
tp1.div.rsrp.p
#P=0.0001

#TP2 Conv.
tp2.conv.rsrp <- subset(otu.meta, Timepoint=="TP2")
tp2.conv.rsrp <- subset(tp2.conv.rsrp, Soil=="Conv.")
tp2.conv.rsrp <- tp2.conv.rsrp[!(tp2.conv.rsrp$Compartment=="bulk"),]
head(tp2.conv.rsrp[,1:10])
dim(tp2.conv.rsrp)

dist.tp2.conv.rsrp <- vegdist(tp2.conv.rsrp[,-(1:5)], method="bray")
set.seed(3)
test49 <- adonis(dist.tp2.conv.rsrp~Compartment, data=tp2.conv.rsrp, strata=tp2.conv.rsrp$Block, permutation=9999)
tp2.conv.rsrp.p <- test49$aov.tab$`Pr(>F)`[1]
tp2.conv.rsrp.p
#P=0.0005, significant

#TP2 Div.
tp2.div.rsrp <- subset(otu.meta, Timepoint=="TP2")
tp2.div.rsrp <- subset(tp2.div.rsrp, Soil=="Div.")
tp2.div.rsrp <- tp2.div.rsrp[!(tp2.div.rsrp$Compartment=="bulk"),]
head(tp2.div.rsrp[,1:10])
dim(tp2.div.rsrp)

dist.tp2.div.rsrp <- vegdist(tp2.div.rsrp[,-(1:5)], method="bray")
set.seed(3)
test50 <- adonis(dist.tp2.div.rsrp~Compartment, data=tp2.div.rsrp, strata=tp2.div.rsrp$Block, permutation=9999)
tp2.div.rsrp.p <- test50$aov.tab$`Pr(>F)`[1]
tp2.div.rsrp.p
#P=0.0005

#TP3 Conv.
head(otu.meta[,1:10])
tp3.conv.rsrp <- subset(otu.meta, Timepoint=="TP3")
tp3.conv.rsrp <- subset(tp3.conv.rsrp, Soil=="Conv.")
tp3.conv.rsrp <- tp3.conv.rsrp[!(tp3.conv.rsrp$Compartment=="bulk"),]
head(tp3.conv.rsrp[,1:10])
dim(tp3.conv.rsrp)

dist.tp3.conv.rsrp <- vegdist(tp3.conv.rsrp[,-(1:5)], method="bray")
set.seed(3)
test51 <- adonis(dist.tp3.conv.rsrp~Compartment, data=tp3.conv.rsrp, strata=tp3.conv.rsrp$Block, permutation=9999)
tp3.conv.rsrp.p <- test51$aov.tab$`Pr(>F)`[1]
tp3.conv.rsrp.p
#P=0.0004

#TP3 Div.
tp3.div.rsrp <- subset(otu.meta, Timepoint=="TP3")
tp3.div.rsrp <- subset(tp3.div.rsrp, Soil=="Div.")
tp3.div.rsrp <- tp3.div.rsrp[!(tp3.div.rsrp$Compartment=="bulk"),]
head(tp3.div.rsrp[,1:10])
dim(tp3.div.rsrp)

dist.tp3.div.rsrp <- vegdist(tp3.div.rsrp[,-(1:5)], method="bray")
set.seed(3)
test52 <- adonis(dist.tp3.div.rsrp~Compartment, data=tp3.div.rsrp, strata=tp3.div.rsrp$Block, permutation=9999)
tp3.div.rsrp.p <- test52$aov.tab$`Pr(>F)`[1]
tp3.div.rsrp.p
#P=0.0005

#TP4 Conv.
head(otu.meta[,1:10])
tp4.conv.rsrp <- subset(otu.meta, Timepoint=="TP4")
tp4.conv.rsrp <- subset(tp4.conv.rsrp, Soil=="Conv.")
tp4.conv.rsrp <- tp4.conv.rsrp[!(tp4.conv.rsrp$Compartment=="bulk"),]
head(tp4.conv.rsrp[,1:10])
dim(tp4.conv.rsrp)

dist.tp4.conv.rsrp <- vegdist(tp4.conv.rsrp[,-(1:5)], method="bray")
set.seed(3)
test53 <- adonis(dist.tp4.conv.rsrp~Compartment, data=tp4.conv.rsrp, strata=tp4.conv.rsrp$Block, permutation=9999)
tp4.conv.rsrp.p <- test53$aov.tab$`Pr(>F)`[1]
tp4.conv.rsrp.p
#P=0.0005

#TP4 Div.
tp4.div.rsrp <- subset(otu.meta, Timepoint=="TP4")
tp4.div.rsrp <- subset(tp4.div.rsrp, Soil=="Div.")
tp4.div.rsrp <- tp4.div.rsrp[!(tp4.div.rsrp$Compartment=="bulk"),]
head(tp4.div.rsrp[,1:10])
dim(tp4.div.rsrp)

dist.tp4.div.rsrp <- vegdist(tp4.div.rsrp[,-(1:5)], method="bray")
set.seed(3)
test54 <- adonis(dist.tp4.div.rsrp~Compartment, data=tp4.div.rsrp, strata=tp4.div.rsrp$Block, permutation=9999)
tp4.div.rsrp.p <- test54$aov.tab$`Pr(>F)`[1]
tp4.div.rsrp.p
#P=0.001

#Adjust p-values for multiple comparisons
p.adj2 <- p.adjust(c(tp1.conv.rsrp.p, tp1.div.rsrp.p, tp2.conv.rsrp.p, tp2.div.rsrp.p,
                     tp3.conv.rsrp.p, tp3.div.rsrp.p, tp4.conv.rsrp.p, tp4.div.rsrp.p), method="bonferroni", n=54)
p.adj2
#all signficant


##############
#Create graphs

# Figure 3b

# create dataframe from nmds ordination
ord.nmds.df <- as.data.frame(ord.nmds$points)
dim(ord.nmds.df)

# add sample metadata
meta <- as.data.frame(sample_data(physeq.want))
dim(meta)
ord.nmds.meta <- merge(meta, ord.nmds.df, by.x=0, by.y=0)
head(ord.nmds.meta)

#Reorder levels of compartment
levels(ord.nmds.meta$Compartment)
ord.nmds.meta$Compartment <- gsub("bulk", "Bulk soil", ord.nmds.meta$Compartment)
ord.nmds.meta$Compartment <- gsub("rhizosphere", "Rhizosphere soil", ord.nmds.meta$Compartment)
ord.nmds.meta$Compartment <- gsub("rhizoplane", "Rhizoplane", ord.nmds.meta$Compartment)
ord.nmds.meta$Compartment <- factor(ord.nmds.meta$Compartment, levels=c("Bulk soil", "Rhizosphere soil", "Rhizoplane"))

#Rename timepoint to dev stage
levels(ord.nmds.meta$Timepoint)
ord.nmds.meta$Timepoint <- revalue(ord.nmds.meta$Timepoint, c("TP1"="V4", "TP2"="V11", "TP3"="R2", "TP4"="R5")) 

#graph
head(ord.nmds.meta)
tp.nmds.plot <- ggplot(ord.nmds.meta, aes(x=MDS1, y=MDS2, color=Compartment, shape=Soil)) +
  geom_point(size=1) +
  facet_wrap(~Timepoint) +
  stat_ellipse(type="t", size=0.5, aes(linetype=Soil)) +
  labs(title="", x="NMDS1", y="NMDS2") +
  theme_bw() +
  theme(axis.title=element_text(size=8),
        axis.text=element_text(size=8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x = element_text(size = 8, margin=margin(3,3,3,3)),
        legend.position="none",
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  geom_hline(yintercept=0, color="grey69", linetype=2) +
  geom_vline(xintercept=0, color="grey69", linetype=2) +
  guides(fill=FALSE)
tp.nmds.plot

save(tp.nmds.plot, file="betadiv.tp.its.graph.RData")

# Figure 4b

#spell out conv and div completely 
ord.nmds.meta$Soil <- revalue(sample_data(physeq.want)$Soil, c("Conv."="Conventional", "Div."="Diversified"))
levels(sample_data(physeq.want)$Soil)

#graph
prox.nmds.plot <- ggplot(ord.nmds.meta, aes(x=MDS1, y=MDS2, color=Timepoint)) +
  geom_point(size=1) +
  facet_grid(Soil~Compartment) +
  stat_ellipse(type="t", size=0.5) +
  labs(title="", x="NMDS1", y="NMDS2") +
  theme_bw() +
  theme(axis.title=element_text(size=8),
        axis.text=element_text(size=8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x = element_text(size = 8, margin=margin(3,3,3,3)),
        strip.text.y = element_text(size = 8, margin=margin(3,3,3,3)),
        legend.position="none",
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  geom_hline(yintercept=0, color="grey69", linetype=2) +
  geom_vline(xintercept=0, color="grey69", linetype=2) +
  guides(fill=FALSE)
prox.nmds.plot

save(prox.nmds.plot, file="betadiv.prox.its.graph.RData")


################
# Create figures

rm(list=ls())

#Fig 3
load("betadiv.tp.16s.graph.RData")
fig3a <- tp.nmds.plot

load("betadiv.tp.its.graph.RData")
fig3b <- tp.nmds.plot

fig3 <- plot_grid(fig3a, NULL, fig3b,  ncol=3, labels=c("A", "", "B"), label_size=14, rel_widths=c(1,0.05,1))
fig3

ggsave("figure3.tiff", plot=fig3, scale=1, width=7, height=3.5, units="in")

#Fig 4
load("betadiv.prox.16s.graph.RData")
fig4a <- prox.nmds.plot

load("betadiv.prox.its.graph.RData")
fig4b <- prox.nmds.plot

fig4 <- plot_grid(fig4a, NULL, fig4b,  ncol=3, labels=c("A", "", "B"), label_size=14, rel_widths=c(1,0.05,1))
fig4

ggsave("figure4.tiff", plot=fig4, scale=1, width=7, height=3.5, units="in")
