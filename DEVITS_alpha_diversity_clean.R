########################
#DEV ITS Alpha Diversity 
#Cassandra Wattenburger

#NOTES: 
#This code analyzes the 16S count data alpha diversity
#Uses data generated in "DEVITS_phyloseq_prep_clean.R"
#Compartment and soil refer to root proximity and system respectively

#Clear workspace and load necessary packages
rm(list=ls())

library("vegan")
library("labdsv")
library("ggplot2")
library("reshape2")
library("plyr")

#set working directory to wherever works for you
setwd("Y:/Cassi/Sequencing Data/Argonne 2017/01 DEV Project/02 Github code")


###################
#Raw count analysis

#read in raw count data
otu <- read.csv("DEVITS_otu.csv", row.names = 1)
head(otu[,1:5])
dim(otu)
head(otu[,214:216])

#Remove blanks, mock samples, taxonomy
otu <- otu[,-216]
dim(otu)
head(otu[,214:215])
remove <- c("B1.A02", "B3.A01", "B3.A02", "B3.A03", "MF1", "MF2", "MF3")
otu <- otu[,!(names(otu) %in% remove)]
head(otu[,1:5])
dim(otu)

otu.trans <- t(otu)
head(otu.trans[,1:5])

#Read in metadata
meta <- read.csv("DEV_metadata.csv", row.names = 1)
head(meta)
dim(meta)

#Remove samples from metadata that were dropped via quality filtering after sequencing
meta.filter <- meta[rownames(meta) %in% names(otu),]
dim(meta.filter)
dim(otu)

#Merge normalized count data and metadata
head(otu.trans[,1:5])
head(meta.filter)
otu.meta <- merge(meta.filter, otu.trans, by= 0)
head(otu.meta[,1:10])
rownames(otu.meta) <- otu.meta[,1]
otu.meta <- otu.meta[,-1]
head(otu.meta[,1:10])

#calculate richness
richness <- specnumber(otu.meta[,-c(1:5)])
head(richness)

#calculate shannon diversity
shannon <- diversity(otu.meta[,-c(1:5)], index="shannon")
head(shannon)

#check normality using histograms
hist(richness)
shapiro.test(richness)
#not normal
hist(shannon)
shapiro.test(shannon)
#left tailed, not normal
hist(log(richness))
hist(log(shannon))
#log transformdidn't help
#use Kruskal-Wallis tests because of non-normality

#combine diversity metrics into one dataframe
alpha <- cbind(richness, shannon)
head(alpha)

#add metadata
alpha.meta <- cbind(meta.filter, alpha)
head(alpha.meta)

#Statistics

#richness
kr1 <- kruskal.test(richness~Timepoint, data=alpha.meta)
kr1
#Time point not significant, P = 0.6553
kr2 <- kruskal.test(richness~Soil, data=alpha.meta)
kr2
#System not significant, P = 0.783
kr3 <- kruskal.test(richness~Compartment, data=alpha.meta)
kr3
#Root proximity not significant, P = 0.3283
kr4 <- kruskal.test(richness~Block, data=alpha.meta)
kr4
#Block significant, P = 6.826*10^-6

#shannon
ks1 <- kruskal.test(shannon~Timepoint, data=alpha.meta)
ks1
#Time point significant, P = 4.914*10^-7
ks2 <- kruskal.test(shannon~Soil, data=alpha.meta)
ks2
#System not significant, P = 0.3634
ks3 <- kruskal.test(shannon~Compartment, data=alpha.meta)
ks3
#Root proximity not significant, P = 0.1875
ks4 <- kruskal.test(shannon~Block, data=alpha.meta)
ks4
#Block significant, P = 9.9097*10^-5

#calculate mean and standard error
alpha.meta1 <- alpha.meta[,-4]
head(alpha.meta1)
alpha.melt <- melt(alpha.meta1, id.vars=c("Timepoint","Soil","Compartment"))
head(alpha.melt)
summary.rich <- ddply(alpha.melt, c("Timepoint", "Soil", "Compartment","variable"), summarize, mean=mean(value), sd=sd(value), se=sd(value)/sqrt(length(value)))
summary.rich


###############################
#Normalized count data analysis
#this data has been CSS normalized and preprocessed to remove unwanted sequences
#files created in "DEV16S_phyloseq_prep_clean.R"

#read in ITS normalized count data
otu.norm <- read.csv("DEVITS_phyloseq_otu.csv", row.names = 1)
head(otu.norm[,1:5])

#metadata
meta.norm <- read.csv("DEVITS_phyloseq_metadata.csv", row.names = 2)
meta.norm <- meta.norm[,-1]
head(meta.norm)
dim(meta.norm)

#Combine normalized otu counts and meta.normdata
otu.norm.trans <- t(otu.norm)
head(otu.norm.trans[,1:5])
head(meta.norm)
dim(meta.norm)
dim(otu.norm.trans)
otu.meta.norm <- merge(meta.norm, otu.norm.trans, by = 0)
head(otu.meta.norm[,1:10])
dim(otu.meta.norm)
rownames(otu.meta.norm) <- otu.meta.norm$Row.names
otu.meta.norm <- otu.meta.norm[,-1]
head(otu.meta.norm[,1:10])

#calculate richness
richness.norm <- specnumber(otu.meta.norm[,-c(1:5)])
head(richness.norm)

#calculate shannon diversity
shannon.norm <- diversity(otu.meta.norm[,-c(1:5)], index="shannon")
head(shannon.norm)

#check normality using histograms
hist(richness.norm)
shapiro.test(richness.norm)
#looks good
hist(shannon.norm)
shapiro.test(shannon.norm)
#left tailed, not normal
hist(log(richness.norm))
hist(log(shannon.norm))
#log transform didn't help
#use Kruskal-Wallis tests because of non-normality

#combine diversity metrics into one dataframe
alpha.norm <- cbind(richness.norm, shannon.norm)
head(alpha.norm)

#add meta.norm data
alpha.meta.norm <- cbind(meta.norm, alpha.norm)
head(alpha.meta.norm)

#Statistics

#richness
kr1.norm <- kruskal.test(richness.norm~Timepoint, data=alpha.meta.norm)
kr1.norm
#Time point not significant, P = 0.3178
kr2.norm <- kruskal.test(richness.norm~Soil, data=alpha.meta.norm)
kr2.norm
#System not significant, P = 0.8138
kr3.norm <- kruskal.test(richness.norm~Compartment, data=alpha.meta.norm)
kr3.norm
#Root proximity not significant, P = 0.08823
kr4.norm <- kruskal.test(richness.norm~Block, data=alpha.meta.norm)
kr4.norm
#Block significant, P = 1.023*10^-6
#significance did not change using pre-processed data

#shannon
ks1.norm <- kruskal.test(shannon.norm~Timepoint, data=alpha.meta.norm)
ks1.norm
#Time point significant, P = 8.223*10^-11
ks2.norm <- kruskal.test(shannon.norm~Soil, data=alpha.meta.norm)
ks2.norm
#System not significant, P = 0.394
ks3.norm <- kruskal.test(shannon.norm~Compartment, data=alpha.meta.norm)
ks3.norm
#Root proximity not significant, P = 0.4544
ks4.norm <- kruskal.test(shannon.norm~Block, data=alpha.meta.norm)
ks4.norm
#Block significant, P = 0.0003386
#significance did not change using pre-processed data


########################
#Graph data for Figure 2

#add days after planting to metadata for graphing purposes
alpha.meta.norm <- mutate(alpha.meta.norm, Days=alpha.meta.norm$Timepoint)
head(alpha.meta.norm)
alpha.meta.norm$Days <- ifelse(alpha.meta.norm$Days=="TP1", paste0("35"), paste0(alpha.meta.norm$Days))
alpha.meta.norm$Days <- ifelse(alpha.meta.norm$Days=="TP2", paste0("55"), paste0(alpha.meta.norm$Days))
alpha.meta.norm$Days <- ifelse(alpha.meta.norm$Days=="TP3", paste0("82"), paste0(alpha.meta.norm$Days))
alpha.meta.norm$Days <- ifelse(alpha.meta.norm$Days=="TP4", paste0("119"), paste0(alpha.meta.norm$Days))
head(alpha.meta.norm)
tail(alpha.meta.norm)
alpha.meta.norm$Days <- as.numeric(alpha.meta.norm$Days)

#summarize data over system and root proximity because time point is the only significant treatment effect
#richness
rich.time <- ddply(alpha.meta.norm, .(Timepoint, Days), summarize, .progress="text", mean=mean(richness.norm), se=(sd(richness.norm)/sqrt(length(richness.norm))))
head(rich.time)
#shannon
shannon.time <- ddply(alpha.meta.norm, .(Timepoint, Days), summarize, .progress="text", mean=mean(shannon.norm), se=(sd(shannon.norm)/sqrt(length(shannon.norm))))
head(shannon.time)

#graph richness
rich.graph <- ggplot(rich.time) + 
  geom_pointrange(aes(x=Days, y=mean, ymin=(mean-se), ymax=(mean+se))) +
  ggtitle("") +
  ylab("Richness (No. Species)") +
  xlab("Days after emergence") +
  scale_x_continuous(limits=c(0, 130)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.key = element_blank(),
        plot.title = element_text(size=20),
        axis.title = element_text(size=12),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(colour="grey20",size=10),
        axis.text.y = element_text(colour="grey20",size=10))
rich.graph

#graph shannon
shannon.graph <- ggplot(shannon.time) + 
  geom_pointrange(aes(x=Days, y=mean, ymin=(mean-se), ymax=(mean+se))) +
  ggtitle("") +
  ylab("Shannon's Diversity Index") +
  xlab("Days after emergence") +
  scale_x_continuous(limits=c(0, 130)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.key = element_blank(),
        plot.title = element_text(size=20),
        axis.title = element_text(size=12),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(colour="grey20",size=10),
        axis.text.y = element_text(colour="grey20",size=10))
shannon.graph


#############################
#Graph Supplementary Figure 1

head(alpha.meta.norm)

#Calculate richness data
alpha.meta.norm <- alpha.meta.norm[,-4] #remove block variable
rich <- ddply(alpha.meta.norm, .(Timepoint, Soil, Compartment, Treatment), summarize, .progress="text", mean=mean(richness.norm), se=(sd(richness.norm)/sqrt(length(richness.norm))))
head(rich)

#Change level order of compartment factor
with(rich, levels(Compartment))
rich <- within(rich, Compartment <- factor(Compartment, levels=c("bulk","rhizosphere","rhizoplane")))
with(rich, levels(Compartment))

#graph richness
rich.supp <- ggplot(rich) + 
  geom_pointrange(aes(x=Soil, y=mean, ymin=(mean-se), ymax=(mean+se), colour=Compartment)) +
  facet_grid(. ~ Timepoint) +
  ylab("Richness (number of species)") +
  xlab("") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_rect(colour="black", size=0.5),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(size=20),
        axis.title = element_text(size=10),
        axis.text.x = element_text(angle = 45, hjust = 1, colour="grey20",size=10),
        axis.text.y = element_text(colour="grey20",size=10))
rich.supp

#Calculate Shannon's Diversity data
shan <- ddply(alpha.meta.norm, .(Timepoint, Soil, Compartment, Treatment), summarize, .progress="text", mean=mean(shannon.norm), se=(sd(shannon.norm)/sqrt(length(shannon.norm))))
head(shan)
#change level order of compartment factor
with(shan, levels(Compartment))
shan <- within(shan, Compartment <- factor(Compartment, levels=c("bulk","rhizosphere","rhizoplane")))
with(shan, levels(Compartment))

#graph shannon
shan.supp <- ggplot(shan) + 
  geom_pointrange(aes(x=Soil, y=mean, ymin=(mean-se), ymax=(mean+se), colour=Compartment)) +
  facet_grid(. ~ Timepoint) +
  ylab("Shannon's Diversity Index") +
  xlab("") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_rect(colour="black", size=0.5),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(size=20),
        axis.title = element_text(size=10),
        axis.text.x = element_text(angle = 45, hjust = 1, colour="grey20",size=10),
        axis.text.y = element_text(colour="grey20",size=10))
shan.supp

