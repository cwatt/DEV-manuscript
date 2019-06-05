
# Cassandra Wattenburger
# 01/24/18

# Cleaned up scripts to generate alpha diversity data and figure 2 and Supp Figure S2

############################################

# load necessary packages
library("phyloseq")
library("ggplot2")
library("tidyr")
library("plyr")
library("reshape2")
library("cowplot")
library("vegan")
library("labdsv")

setwd("Y:/Cassi/Sequencing Data/Argonne 2017/01 DEV Project/03 Final Scripts and figures")


##########
# 16S data

rm(list=ls())


###################
#Raw count analysis

#read in 16S raw count data
otu <- read.csv("DEV16S_otu.csv", row.names = 1)
head(otu[,1:5])
dim(otu)

#remove blanks, mock samples, and taxonomy information
head(otu[,220:224])
otu <- otu[,-224]
remove <- c("B1.A2", "B1.A1", "B3.A1", "B2.A3", "B3.A2", "B1.A3", "B3.A3", "MB1", "MB2", "MB3")
otu <- otu[,!(names(otu) %in% remove)]
head(otu[,1:5])
dim(otu)

#Read in metadata
meta <- read.csv("DEV_metadata.csv", row.names = 1)
head(meta)
dim(meta)

#remove sample MDP4 from metadata (sample was dropped from count table during quality filtering)
meta["MDP4",]
otu[,"MDP4"]
meta.filter <- meta[rownames(meta) %in% names(otu),]
str(meta.filter)

#Combine otu and meta
head(meta.filter)
dim(meta.filter)
head(otu[,1:5])
dim(otu)
otu.trans <- t(otu)
dim(otu.trans)
head(otu.trans[,1:5])
otu.meta <- merge(meta.filter, otu.trans, by = 0)
head(otu.meta[,1:10])
dim(otu.meta)
rownames(otu.meta) <- otu.meta[,1]
otu.meta <- otu.meta[,-1]
head(otu.meta[,1:10])

#calculate richness
richness <- specnumber(otu.meta[,-c(1:4)])
head(richness)

#calculate shannon diversity
shannon <- diversity(otu.meta[,-c(1:4)], index="shannon")
head(shannon)

#combine diversity indices into one dataframe
alpha <- cbind(richness, shannon)
head(alpha)

#add metadata
alpha.meta <- cbind(meta.filter, alpha)
head(alpha.meta)

###########
#Statistics

#check normality using histograms
hist(richness)
shapiro.test(richness)
#not normal
hist(shannon)
shapiro.test(shannon)
#left tailed, not normal
hist(log(richness))
hist(log(shannon))
#log transform didn't help
#use Kruskal-Wallis tests because of non-normality

#richness

kr1 <- kruskal.test(richness~Timepoint, data=alpha.meta)
kr1
# time point significant, P = 1.013*10^-11

kr2 <- kruskal.test(richness~Soil, data=alpha.meta)
kr2
# System not significant, P = 0.9247

kr3 <- kruskal.test(richness~Compartment, data=alpha.meta)
kr3
#Root proximity not significant, P = 0.6017

kr4 <- kruskal.test(richness~Block, data=alpha.meta)
kr4
#Block significant, P = 0.002779

#shannon

ks1 <- kruskal.test(shannon~Timepoint, data=alpha.meta)
ks1
#Time point significant, P = 2.962*10^-12

ks2 <- kruskal.test(shannon~Soil, data=alpha.meta)
ks2
#System not significant, P = 0.4837

ks3 <- kruskal.test(shannon~Compartment, data=alpha.meta)
ks3
#Root proximity not significant, P = 0.9872

ks4 <- kruskal.test(shannon~Block, data=alpha.meta)
ks4
#Block significant, P = 0.02251

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
#files created in "________________________________________________________________________________________.R"

#read in 16S normalized count data
otu.norm <- read.csv("DEV16S_phyloseq_otu.csv", row.names = 1)
head(otu.norm[,1:5])
#metadata
meta.norm <- read.csv("DEV16S_phyloseq_metadata.csv", row.names = 2)
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

#combine diversity metrics into one dataframe
alpha.norm <- cbind(richness.norm, shannon.norm)
head(alpha.norm)
#add meta.normdata
alpha.meta.norm <- cbind(meta.norm, alpha.norm)
head(alpha.meta.norm)

###########
#Statistics

#check normality using histograms

hist(richness.norm)
shapiro.test(richness.norm)
#not normal

hist(shannon.norm)
shapiro.test(shannon.norm)
#left tailed, not normal

hist(log(richness.norm))
hist(log(shannon.norm))
#log transform didn't help
#use Kruskal-Wallis tests because of non-normality

#richness

kr1.norm <- kruskal.test(richness.norm~Timepoint, data=alpha.meta.norm)
kr1.norm
#Time point significant, P = 5.406*10^-12

kr2.norm <- kruskal.test(richness.norm~Soil, data=alpha.meta.norm)
kr2.norm
#System not significant, P = 0.8772

kr3.norm <- kruskal.test(richness.norm~Compartment, data=alpha.meta.norm)
kr3.norm
#Root proximity not significant, P = 0.5634

kr4.norm <- kruskal.test(richness.norm~Block, data=alpha.meta.norm)
kr4.norm
#Block significant, P = 0.002423

#significance did not change compared to raw data

#shannon

ks1.norm <- kruskal.test(shannon.norm~Timepoint, data=alpha.meta.norm)
ks1.norm
#Time point significant, P = 2.034*10^-10

ks2.norm <- kruskal.test(shannon.norm~Soil, data=alpha.meta.norm)
ks2.norm
#System not significant, P = 0.5019

ks3.norm <- kruskal.test(shannon.norm~Compartment, data=alpha.meta.norm)
ks3.norm
#Root proximity not significant, P = 0.8478

ks4.norm <- kruskal.test(shannon.norm~Block, data=alpha.meta.norm)
ks4.norm
#Block significant, P = 0.02523
#significance did not change using pre-processed data

######################
# Create graphing data

# Figure 2A

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

# add metadata
head(rich.time)
rich.time$index <- rep("Richness (no. OTUs)", nrow(rich.time))

head(shannon.time)
shannon.time$index <- rep("Shannon's diversity index", nrow(shannon.time))

fig2.data.16s <- rbind(rich.time, shannon.time)

save(fig2.data.16s, file="fig2.data.16s.RData")


# Supp Fig 2A

head(alpha.meta.norm)

#  summarize data
alpha.meta.norm <- alpha.meta.norm[,-4] #remove block variable
rich <- ddply(alpha.meta.norm, .(Timepoint, Soil, Compartment, Treatment), summarize, .progress="text", mean=mean(richness.norm), se=(sd(richness.norm)/sqrt(length(richness.norm))))
head(rich)

# rename timepoints to developmental stage
levels(rich$Timepoint)
rich$Timepoint <- revalue(rich$Timepoint, c("TP1"="V4", "TP2"="V11", "TP3"="R2", "TP4"="R5")) 
levels(rich$Timepoint)
head(rich)

# add index metadata
rich$index <- rep("Richness (no. OTUs)", nrow(rich))
head(rich)

# summarize Shannon's Diversity data
shan <- ddply(alpha.meta.norm, .(Timepoint, Soil, Compartment, Treatment), summarize, .progress="text", mean=mean(shannon.norm), se=(sd(shannon.norm)/sqrt(length(shannon.norm))))
head(shan)

# rename timepoints to developmental stage
levels(shan$Timepoint)
shan$Timepoint <- revalue(shan$Timepoint, c("TP1"="V4", "TP2"="V11", "TP3"="R2", "TP4"="R5")) 
levels(shan$Timepoint)

# add index metadata
shan$index <- rep("Shannon's diversity index", nrow(shan))
head(shan)

# merge dataframe
suppfig2.data.16s <- rbind(rich, shan)
head(suppfig2.data.16s)

save(suppfig2.data.16s, file="suppfig2.data.16s.RData")


##########
# ITS Data

rm(list=ls())


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

#combine diversity metrics into one dataframe
alpha <- cbind(richness, shannon)
head(alpha)

#add metadata
alpha.meta <- cbind(meta.filter, alpha)
head(alpha.meta)

###########
#Statistics

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

#combine diversity metrics into one dataframe
alpha.norm <- cbind(richness.norm, shannon.norm)
head(alpha.norm)

#add meta.norm data
alpha.meta.norm <- cbind(meta.norm, alpha.norm)
head(alpha.meta.norm)

###########
#Statistics

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


######################
# Create graphing data

# Fig 2B

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

# add index metadata
head(rich.time)
rich.time$index <- rep("Richness (no. OTUs)", nrow(rich.time))

head(shannon.time)
shannon.time$index <- rep("Shannon's diversity index", nrow(shannon.time))

fig2.data.its <- rbind(rich.time, shannon.time)

save(fig2.data.its, file="fig2.data.its.RData")


#Supp Fig 2B

head(alpha.meta.norm)

#Calculate richness data
alpha.meta.norm <- alpha.meta.norm[,-4] #remove block variable
rich <- ddply(alpha.meta.norm, .(Timepoint, Soil, Compartment, Treatment), summarize, .progress="text", mean=mean(richness.norm), se=(sd(richness.norm)/sqrt(length(richness.norm))))
head(rich)

#rename timepoints to developmental stage
levels(rich$Timepoint)
rich$Timepoint <- revalue(rich$Timepoint, c("TP1"="V4", "TP2"="V11", "TP3"="R2", "TP4"="R5")) 
levels(rich$Timepoint)

# add index metadata
rich$index <- rep("Richness (no. OTUs)", nrow(rich))
head(rich)

#Calculate Shannon's Diversity data
shan <- ddply(alpha.meta.norm, .(Timepoint, Soil, Compartment, Treatment), summarize, .progress="text", mean=mean(shannon.norm), se=(sd(shannon.norm)/sqrt(length(shannon.norm))))
head(shan)

#rename timepoints to developmental stage
levels(shan$Timepoint)
shan$Timepoint <- revalue(shan$Timepoint, c("TP1"="V4", "TP2"="V11", "TP3"="R2", "TP4"="R5")) 
levels(shan$Timepoint)

# add index metadata
shan$index <- rep("Shannon's diversity index", nrow(shan))
head(shan)

# merge dataframes
suppfig2.data.its <- rbind(rich, shan)

save(suppfig2.data.its, file="suppfig2.data.its.RData")


################
# Create Figures

# Create figure 2

rm(list=ls())

load("fig2.data.16s.RData")
load("fig2.data.its.RData")

fig2a <- ggplot(fig2.data.16s, aes(y=mean, x=Days)) +
  geom_pointrange(aes(x=Days, y=mean, ymin=(mean-se), ymax=(mean+se)), size=0.2) +
  facet_wrap(~index, ncol=1, scales=("free_y"), strip.position="left") +
  scale_x_continuous(limits=c(0, 130)) +
  xlab("Days after planting") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.key = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour="grey20",size=8),
        axis.text.y = element_text(colour="grey20",size=8),
        axis.title.x = element_text(size=8),
        strip.text.y = element_text(size = 8, margin=margin(3,3,3,3)),
        plot.margin = unit(c(0.5, 0.25, 0, 1), "cm"))
fig2a

fig2b <- ggplot(fig2.data.its, aes(y=mean, x=Days)) +
  geom_pointrange(aes(x=Days, y=mean, ymin=(mean-se), ymax=(mean+se)), size=0.2) +
  facet_wrap(~index, ncol=1, scales=("free_y"), strip.position="left") +
  scale_x_continuous(limits=c(0, 130)) +
  xlab("Days after planting") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.key = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour="grey20",size=8),
        axis.text.y = element_text(colour="grey20",size=8),
        axis.title.x = element_text(size=8),
        strip.text.y = element_text(size = 8, margin=margin(3,3,3,3)),
        plot.margin = unit(c(0.5, 0.25, 0, 1), "cm"))
fig2b

fig2 <- plot_grid(fig2a, NULL, fig2b, nrow=1, labels=c("A", "","B"), label_size=14, rel_widths=c(1,0.01,1))
fig2

ggsave("figure2.tiff", plot=fig2, scale=1, width=5, height=3.5, units="in")


# Create Supp figure S2

load("suppfig2.data.16s.RData")
load("suppfig2.data.its.RData")

suppfig2a <- ggplot(suppfig2.data.16s, aes(y=mean, x=Soil)) +
  geom_pointrange(aes(x=Soil, y=mean, ymin=(mean-se), ymax=(mean+se), shape=Compartment), size=0.2) +
  facet_grid(index~Timepoint, scales="free_y") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_rect(colour="black", size=0.5),
        legend.title = element_blank(),
        legend.position = "none",
        strip.text = element_text(size=8),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, colour="grey20",size=8),
        axis.text.y = element_text(colour="grey20",size=8),
        strip.text.y = element_text(size = 8, margin=margin(3,3,3,3)),
        strip.text.x = element_text(size = 8, margin=margin(3,3,3,3)))
suppfig2a


suppfig2b <- ggplot(suppfig2.data.its, aes(y=mean, x=Soil)) +
  geom_pointrange(aes(x=Soil, y=mean, ymin=(mean-se), ymax=(mean+se), shape=Compartment), size=0.2) +
  facet_grid(index~Timepoint, scales="free_y") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_rect(colour="black", size=0.5),
        legend.title = element_blank(),
        legend.position = "none",
        strip.text = element_text(size=8),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, colour="grey20",size=8),
        axis.text.y = element_text(colour="grey20",size=8),
        strip.text.y = element_text(size = 8, margin=margin(3,3,3,3)),
        strip.text.x = element_text(size = 8, margin=margin(3,3,3,3)))
suppfig2b

suppfig2 <- plot_grid(suppfig2a, suppfig2b, nrow=1, labels=c("A","B"), label_size=14)
suppfig2

ggsave("suppfig2.tiff", plot=suppfig2, width=7, height=4, units="in")


