###########################
#DEV ITS CSS normalization
#Cassandra Wattenburger

#NOTES: 
#This code normalizes the count data to avoid sequencing depth bias
#CSS normalization code was kindly provided by Allison Thompson (allison.thompson@pnnl.gov

#Clear workspace and load necessary packages
rm(list=ls())

#load in necessary packages
library("vegan")
library("labdsv")
library("ggplot2")

#set working directory to wherever works for you
setwd("Y:/Cassi/Sequencing Data/Argonne 2017/01 DEV Project/02 Github code")


#######################
#Read in raw count data
otu <- read.csv("DEVITS_otu.csv", row.names=1)
str(otu)
head(otu[,1:10])
head(otu[,210:216])

#remove taxonomy column
otu <- otu[,-216]
str(otu)
head(otu[210:215])

#Remove blanks
#Likely still included due to contamination or well-hoping during Illumina sequencing
blanks <- c("B1.A02", "B3.A01", "B3.A02", "B3.A03")
otu <- otu[,!(names(otu) %in% blanks)]
str(otu)

###################################
#Quality check data for bad samples
#abundance x richness graph for each sample

#transpose data frame
otu.trans <- t(otu)
head(otu.trans[,1:10])

#calculate richness per sample
richness <- specnumber(otu.trans)
head(richness)
min(richness)
max(richness)
hist(richness)

#calculate total abundance per sample
abund <- rowSums(otu.trans)
head(abund)
min(abund)
max(abund)
hist(abund)

#format for graphing
rich.abund <- rbind(richness, abund)
head(rich.abund[,1:5])
rich.abund.trans <- t(rich.abund)
head(rich.abund.trans)
str(rich.abund.trans)
rich.abund.trans.data <- data.frame(rich.abund.trans)
rich.abund.trans.data$Sample <- row.names(rich.abund.trans.data)
str(rich.abund.trans.data)

#graph
ggplot(rich.abund.trans.data, aes(x=richness, y=abund)) +
  geom_point() +
  geom_vline(xintercept = 10, color="red") +
  geom_hline(yintercept = 1000, color="red") +
  geom_text(aes(label=Sample),hjust=0, vjust=0) +
  ggtitle("ITS Raw Counts Sample Abundance vs Richness")
#MDP71, MDP66, MDP25 fall below hard cut-off of 10 richness and/or 1000 abundance
#Mock communities (MF1-3) have low richness, as expected

#remove MDP71, MDP66, MDP25
remove <- c("MDP71", "MDP66", "MDP25")
otu.filter <- otu[,!names(otu) %in% remove]
str(otu.filter)
head(otu.filter[,1:5])


##############
#CSS normalize

#remove mock communities before CSS normalization
#see "DEV16S_css_normalization_clean.R" to see how mock communities weirdly affect normalization
mock <- c("MF1", "MF2", "MF3")
otu.filter <- otu.filter[,!names(otu.filter) %in% mock]
str(otu.filter)

#Remove singletons and absent taxa
otu.filter.trans <- t(otu.filter)
head(otu.filter.trans[,1:5])
otu.nosingles <- dropspc(otu.filter.trans, 1)
str(otu.nosingles)
str(otu.filter.trans)
#removed 445 taxa that only occurred in one or fewer samples

#reformat
head(otu.nosingles[,1:10])
otu.nosingles.trans <- t(otu.nosingles)
head(otu.nosingles.trans[,1:10])
otu.nosingles.trans <- data.frame(OTU.ID = row.names(otu.nosingles.trans), otu.nosingles.trans)
head(otu.nosingles.trans[,1:10])

e_data <- otu.nosingles.trans
edata_id <- "OTU.ID"

#CSS normalization function
#code from Allison Thompson (allison.thompson@pnnl.gov)
CSS_Norm <- function(e_data, edata_id, q=0.75, qg="median"){
  e_data_norm <- e_data[,-which(colnames(e_data)==edata_id)]
  e_data_norm[e_data_norm==0] <- NA
  
  # calculate the q quantile of data, per sample and globally #
  col.q = apply(e_data_norm, 2, function(x) sum(x[x<=quantile(x, probs = q, na.rm=TRUE)], na.rm=TRUE))
  
  if(qg == 1000){
    g.q = 1000
  }else if(qg == "median"){
    g.q = median(col.q, na.rm=TRUE)
  }else{
    stop("Invalid value for qg")
  }
  #g.q = 1000
  #g.q = median(col.q, na.rm=TRUE)
  
  # normalize omics_data data by q quantile and transform back to count data #
  for(i in 1:ncol(e_data_norm)){
    e_data_norm[,i] = (e_data_norm[,i] / col.q[i]) * g.q
  }
  
  e_data <- data.frame(e_data[,which(colnames(e_data)==edata_id)],e_data_norm)
  colnames(e_data)[1] <- edata_id
  e_data[is.na(e_data)] <- 0
  return(e_data)
}

#Normalize the data
otu_cssnorm <- CSS_Norm(e_data,edata_id)
head(otu_cssnorm[,1:5])
dim(otu_cssnorm)

#calculate richness per sample
otu_cssnorm <- otu_cssnorm[,-1]
head(otu_cssnorm[,1:5])
otu_cssnorm.trans <- t(otu_cssnorm)
head(otu_cssnorm.trans[,1:5])
richness.norm <- specnumber(otu_cssnorm.trans)
head(richness.norm)
min(richness.norm)
max(richness.norm)
hist(richness.norm)
#not changed from raw, as expected

#calculate total abundance per sample
abund.norm <- rowSums(otu_cssnorm.trans)
head(abund.norm)
min(abund.norm)
max(abund.norm)
hist(abund.norm)
#very right tailed...

#make dataframe for graphing
rich.abund.norm <- rbind(richness.norm, abund.norm)
head(rich.abund.norm[,1:5])
rich.abund.norm.trans <- t(rich.abund.norm)
head(rich.abund.norm.trans)
rich.abund.norm.trans.data <- data.frame(rich.abund.norm.trans)
head(rich.abund.norm.trans.data)
rich.abund.norm.trans.data$Sample <- row.names(rich.abund.norm.trans.data)
head(rich.abund.norm.trans.data)

#graph
ggplot(rich.abund.norm.trans.data, aes(x=richness.norm, y=abund.norm)) +
  geom_point() +
  geom_vline(xintercept = 10, color="red") +
  geom_hline(yintercept = 1000, color="red") +
  geom_text(aes(label=Sample),hjust=0, vjust=0) +
  ggtitle("ITS CSS Normalized Counts Sample Abundance vs Richness")
#seems to have removed some positive correlation skew

#save file with normalized data
write.csv(otu_cssnorm, "DEVITS_cssnorm.csv")
