##########################
#DEV 16S CSS normalization 
#Cassandra Wattenburger

#NOTES: 
#This code normalizes the count data to avoid sequencing depth bias
#CSS normalization code was kindly provided by Allison Thompson (allison.thompson@pnnl.gov)

#Clear workspace and load necessary packages
rm(list=ls())

library("vegan")
library("labdsv")
library("ggplot2")

#set working directory to wherever works for you
setwd("Y:/Cassi/Sequencing Data/Argonne 2017/01 DEV Project/02 Github code")


#######################
#Read in raw count data
otu <- read.csv("DEV16S_otu.csv", row.names=1)
str(otu)
head(otu[,1:10])
head(otu[,220:224])

#remove taxonomy column
otu <- otu[,-224]
dim(otu)
head(otu[220:223])

#Remove blanks
#Likely still included due to contamination or well-hopping during Illumina sequencing
blanks <- c("B1.A2", "B1.A1", "B3.A1", "B2.A3", "B3.A2", "B1.A3", "B3.A3")
otu <- otu[,!(names(otu) %in% blanks)]
dim(otu)


######################################################
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

#calculate total abundance (counts) per sample
abund <- rowSums(otu.trans)
head(abund)
min(abund)
max(abund)
hist(abund)

#format for graphing
rich.abund <- rbind(richness, abund)
head(rich.abund)
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
  ggtitle("16S Raw Counts Sample Abundance vs Richness")
#definitely a positive correlation between sampling depth and richness
#nothing falls below hard cut-offs of 10 richness and 1000 abundance
#Mock communities (MB1-3) are lowest richness, as expected


##############
#CSS normalize

#Remove singletons and absent taxa
otu.nosingles <- dropspc(otu.trans, 1)
str(otu.nosingles)
str(otu.trans)
#removed 1781 taxa that only occurred in one or fewer samples

#reformat count data
head(otu.nosingles[,1:10])
otu.nosingles.trans <- t(otu.nosingles)
head(otu.nosingles.trans[,1:10])
otu.nosingles.trans <- data.frame(OTU.ID = row.names(otu.nosingles.trans), otu.nosingles.trans)
head(otu.nosingles.trans[,1:10])

e_data <- otu.nosingles.trans
edata_id <- "OTU.ID"

#create CSS normalization function
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

#normalize
otu_cssnorm <- CSS_Norm(e_data, edata_id)
head(otu_cssnorm[1:10,1:10])
dim(otu_cssnorm)

#graph the richness x abundance of normalized data
#calculate richness per sample
otu_cssnorm <- otu_cssnorm[,-1]
head(otu_cssnorm[,1:5])
otu_cssnorm.trans <- t(otu_cssnorm)
head(otu_cssnorm.trans[,1:10])
richness.norm <- specnumber(otu_cssnorm.trans[,-1])
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
#way less normal...?

#format count data for graphing
rich.abund.norm <- rbind(richness.norm, abund.norm)
head(rich.abund.norm[,1:5])
rich.abund.norm.trans <- t(rich.abund.norm)
head(rich.abund.norm.trans)
rich.abund.norm.trans.data <- data.frame(rich.abund.norm.trans)
head(rich.abund.norm.trans.data)
rich.abund.norm.trans.data$Sample <- row.names(rich.abund.norm.trans.data)
head(rich.abund.norm.trans.data)

#graph richness x abundance
ggplot(rich.abund.norm.trans.data, aes(x=richness.norm, y=abund.norm)) +
  geom_point() +
  geom_vline(xintercept = 10, color="red") +
  geom_hline(yintercept = 1000, color="red") +
  geom_text(aes(label=Sample),hjust=0, vjust=0) +
  ggtitle("16S CSS Normalized Counts Sample Abundance vs Richness")
#woah, remove those mock communities
#likely due to BIG difference in richness between samples and mocks


#############################################
#CSS normalize without mock community samples

#Remove mocks form original data (no blanks included)
head(otu[,1:5])
dim(otu)
mocks <- c("MB1", "MB2", "MB3")
otu.nomock <- otu[,!(names(otu) %in% mocks)]
dim(otu.nomock)
#Remove singletons and OTUs that don't occur
otu.nomock.trans <- t(otu.nomock)
head(otu.nomock.trans[,1:5])
otu.nomocksingles <- dropspc(otu.nomock.trans, 1)
dim(otu.nomocksingles)
#Removed 1783 singleton OTUs

#reformat
head(otu.nomocksingles[,1:5])
otu.nomocksingles.trans <- t(otu.nomocksingles)
head(otu.nomocksingles.trans[,1:5])
otu.nomocksingles.trans <- data.frame(OTU.ID = row.names(otu.nomocksingles.trans), otu.nomocksingles.trans)
head(otu.nomocksingles.trans[,1:5])

#Be careful, I didn't change the function so these object names are the same as for the previous CSS normalization
e_data <- otu.nomocksingles.trans
edata_id <- "OTU.ID"

#Perform normalization
otu_cssnorm.nomock <- CSS_Norm(e_data,edata_id)
head(otu_cssnorm.nomock[1:10,1:10])
dim(otu_cssnorm.nomock)

#graph richness x abundance without mocks

#calculate richness per sample
otu_cssnorm.nomock <- otu_cssnorm.nomock[,-1]
head(otu_cssnorm.nomock[,1:5])
otu_cssnorm.nomock.trans <- t(otu_cssnorm.nomock)
head(otu_cssnorm.nomock.trans[,1:10])
richness.norm.nomock <- specnumber(otu_cssnorm.nomock.trans)
head(richness.norm.nomock)
min(richness.norm.nomock)
max(richness.norm.nomock)
hist(richness.norm.nomock)
#not changed from raw, as expected

#calculate total abundance per sample
abund.norm.nomock <- rowSums(otu_cssnorm.nomock.trans)
head(abund.norm.nomock)
min(abund.norm.nomock)
max(abund.norm.nomock)
hist(abund.norm.nomock)
#not bad

#format data
rich.abund.norm.nomock <- rbind(richness.norm.nomock, abund.norm.nomock)
head(rich.abund.norm.nomock[,1:5])
rich.abund.norm.nomock.trans <- t(rich.abund.norm.nomock)
head(rich.abund.norm.nomock.trans)
rich.abund.norm.nomock.trans.data <- data.frame(rich.abund.norm.nomock.trans)
head(rich.abund.norm.nomock.trans.data)
rich.abund.norm.nomock.trans.data$Sample <- row.names(rich.abund.norm.nomock.trans.data)
head(rich.abund.norm.nomock.trans.data)

#graph
ggplot(rich.abund.norm.nomock.trans.data, aes(x=richness.norm.nomock, y=abund.norm.nomock)) +
  geom_point() +
  geom_vline(xintercept = 10, color="red") +
  geom_hline(yintercept = 1000, color="red") +
  geom_text(aes(label=Sample),hjust=0, vjust=0) +
  ggtitle("16S CSS Normalized Counts without Mocks Sample Abundance vs Richness")
#positive correlation is reduced

#save file with normalized data
write.csv(otu_cssnorm.nomock, "DEV16S_cssnorm.csv")

