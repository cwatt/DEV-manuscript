#######################
#DEV 16S mock community
#Cassandra Wattenburger

#NOTES: 
#This code explores the mock community data
#Compartment and soil refer to root proximity and system respectively

#Clear workspace and load necessary packages
rm(list=ls())

library("labdsv")

#set working directory to wherever works for you
setwd("Y:/Cassi/Sequencing Data/Argonne 2017/01 DEV Project/02 Github code")


############################
#Prep mock sample count data

#Read in raw count table
raw <- read.csv("DEV16S_otu.csv", row.names = 1)
head(raw[,1:5])
dim(raw)

#Subset only mock samples with OTU labels
mock <- raw[,c("MB1", "MB2", "MB3")]
head(mock)

#Remove 1 and 0 abundance OTUs
mock.trans <- t(mock)
head(mock.trans[,1:5])
mock.nosingle <- dropspc(mock.trans, 1)
dim(mock)
dim(mock.nosingle)
head(mock.nosingle)
tail(mock.nosingle)
#33 OTUs remaning, this is more species than are in the actual mock community
#likely due to uncertainty in taxonomic assignment, contamination, and/or well hopping

#Add in taxonomic info
mock.nosingle.trans <- t(mock.nosingle)
head(mock.nosingle.trans)
#I include column 223 here because I want it to be considered a dataframe for adding row names
tax <- raw[,c(223:224)]
head(tax)
otu.labels <- row.names(raw)
head(otu.labels)
rownames(tax) <- otu.labels
head(tax)
mock.tax <- merge(mock.nosingle.trans, tax, by = 0)
head(mock.tax)
dim(mock.tax)
mock.tax <- mock.tax[,-5]
head(mock.tax)
rownames(mock.tax) <- mock.tax[,1]
mock.tax <- mock.tax[,-1]
head(mock.tax)

#Separate taxonomic levels
mock.tax <- separate(data = mock.tax, 
                     col = taxonomy, 
                     into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                     sep=",")
head(mock.tax)

#Attempt to find the taxa present in the mock communities
#From ZymoBIOMICS Microbial Community DNA Standard Instruction Manual these are:
#1. Pseudomonas aeruginosa, theoretical comp: 4.6
#2. Escherichia coli, theoretical comp: 10
#3. Salmonella enterica, theoretical comp: 11.3
#4. Lactobacillus fermentum, theoretical comp: 18.8
#5. Enterococcus faecalis, theoretical comp: 10.4
#6. Staphylococcus aureus, theoretical comp: 13.3
#7. Listeria monocytogenes, theoretical comp: 15.9
#8. Bacillus subtilis, theoretical comp: 15.7
#And two yeasts (not present)

#Subset Pseudomonas aeruginosa
mock1g <- subset(mock.tax, Genus == "g__Pseudomonas")
mock1g
mock1f <- subset(mock.tax, Family == "f__Pseudomonadaceae")
mock1f
mock1o <- subset(mock.tax, Order == "o__Pseudomonadales")
mock1o
#OTU_155, Pseudomonas auruginosa assigned at order level

#Subset Excherichia coli
mock2g <- subset(mock.tax, Genus=="g__Escherichia")
mock2g
mock2f <- subset(mock.tax, Family == "f__Enterobacteriaceae")
mock2f
#uncertainty in E.coli placement

#Subset Salmonella enterica
mock3g <- subset(mock.tax, Genus=="g__Salmonella")
mock3g
mock3f <- subset(mock.tax, Family == "f__Enterobacteriaceae")
mock3f
#uncertainty

#E.coli and S. enterica are closely related, the two most abundant OTUs in mock2f or mock3f are likely
#them but indistinguishable.
#OTU_119 ... g__Escherichia-Shigella is likely E.coli (Shigella and Escherichia are indistinguishable genetically)
#OTU_9 ... g__Enterobacter is likely S.enterica
mock23f <- subset(mock2f, Genus=="g__Escherichia-Shigella" | Genus=="g__Enterobacter")
mock23f
mock23f <- mock23f[-2,]
mock23f

#Subset Lactobacillus fermentum
mock4g <- subset(mock.tax, Genus=="g__Lactobacillus")
mock4g
#OTU_170, known at genus level

#Subset Enterococcus faecalis
mock5g <- subset(mock.tax, Genus=="g__Enterococcus")
mock5g
#OTU_80, known at genus level

#Subset Staphylococcus aureus
mock6g <- subset(mock.tax, Genus=="g__Staphylococcus")
mock6g
#OTU_39, known at genus level

#Subset Listeria monocytogenes
mock7g <- subset(mock.tax, Genus=="g__Listeria")
mock7g
#OTU_44, known at genus level

#Subset Bacillus subtilis
mock8g <- subset(mock.tax, Genus=="g__Bacillus")
mock8g
#known at genus level, probably OTU_61 (most abundant)
mock8g.otu61 <- mock8g[1,]

#Combine OTUs into single dataframe
mock.otus <- rbind(mock1o, mock23f, mock4g, mock5g, mock6g, mock7g, mock8g.otu61)
head(mock.otus)
dim(mock.otus)
mock.otus.notax <- mock.otus[,-(4:10)]
head(mock.otus.notax)
dim(mock.otus.notax)

#Calculate proportion of total of each OTU
prop <- scale(mock.otus.notax, center = FALSE, 
              scale = colSums(mock.otus.notax))
prop
prop <- prop*100
prop

#Average each OTU across mock samples
mean <- rowMeans(prop)
mean

#Make a nice summary table
theoretical <- c(4.6, 10.0, 11.3, 18.8, 10.4, 13.3, 15.9, 15.7)
theoretical
compare <- rbind(mean, theoretical)
compare
compare.trans <- t(compare)
compare.trans

#The taxonomic level correctly identified at
taxlevel <- c("Order", "Family", "Family", "Family", "Genus", "Genus", "Genus", "Genus")
as.data.frame(taxlevel)

#Species names
species <- c("Pseudomonas aeruginosa", 
             "Escherichia coli", 
             "Salmonella enterica", 
             "Lactobacillus fermentum", 
             "Enterococcus faecalis", 
             "Staphylococcus aureus", 
             "Listeria monocytogenes", 
             "Bacillus subtilis")
species <- as.data.frame(species)
species

summary <- cbind(compare.trans, species, taxlevel)
summary

labels <- c("Actual", "Theoretical", "Species", "Identified at")

colnames(summary) <- labels
summary

