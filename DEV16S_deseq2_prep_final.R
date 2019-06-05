####################
#DEV 16S DESeq2 prep
#Cassandra Wattenburger

#NOTES: 
#This code creates DESeq2 ready data to analyze differential abundance
#DESeq2 requires raw, un-normalized count data
#Compartment and soil refer to root proximity and system respectively

#Clear workspace and load necessary packages
rm(list=ls())

library("phyloseq")
library("DESeq2")
library("labdsv")
library("tidyr")

#set working directory to wherever works for you
setwd("Y:/Cassi/Sequencing Data/Argonne 2017/01 DEV Project/02 Github code")


######################
#Format raw count data

#Read in raw otu data
raw <- read.csv("DEV16S_otu.csv")
head(raw[,1:5])
dim(raw)
raw[,225]

#taxonomy info
tax <- raw[,c(1,225)]
dim(tax)
head(tax)

#Remove taxonomy from otu count data
raw <- raw[,-225]
dim(raw)
head(raw[,1:5])
rownames(raw) <- raw[,1]
head(raw[,1:5])
raw <- raw[,-1]
head(raw[,1:5])
dim(raw)

#Remove mocks and blank samples from raw
remove <- c("MB1", "MB2", "MB3", "B1.A2", "B1.A1", "B3.A1", "B2.A3", "B3.A2", "B1.A3", "B3.A3")
raw.filter <- raw[,!(names(raw) %in% remove)]
dim(raw.filter)

#Remove singletons and OTUs that don't occur
raw.trans <- t(raw.filter)
head(raw.trans[,1:5])
dim(raw.trans)
raw.nosingles <- dropspc(raw.trans, 1)
dim(raw.nosingles)
#1783 singleton OTUs removed

#Remove singleton OTUs from tax as well
head(tax)
head(raw.nosingles[,1:5])
raw.nosingles.trans <- t(raw.nosingles)
head(raw.nosingles.trans[,1:5])
tax.nosingles <- tax[(tax$X.OTU.ID %in% rownames(raw.nosingles.trans)),]
head(tax.nosingles)
dim(tax.nosingles)
rownames(tax.nosingles) <- tax.nosingles$X.OTU.ID
head(tax.nosingles)

#Merge tax and raw.nosingles.trans to put OTUs in same order
head(raw.nosingles.trans[,1:5])
nosingles.tax <- merge(raw.nosingles.trans, tax.nosingles, by = 0)
head(nosingles.tax[,1:5])
dim(nosingles.tax)
head(nosingles.tax[,215:216])
nosingles.tax <- nosingles.tax[,-215]
dim(nosingles.tax)

#Separate taxonomy and count data
tax.split <- nosingles.tax[,c(1,215)]
head(tax.split)
raw.split <- nosingles.tax[,c(1:214)]
head(raw.split[,1:5])
dim(raw.split)

#Create phyloseq and DESeq2 ready raw count data file
rownames(raw.split) <- raw.split$Row.names
head(raw.split[,1:5])
raw.split <- raw.split[,-1]
head(raw.split[,1:5])
dim(raw.split)

#save file
write.csv(raw.split, "DEV16S_deseq2_otu.csv")

#Create phyloseq and DESeq2 ready taxonomy file
#separate taxonomic levels into separate columns
head(tax.split)
tax.levels <- separate(data = tax.split, 
                       col = taxonomy, 
                       into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                       sep=",")
head(tax.levels)
rownames(tax.levels) <- tax.levels$Row.names
head(tax.levels)
tax.levels <- tax.levels[,-1]
head(tax.levels)

#save file
write.csv(tax.levels, "DEV16S_deseq2_taxonomy.csv")


#######################
#Create phyloseq object

#Raw OTU count data, file: DEV16S_deseq2_otu.csv
otu <- read.csv("DEV16S_deseq2_otu.csv", row.names = 1)
head(otu[,1:5])
class(otu)
otu.m <- as.matrix(otu)
class(otu.m)

OTU = otu_table(otu.m, taxa_are_rows=TRUE)

#Taxonomy data, file: DEV16S_deseq2_taxonomy.csv
tax <- read.csv("DEV16S_deseq2_taxonomy.csv", row.names = 1)
head(tax)
class(tax)
tax.m <- as.matrix(tax)
class(tax.m)

TAX = tax_table(tax.m)

#Metdata, file: DEV16S_phyloseq_metadata.csv
#Metadata will be the same for raw and normalized count data
meta <- read.csv("DEV16S_phyloseq_metadata.csv", row.names = 2)
meta <- meta[,-1]
head(meta)
#Need to carefully set factor levels so we understand differentual abundances
with(meta, levels(Soil))
#soil looks good, Conv. soil will be baseline for comparisons
with(meta, levels(Compartment))
#change so that level order is bulk, rhizosphere, rhizoplane
meta <- within(meta, Compartment <- factor(Compartment, levels=c("bulk","rhizosphere","rhizoplane")))
with(meta, levels(Compartment))
#now bulk will be baseline for comparison unless otherwise specified
with(meta, levels(Timepoint))
#timepoint levels make sense, TP1 is first and they are in chronological order
#Group the variables so that you can do contrasts later
#group variables
meta$group <- factor(paste0(meta$Timepoint, meta$Soil, meta$Compartment))
with(meta, levels(group))
#Create treatment variable without block for easier labelling
meta$Treatment <- with(meta, paste0(Timepoint, sep=" ", Soil, sep=" ", Compartment))
head(meta)

SAM = sample_data(meta, errorIfNULL=TRUE)

#Create phyloseq object from components for DESeq2
deseq = phyloseq(OTU, TAX, SAM)
deseq


###############
#Pre-processing
#Same process as done for alpha and beta diversity analyses

#Remove taxa with unknown phylum-level taxonomy
badTaxa = subset_taxa(deseq, Phylum=="p__?")
badTaxa
#967 unknown sequences were found
#find chloroplast sequences
badTaxa1 = subset_taxa(deseq, Class=="c__Chloroplast")
badTaxa1
#66 chloroplast sequences need to be removed
#find mitochondrial sequences
badTaxa2 = subset_taxa(deseq, Family=="f__Mitochondria")
badTaxa2
#86 mitochondrial sequences need to be removed
#merge all unwanted data together
badTaxa.all = merge_phyloseq(badTaxa, badTaxa1, badTaxa2)
badTaxa.all
#1119 total taxa need to be removed before analysis

#remove unwanted taxa
removeTaxa = taxa_names(badTaxa.all)
allTaxa = taxa_names(deseq)
keepTaxa = allTaxa[!(allTaxa %in% removeTaxa)]
keepTaxa
#keepTaxa contains the names of all the OTUs that we want to keep for analysis
#prune phyloseq object to keep only wanted OTUs
deseq.want = prune_taxa(keepTaxa, deseq)
deseq.want
#sucessfully removed sequences, 16941 OTUs remain

#save phyloseq object
save(deseq.want, file="DEV16S.deseq.want.RData")
