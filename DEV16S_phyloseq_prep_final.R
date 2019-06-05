###################
#DEV 16S preprocess 
#Cassandra Wattenburger

#NOTES: 
#This code preprocesses the count data for analysis and formats the data for phyloseq
#Uses files generated in "DEV16S_css_normalize.R"
#Compartment and soil refer to root proximity and system respectively

#Clear workspace and load necessary packages
rm(list=ls())

library("phyloseq")
library("tidyr")

#set working directory to wherever works for you
setwd("Y:/Cassi/Sequencing Data/Argonne 2017/01 DEV Project/02 Github code")


##################################
#Read in CSS normalized count data
#File generated in "DEV16S_css_normalize.R"
otu <- read.csv("DEV16S_cssnorm.csv")
head(otu[,1:10])
dim(otu)

#Read in raw otu count data which includes taxon assignments for every OTU
raw <- read.csv("DEV16S_otu.csv")
head(raw[,1:5])
dim(raw)
head(raw[,220:225])
#only need OTU ID and taxonomy
tax <- raw[,c(1,225)]
head(tax)
dim(tax)

#Remove excess OTUs in tax that do not occur in otu (1781 taxa)
head(tax)
head(otu[,1:5])
tax.filter <- tax[(tax$X.OTU.ID %in% otu$X),]
head(tax.filter)
dim(tax.filter)

#Merge tax.filter and otu so that both have OTUs in same order
otu.tax <- merge(otu, tax.filter, by = 1)
head(otu.tax[,1:5])
dim(otu.tax)
head(otu.tax[,210:215])

#Split taxonomy and otu again
tax.split <- otu.tax[,c(1,215)]
head(tax.split)
dim(tax.split)
otu.split <- otu.tax[,1:214]
head(otu.split[,1:5])
dim(otu.split)

#split up taxonomy at each level
head(tax.split)
dim(tax.split)
tax.levels <- separate(data = tax.split, 
                       col = taxonomy, 
                       into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                       sep=",")
head(tax.levels)
dim(tax.levels)

#Metadata needs to also have same number of samples (some were dropped from data during quality filtering)
#Remove samples that are not present in css normalized count data
meta <- read.csv("DEV_metadata.csv")
head(meta)
dim(meta)

#remove MDP4 from metadata (sample was dropped from count table during quality filtering)
meta.filter <- meta[meta$Sample %in% names(otu),]
head(meta.filter)
dim(meta.filter)
head(otu[,1:5])
dim(otu)

#reorder levels of Root proximity column to bulk, rhizosphere, rhizoplane
levels(meta.filter$Compartment)
meta.filter$Compartment <- factor(meta.filter$Compartment, levels=c("bulk", "rhizosphere", "rhizoplane"))
levels(meta.filter$Compartment)

#create a "Treatment" column for easier labelling later on
meta.filter$Treatment <- with(meta.filter, paste0(Timepoint, sep = " ", Soil, sep = " ", Compartment))
head(meta.filter)
class(meta.filter$Treatment)
meta.filter$Treatment <- as.factor(meta.filter$Treatment)

#write metadata file for phyloseq
write.csv(meta.filter, "DEV16S_phyloseq_metadata.csv")


#######################
#Create phyloseq object

#OTU count data
head(otu.split[,1:5])
rownames(otu.split) <- otu.split$X
head(otu.split[,1:5])
otu.split <- otu.split[,-1]
head(otu.split[,1:5])
class(otu.split)
otu.m <- as.matrix(otu.split)
head(otu.m)
class(otu.m)

OTU = otu_table(otu.m, taxa_are_rows=TRUE)

#Taxonomy data
head(tax.levels)
rownames(tax.levels) <- tax.levels$X
head(tax.levels)
tax.levels <- tax.levels[-1]
head(tax.levels)
class(tax.levels)
tax.m <- as.matrix(tax.levels)
class(tax.m)

TAX = tax_table(tax.m)

#Metdata
head(meta.filter)
rownames(meta.filter) <- meta.filter$Sample
head(meta.filter)
meta.filter <- meta.filter[,-1]
head(meta.filter)

SAM = sample_data(meta.filter, errorIfNULL=TRUE)

#Create phyloseq object from components
physeq = phyloseq(OTU, TAX, SAM)
physeq


##############
#Preprocessing

#Remove taxa with unknown phylum-level taxonomy
badTaxa = subset_taxa(physeq, Phylum=="p__?")
badTaxa
#967 unknown sequences were found
#find chloroplast sequences
badTaxa1 = subset_taxa(physeq, Class=="c__Chloroplast")
badTaxa1
#66 chloroplast sequences need to be removed
#find mitochondrial sequences
badTaxa2 = subset_taxa(physeq, Family=="f__Mitochondria")
badTaxa2
#86 mitochondrial sequences need to be removed
#merge all unwanted data together
badTaxa.all = merge_phyloseq(badTaxa, badTaxa1, badTaxa2)
badTaxa.all
#1119 total taxa need to be removed before analysis

#remove unwanted taxa
removeTaxa = taxa_names(badTaxa.all)
allTaxa = taxa_names(physeq)
keepTaxa = allTaxa[!(allTaxa %in% removeTaxa)]
keepTaxa
#keepTaxa contains the names of all the OTUs that we want to keep for analysis
#prune phyloseq object to keep only wanted OTUs
physeq.want = prune_taxa(keepTaxa, physeq)
physeq.want
#sucessfully removed sequences, 16940 OTUs remain

#save phyloseq object
save(physeq.want, file="DEV16S.physeq.want.RData")

#write file with otu count data of preprocessed data from phyloseq object
OTU.want = as(otu_table(physeq.want), "matrix")
head(OTU.want[,1:5])
write.csv(OTU.want, "DEV16S_phyloseq_otu.csv")
