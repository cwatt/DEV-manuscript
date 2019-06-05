#######################
#DEV ITS mock community
#Cassandra Wattenburger

#NOTES: 
#This code explores the mock community data
#Compartment and soil refer to root proximity and system respectively

#Clear workspace and load necessary packages
rm(list=ls())

library("labdsv")

#set working directory to wherever works for you
setwd("Y:/Cassi/Sequencing Data/Argonne 2017/01 DEV Project/02 Github code")


##############################################

#Read in raw count table, file: DEVITS_otu.csv
raw <- read.csv("DEVITS_otu.csv", row.names = 1)
head(raw[,1:5])
dim(raw)

#Subset only mock samples with OTU labels
mock <- raw[,c("MF1", "MF2", "MF3")]
head(mock)

#Remove 1 and 0 abundance OTUs
mock.trans <- t(mock)
head(mock.trans[,1:5])
mock.nosingle <- dropspc(mock.trans, 1)
dim(mock)
dim(mock.nosingle)
head(mock.nosingle)
tail(mock.nosingle)
#22 OTUs remaning, this is more species than are in the actual mock community
#likely due to uncertainty in taxonomic assignment, contamination, and/or well hopping

#Add in taxonomic info
mock.nosingle.trans <- t(mock.nosingle)
head(mock.nosingle.trans)
dim(mock.nosingle.trans)
dim(raw)
#I include column 215 here because I want it to be considered a dataframe for adding row names
tax <- raw[,c(215:216)]
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

write.csv(mock.tax, "DEVITS_mock.csv")

#Attempt to find the taxa present in the mock communities
#From Matt Bakker (Staggered A):
#1 Alternia alternata, relative abundance 1
#2 Aspergillus falvus, 10
#3 Neosartorya fischeri, 10
#4 Penicillium expansum, 10
#5 Candida apicola, 10
#6 Saccharomyces cerevisiae, 100
#7 Claviceps purpurea, 10
#8 Trichoderma reesei, 10
#9 Sarocladium zeae, 10
#10 Fusarium graminearum, 100
#11 Fusarium oxysporum, 10
#12 Fusarium vesticillioides, 500
#13 Soitoella complicata, 2000
#14 Rhizoctonia solani, 5
#15 Naganishia albida, 1000
#16 Chytriomyces hyalinus, 10
#17 Rhizophagus irregularis, 1
#18 Mortierella verticillata, 10
#19 Rhizomucor miehei, 5

#Subset Alternia alternata
mock1g <- subset(mock.tax, Genus=="g__Alternaria")
mock1g
#OTU_12, identified at species level

#Subset Aspergillus flavus
mock2g <- subset(mock.tax, Genus=="g__Aspergillus")
mock2g
mock2f <- subset(mock.tax, Family=="f__Aspergillaceae")
mock2f
mock2o <- subset(mock.tax, Order=="o__Eurotiales")
mock2o
mock2c <- subset(mock.tax, Class=="c__Eurotiomycetes")
mock2c
#not found

#Subset Neosartorya fischeri
mock3g <- subset(mock.tax, Genus=="g__Neosartorya")
mock3g
mock3f <- subset(mock.tax, Family=="f__Trichocomaceae")
mock3f
mock3o <- subset(mock.tax, Order=="o__Eurotiales")
mock3o
mock3c <- subset(mock.tax, Class=="c__Eurotiomycetes")
mock3c
#not found
#Look for synonymous names
mock3g.2 <- subset(mock.tax, Genus=="g__Sartorya")
mock3g.2
#not found

#Subset Penicillium expansum
mock4g <- subset(mock.tax, Genus=="g__Penicillium")
mock4g
mock4f <- subset(mock.tax, Family=="f__Trichocomaceae")
mock4f
mock4o <- subset(mock.tax, Order=="o__Eurotiales")
mock4o
mock4c <- subset(mock.tax, Class=="c__Eurotiomycetes")
mock4c
#not found
#Look for synonym
mock4g.2 <- subset(mock.tax, Genus=="g__Coremium")
mock4g.2
#not found

#Subset Candida apicola
mock5g <- subset(mock.tax, Genus=="g__Candida")
mock5g
mock5f <- subset(mock.tax, Family=="f__incertae_sedis")
mock5f
mock5o <- subset(mock.tax, Order=="o__Saccharomycetales")
mock5o
#not found

#Subset Saccharomyces cerevisiae
mock6g <- subset(mock.tax, Genus=="g__Saccharomyces")
mock6g
#OTU_580, correctly identified at genus level

#Subset Claviceps purpurea
mock7g <- subset(mock.tax, Genus=="g__Claviceps")
mock7g
mock7f <- subset(mock.tax, Family=="f__Clavicipitaceae")
mock7f
#OTU_1141, identified at family level

#Subset Trichoderma reesei, Sarocladium zeae, Fusarium graminearum, Fusarium oxysporum, 
#Fusarium vesticillioides (similar taxonomies)
mock8.12o <- subset(mock.tax, Order=="o__Hypocreales")
mock8.12o
#Fusarium verticillioides = OTU_3, likely, at genus level
#Fusarium oxysporum = OTU_119, likely, at genus level
#Fusarium graminearum not found?
##above three likely binned together
#Trichoderma reesei = OTU_433?, identified at family level
#Sarocladium zeae not found?
##NOTE: Sarocladium zeae was excluded from the mock community that I was given
#How does that change theoretical relative abundance?

#I'm going to combine the Fusariums into one because they were not distinguished between one another
#Combine ORU_119 and OTU_3
fusarium <- mock8.12o[-c(1,4),-c(4:10)]
fusarium
fus.tax <- mock8.12o[-c(1,3,4),4:10]
fus.tax
fus.sum <- colSums(fusarium)
fus.sum <- as.data.frame(fus.sum)
fus.sum <- t(fus.sum)
fus.sum
fus.sum <- cbind(fus.sum, fus.tax)
fus.sum
row.names(fus.sum) <- "Fusarium"
fus.sum

#Subset reesei
reesei <- mock8.12o[4,]
reesei

#Subset Saitoella complicata
mock13g <- subset(mock.tax, Genus=="g__Saitoella")
mock13g
#OTU_58, identified at species level

#Subset Rhizoctonia solani
mock14g <- subset(mock.tax, Genus=="g__Rhizoctonia")
mock14g
mock14f <- subset(mock.tax, Family=="f__Ceratobasidiaceae")
mock14f
#OTU_900, identified at Family level

#Subset Naganishia albida
mock15g <- subset(mock.tax, Genus=="g__Naganishia")
mock15g
mock15f <- subset(mock.tax, Family=="f__Filobasidiaceae")
mock15f
mock15o <- subset(mock.tax, Order=="o__Filobasidiales")
mock15o
mock15c <- subset(mock.tax, Class=="c__Tremellomycetes")
mock15c
#OTU_21, misidentified as Cryptococcus adeliensis, first correctly identified at Genus level
#Cryptococcus is same genus, albida and adeliensis are closely related
mock15g.2 <- subset(mock.tax, Genus=="g__Cryptococcus")
mock15g.2 <- mock15g.2[-2,]
mock15g.2

#Subset Chytriomyces hyalinus 
mock16g <- subset(mock.tax, Genus=="g__Chytriomyces")
mock16g
mock16f <- subset(mock.tax, Family=="f__Chytridiaceae")
mock16f
mock16o <- subset(mock.tax, Order=="o__Chytridiales")
mock16o
mock16c <- subset(mock.tax, Class=="c__Chytridiomycetes")
mock16c
mock16p <- subset(mock.tax, Phylum=="p__Chytridiomycota")
mock16p
#not found

#Subset Rhizophagus irregularis
mock17g <- subset(mock.tax, Genus=="g__Rhizophagus")
mock17g
mock17f <- subset(mock.tax, Family=="f__Glomeraceae")
mock17f
mock17o <- subset(mock.tax, Order=="o__Glomerales")
mock17o
mock17c <- subset(mock.tax, Class=="c__Glomeromycetes")
mock17c
mock17p <- subset(mock.tax, Phylum=="p__Glomeromycota")
mock17p
#not found

#Subset Mortierella verticillata
mock18g <- subset(mock.tax, Genus=="g__Mortierella")
mock18g
#OTU_427, misidentified as Mortierella humilis, first correctly identified at genus level
mock18g <- mock18g[-1,]
mock18g

#Subset Rhizomucor miehei
mock19g <- subset(mock.tax, Genus=="g__Rhizomucor")
mock19g
mock19f <- subset(mock.tax, Family=="f__Lichtheimiaceae")
mock19f
#OTU_3005, identified at family

#Create dataframe for species that were not found
notfound <- data.frame(MF1 = 0, MF2 = 0, MF3 = 0, Kingdom = "no", Phylum = "no", Class = "no", Order = "no", Family = "no", Genus = "no", Species = "no")
notfound

#Combine OTUs into single dataframe
mock.otus <- rbind(mock1g,
                   notfound,
                   notfound,
                   notfound,
                   notfound,
                   mock6g,
                   mock7f,
                   reesei,
                   fus.sum,
                   mock13g,
                   mock14f,
                   mock15g.2,
                   notfound,
                   notfound,
                   mock18g,
                   mock19f)
mock.otus
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
mean <- as.data.frame(mean)
mean

#Make a nice summary table
#Theoretical values have changed from what Dr. Bakker gave me because Sarocladium zeae was excluded and 
#Fusarium had to be combined
#Calculated these outside of R in excel because I don't have that kind of time
theoretical <- c(0.0263, 0.263, 0.263, 0.263, 0.263, 2.63, 0.263, 0.263, 
                 16.04, 52.60, 0.132, 26.3, 0.263, 0.0263, 0.263, 0.132)
theoretical
theoretical <- as.data.frame(theoretical)

dim(theoretical)
dim(mean)
compare <- cbind(mean, theoretical)
compare

#The taxonomic level correctly identified at
taxlevel <- c("Species", "Not Found", "Not Found", "Not Found", "Not Found", "Genus",
              "Family", "Family", "Genus", "Species", "Family", "Genus", "Not Found", 
              "Not Found", "Genus", "Family")
taxlevel <- as.data.frame(taxlevel)
taxlevel

#Species names
species <- c("Alternia alternata",
             "Aspergillus flavus",
             "Neosartorya fischeri",
             "Penicillium expansum",
             "Candida apicola",
             "Saccharomyces cerevisiae",
             "Claviceps purpurea",
             "Trichoderma reesei",
             "Fusarium graminearum, oxysporum, and vesticilloides",
             "Soitoella complicata",
             "Rhizoctonia solani",
             "Naganishia albida",
             "Chytriomyces hyalinus",
             "Rhizophagus irregularis",
             "Mortierella verticillata",
             "Rhizomucor miehei")
species <- as.data.frame(species)
species

summary <- cbind(compare, species, taxlevel)
summary

labels <- c("Actual", "Theoretical", "Species", "Identified at")

colnames(summary) <- labels
summary

#Save the summary table
write.csv(summary, "DEVITS_mock_summary.csv")

