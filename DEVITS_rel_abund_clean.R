###########################
#DEV ITS relative abundance
#Cassandra Wattenburger

#NOTES: 
#This code produces relative abundance data at the phylum level
#Uses data generated in "DEVITS_phyloseq_prep_clean.R"
#Compartment and soil refer to root proximity and system respectively

#Clear workspace and load necessary packages
rm(list=ls())

#load in necessary packages
library("phyloseq")
library("ggplot2")
library("tidyr")
library("plyr")
library("reshape2")

#set working directory to wherever works for you
setwd("Y:/Cassi/Sequencing Data/Argonne 2017/01 DEV Project/02 Github code")


############################
#Relative Abundance of phyla

load("DEVITS.physeq.want.RData")

#Aggregate taxa to phylum level
physeq.pglom = tax_glom(physeq.want, taxrank="Phylum")
physeq.pglom

#Calculate relative abundance
relabund  = transform_sample_counts(physeq.pglom, function(x) x / sum(x))

#Convert aggregated data to non-phyloseq objects (easier graphing)

#OTU relative abundance
relabund.otu = as(otu_table(relabund), "matrix")
head(relabund.otu[,1:5])
relabund.trans <- t(relabund.otu)
head(relabund.trans[,1:5])
dim(relabund.trans)

#Taxonomy
relabund.tax = as.data.frame(tax_table(relabund))
head(relabund.tax)
dim(relabund.tax)
relabund.tax <- relabund.tax[,-(3:7)]
head(relabund.tax)

#Metadata
relabund.meta = sample_data(relabund)
relabund.meta <- as.data.frame(relabund.meta)
head(relabund.meta)
dim(relabund.meta)
relabund.meta$SoilComp <- with(relabund.meta, paste0(Soil, sep = " ", Compartment))
head(relabund.meta)

#Combine metadata and otu relative abundance
head(relabund.trans[,1:5])
head(relabund.meta)
rel.otu.meta <- merge(relabund.meta, relabund.trans, by=0)
head(rel.otu.meta[,1:10])
row.names(rel.otu.meta) <- rel.otu.meta$Row.names
head(rel.otu.meta[,1:10])
rel.otu.meta <- rel.otu.meta[,-1]
head(rel.otu.meta[,1:10])

#Convert to long format by melting data frame
#Remove block variable, don't need
rel.otu.meta <- rel.otu.meta[,-4]
head(rel.otu.meta[,1:10])
rel.melt <- melt(rel.otu.meta, id.cars=c("Soil", "Compartment", "Timepoint", "Treatment", "SoilComp"))
head(rel.melt)

#add taxonomic data
head(rel.melt)
head(relabund.tax)
relabund.tax$otu <- row.names(relabund.tax)
rel.melt.tax <- merge(relabund.tax, rel.melt, by.x=3, by.y=6)
head(rel.melt.tax)

#Calculate average relative abundance for each soil
soil.summary <- ddply(rel.melt.tax, c("Soil", "Phylum", "otu"), summarize, mean=mean(value), sd=sd(value), se=sd/sqrt(length(value)))
head(soil.summary)

#Calculate average relative abundance for each compartment
comp.summary <- ddply(rel.melt.tax, c("Compartment", "Phylum", "otu"), summarize, mean=mean(value), sd=sd(value), se=sd/sqrt(length(value)))
head(comp.summary)

#Calculate average relative abundance for each treatment
trt.summary <- ddply(rel.melt.tax, c("Treatment", "Phylum", "otu"), summarize, mean=mean(value), sd=sd(value), se=sd/sqrt(length(value)))
head(trt.summary)


####################
#Graphs for Figure 1

#Remove p__ from phylum labels
trt.summary$Phylum <- gsub("p__", "", trt.summary$Phylum)

#Bulk soil
head(trt.summary)
#Remove rows that aren't bulk soil
trt.summary.b1 <- subset(trt.summary, Treatment=="TP1 Conv. bulk")
trt.summary.b2 <- subset(trt.summary, Treatment=="TP2 Conv. bulk")
trt.summary.b3 <- subset(trt.summary, Treatment=="TP3 Conv. bulk")
trt.summary.b4 <- subset(trt.summary, Treatment=="TP4 Conv. bulk")
trt.summary.b5 <- subset(trt.summary, Treatment=="TP1 Div. bulk")
trt.summary.b6 <- subset(trt.summary, Treatment=="TP2 Div. bulk")
trt.summary.b7 <- subset(trt.summary, Treatment=="TP3 Div. bulk")
trt.summary.b8 <- subset(trt.summary, Treatment=="TP4 Div. bulk")
trt.summary.b <- rbind(trt.summary.b1,trt.summary.b2,trt.summary.b3,trt.summary.b4,
                       trt.summary.b5,trt.summary.b6,trt.summary.b7,trt.summary.b8)
head(trt.summary.b)

#created large palette using http://tools.medialab.sciences-po.fr/iwanthue/
palette <- c("#00a8e7",
             "#9c3700",
             "#0165bf",
             "#76b330",
             "#ad4dbe",
             "#c0df8d",
             "#ff947f")

#graph bulk
trt.b.graph <- ggplot(trt.summary.b, aes(x=Treatment, y=mean, fill=reorder(Phylum, mean))) +
  geom_col(color="Black", width=0.7) +
  ggtitle("") +
  labs(x="", y="") +
  scale_fill_manual(values=palette) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  guides(fill=guide_legend(title="Phyla", reverse=TRUE)) +
  theme_classic() +
  theme(plot.title = element_text(size=18),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14),
        legend.text=element_text(size=10),
        legend.justification = "top")
trt.b.graph

#Rhizosphere soil

head(trt.summary)
#Remove rows that aren't rhizosphere soil
trt.summary.rs1 <- subset(trt.summary, Treatment=="TP1 Conv. rhizosphere")
trt.summary.rs2 <- subset(trt.summary, Treatment=="TP2 Conv. rhizosphere")
trt.summary.rs3 <- subset(trt.summary, Treatment=="TP3 Conv. rhizosphere")
trt.summary.rs4 <- subset(trt.summary, Treatment=="TP4 Conv. rhizosphere")
trt.summary.rs5 <- subset(trt.summary, Treatment=="TP1 Div. rhizosphere")
trt.summary.rs6 <- subset(trt.summary, Treatment=="TP2 Div. rhizosphere")
trt.summary.rs7 <- subset(trt.summary, Treatment=="TP3 Div. rhizosphere")
trt.summary.rs8 <- subset(trt.summary, Treatment=="TP4 Div. rhizosphere")
trt.summary.rs <- rbind(trt.summary.rs1,trt.summary.rs2,trt.summary.rs3,trt.summary.rs4,
                        trt.summary.rs5,trt.summary.rs6,trt.summary.rs7,trt.summary.rs8)
head(trt.summary.rs)

palette2 <- c("#00a8e7",
              "#0165bf",
              "#9c3700",
              "#76b330",
              "#c0df8d",
              "#ad4dbe",
              "#ff947f")
#needed to change order to keep coloring and taxa consistent

#graph rhizosphere
trt.rs.graph <- ggplot(trt.summary.rs, aes(x=Treatment, y=mean, fill=reorder(Phylum, mean))) +
  geom_col(color="Black", width=0.7) +
  ggtitle("") +
  labs(x="", y="") +
  scale_fill_manual(values=palette2) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  guides(fill=guide_legend(title="Phyla", ncol=1, reverse=TRUE)) +
  theme_classic() +
  theme(plot.title = element_text(size=18),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14),
        legend.text=element_text(size=10))
trt.rs.graph

#Rhizoplane

head(trt.summary)
#Remove rows that aren't bulk soil
trt.summary.rp1 <- subset(trt.summary, Treatment=="TP1 Conv. rhizoplane")
trt.summary.rp2 <- subset(trt.summary, Treatment=="TP2 Conv. rhizoplane")
trt.summary.rp3 <- subset(trt.summary, Treatment=="TP3 Conv. rhizoplane")
trt.summary.rp4 <- subset(trt.summary, Treatment=="TP4 Conv. rhizoplane")
trt.summary.rp5 <- subset(trt.summary, Treatment=="TP1 Div. rhizoplane")
trt.summary.rp6 <- subset(trt.summary, Treatment=="TP2 Div. rhizoplane")
trt.summary.rp7 <- subset(trt.summary, Treatment=="TP3 Div. rhizoplane")
trt.summary.rp8 <- subset(trt.summary, Treatment=="TP4 Div. rhizoplane")
trt.summary.rp <- rbind(trt.summary.rp1,trt.summary.rp2,trt.summary.rp3,trt.summary.rp4,
                        trt.summary.rp5,trt.summary.rp6,trt.summary.rp7,trt.summary.rp8)
head(trt.summary.rp)

palette3 <- c("#9c3700",
              "#0165bf",
              "#76b330",
              "#00a8e7",
              "#c0df8d",
              "#ad4dbe",
              "#ff947f")
#changed color order for consistency

trt.rp.graph <- ggplot(trt.summary.rp, aes(x=Treatment, y=mean, fill=reorder(Phylum, mean))) +
  geom_col(color="Black", width=0.7) +
  ggtitle("") +
  labs(x="", y="") +
  scale_fill_manual(values=palette3) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  guides(fill=guide_legend(title="Phyla", ncol=1, reverse=TRUE)) +
  theme_classic() +
  theme(plot.title = element_text(size=18),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14),
        legend.text=element_text(size=10))
trt.rp.graph


#############################
#Graph Supplementary Figure 2

#System graph
soil.summary$Phylum <- gsub("p__", "", soil.summary$Pylum) #why doesn't this work???

relabund.soil.supp <- ggplot(soil.summary, aes(x=Soil, y=mean, fill=reorder(Phylum, mean))) +
  geom_col(color="Black", width=0.7) +
  ggtitle("") +
  labs(x="", y="") +
  scale_fill_manual(values=palette) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  guides(fill=guide_legend(title="Phyla", ncol=1, reverse=TRUE)) +
  theme_classic() +
  theme(plot.title = element_text(size=18),
        axis.text.x = element_text(size=12, angle = 45, hjust = 1),
        axis.title=element_text(size=14),
        legend.text=element_text(size=10))
relabund.soil.supp

#Root proximity graph

comp.summary$Phylum <- gsub("p__", "", comp.summary$Pylum) #????????????????????

palette4 <- c("#00a8e7",
              "#9c3700",
              "#76b330",
              "#0165bf",
              "#ad4dbe",
              "#c0df8d",
              "#ff947f")
#for color coding consistency

relabund.comp.supp <- ggplot(comp.summary, aes(x=Compartment, y=mean, fill=reorder(Phylum, mean))) +
  geom_col(color="Black", width=0.7) +
  ggtitle("") +
  labs(x="", y="") +
  scale_fill_manual(values=palette4) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  guides(fill=guide_legend(title="Phyla", ncol=1, reverse=TRUE)) +
  theme_classic() +
  theme(plot.title = element_text(size=18),
        axis.text.x = element_text(size=12, angle = 45, hjust = 1),
        axis.title=element_text(size=14),
        legend.text=element_text(size=10))
relabund.comp.supp

