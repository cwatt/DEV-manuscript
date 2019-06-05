
# Cassandra Wattenburger
# 01/18/19

# Cleaned up scripts to generate relative abundance data and figure 1.

############################################

# load necessary packages
library("phyloseq")
library("ggplot2")
library("tidyr")
library("plyr")
library("reshape2")
library("cowplot")

setwd("Y:/Cassi/Sequencing Data/Argonne 2017/01 DEV Project/03 Final Scripts and figures")


##########
# 16S data

rm(list=ls())

########
# Format

load("DEV16S.physeq.want.RData")
# file generated in "______________"

# Aggregate taxa to phylum level
physeq.pglom = tax_glom(physeq.want, taxrank="Phylum")
physeq.pglom

# Calculate relative abundance
relabund  = transform_sample_counts(physeq.pglom, function(x) x / sum(x))

# Convert aggregated data to non-phyloseq objects (easier graphing):

# phylum relative abundance data
relabund.otu = as(otu_table(relabund), "matrix")
head(relabund.otu[,1:5])
relabund.trans <- t(relabund.otu)
head(relabund.trans[,1:5])
dim(relabund.trans)

# taxonomic data
relabund.tax = as.data.frame(tax_table(relabund))
head(relabund.tax)
dim(relabund.tax)
relabund.tax <- relabund.tax[,-(3:7)]
head(relabund.tax)

# metadata
relabund.meta = sample_data(relabund)
relabund.meta <- as.data.frame(relabund.meta)
head(relabund.meta)
dim(relabund.meta)
head(relabund.meta)

# combine metadata and relative abundance data
head(relabund.trans[,1:5])
head(relabund.meta)
rel.otu.meta <- merge(relabund.meta, relabund.trans, by=0)
head(rel.otu.meta[,1:10])
row.names(rel.otu.meta) <- rel.otu.meta$Row.names
head(rel.otu.meta[,1:10])
rel.otu.meta <- rel.otu.meta[,-1]
head(rel.otu.meta[,1:10])

# convert to long format by melting data frame
# remove block variable, don't need
rel.otu.meta <- rel.otu.meta[,-4]
head(rel.otu.meta[,1:10])
rel.melt <- melt(rel.otu.meta, id.cars=c("Treatment"))
head(rel.melt)

# add taxonomic data
head(rel.melt)
head(relabund.tax)
relabund.tax$otu <- row.names(relabund.tax)
rel.melt.tax <- merge(relabund.tax, rel.melt, by.x=0, by.y=5)
head(rel.melt.tax)

# calculate average relative abundance for each treatment
trt.summary <- ddply(rel.melt.tax, c("Treatment", "Phylum", "otu"), summarize, mean=mean(value), sd=sd(value), se=sd/sqrt(length(value)))
head(trt.summary)

# remove "p__" etc. in phylum names
trt.summary$Phylum <- gsub("p__", "", trt.summary$Phylum)
trt.summary$Phylum <- gsub("_cls_", "", trt.summary$Phylum)
trt.summary$Phylum <- gsub("_", " ", trt.summary$Phylum)
trt.summary$Phylum <- gsub("\\?", "Unknown", trt.summary$Phylum)
head(trt.summary)

# create dev stage label
trt.summary$xlab1 <- trt.summary$Treatment
trt.summary$xlab1 <- gsub("TP1 .+", "V4", trt.summary$xlab1)
trt.summary$xlab1 <- gsub("TP2 .+", "V11", trt.summary$xlab1)
trt.summary$xlab1 <- gsub("TP3 .+", "R2", trt.summary$xlab1) 
trt.summary$xlab1 <- gsub("TP4 .+", "R5", trt.summary$xlab1)
head(trt.summary)
tail(trt.summary)

# create system label
trt.summary$xlab2 <- trt.summary$Treatment
trt.summary$xlab2 <- gsub("TP[0-9] Conv\\. .+", "C", trt.summary$xlab2)
trt.summary$xlab2 <- gsub("TP[0-9] Div\\. .+", "D", trt.summary$xlab2)
head(trt.summary)
tail(trt.summary)

# create x label
trt.summary$xlab <- paste(trt.summary$xlab1, trt.summary$xlab2)
head(trt.summary)
tail(trt.summary)

#change levels of x label
trt.summary$xlab <- factor(trt.summary$xlab, levels = c("V4 C", "V4 D", "V11 C", "V11 D", "R2 C", "R2 D", "R5 C", "R5 D"))
levels(trt.summary$xlab)

# create root compartment label
head(trt.summary)
trt.summary$compartment <- trt.summary$Treatment
trt.summary$compartment <- gsub(".+ bulk", "Bulk soil", trt.summary$compartment)
trt.summary$compartment <- gsub(".+ rhizosphere", "Rhizosphere soil", trt.summary$compartment)
trt.summary$compartment <- gsub(".+ rhizoplane", "Rhizoplane", trt.summary$compartment)
tail(trt.summary)

# change levels of root compartment
trt.summary$compartment <- factor(trt.summary$compartment, levels=c("Bulk soil", "Rhizosphere soil", "Rhizoplane"))

#######
# Graph

#Palette, created using http://tools.medialab.sciences-po.fr/iwanthue/
#I'm not including all 54 prokaryotic phyla in legend because that's crazy, hence many white placeholders in palette
palette <- c("#ffffff","#ffffff","#ffffff","#ffffff","#ffffff","#ffffff","#ffffff","#ffffff",
             "#ffffff","#ffffff","#ffffff","#ffffff","#ffffff","#ffffff","#ffffff","#ffffff",
             "#ffffff","#ffffff","#ffffff","#ffffff","#ffffff","#ffffff","#ffffff","#ffffff",
             "#ffffff","#ffffff","#ffffff","#ffffff","#ffffff","#ffffff","#ffffff","#ffffff",
             "#ffffff","#ffffff","#ffffff","#ffffff","#ffffff","#ffffff","#ffffff",
             "#751800",
             "#c471e7",
             "#95da5a",
             "#2468a8",
             "#ff9b63",
             "#00ebc8",
             "#916200",
             "#2052c5",
             "#7f9c00",
             "#60004b",
             "#bbe091",
             "#f26ad1",
             "#00913d",
             "#fc4e8a",
             "#63b1ff")

list <- c("Proteobacteria", "Acidobacteria", "Bacteroidetes", "Verrucomicrobia", "Actinobacteria", "Thaumarchaeota", "Planctomycetes", "Chloroflexi", "Gemmatimonadetes", "Nitrospirae", "Firmicutes", "Armatimonadetes", "Latescibacteria", "Elusimicrobia", "Chlorobi")

relabund.16s.graph <- ggplot(trt.summary, aes(x=xlab, y=mean, fill=reorder(Phylum, mean))) +
  geom_col(color="Black", width=0.7) +
  facet_wrap(~compartment) +
  labs(x="", y="Relative abundance") +
  scale_fill_manual(values=palette, breaks=list) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  guides(fill=guide_legend(title="Phylum")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size=8),
        axis.title=element_text(size=8),
        strip.text.x = element_text(size = 8, margin=margin(3,3,3,3)),
        legend.text=element_text(size=8),
        legend.title=element_text(size=8),
        panel.spacing = unit(1, "lines"),
        plot.margin = unit(c(0, 0, 0, 0.5), "cm"))
relabund.16s.graph

save(relabund.16s.graph, file="relabund.16s.graph.RData")


##########
# ITS data

rm(list=ls())

########
# Format

load("DEVITS.physeq.want.RData")

# Aggregate taxa to class level
physeq.cglom = tax_glom(physeq.want, taxrank="Class")
physeq.cglom

# Calculate relative abundance
relabund  = transform_sample_counts(physeq.cglom, function(x) x / sum(x))

# Convert aggregated data to non-phyloseq objects (easier graphing):

# class relative abundance data
relabund.otu = as(otu_table(relabund), "matrix")
head(relabund.otu[,1:5])
relabund.trans <- t(relabund.otu)
head(relabund.trans[,1:5])
dim(relabund.trans)

# taxonomic data
relabund.tax = as.data.frame(tax_table(relabund))
head(relabund.tax)
dim(relabund.tax)
relabund.tax <- relabund.tax[,-(4:7)]
head(relabund.tax)

# metadata
relabund.meta = sample_data(relabund)
relabund.meta <- as.data.frame(relabund.meta)
head(relabund.meta)
dim(relabund.meta)

# combine metadata and relative abundance data
head(relabund.trans[,1:5])
head(relabund.meta)
rel.otu.meta <- merge(relabund.meta, relabund.trans, by=0)
head(rel.otu.meta[,1:10])
row.names(rel.otu.meta) <- rel.otu.meta$Row.names
head(rel.otu.meta[,1:10])
rel.otu.meta <- rel.otu.meta[,-1]
head(rel.otu.meta[,1:10])

# convert to long format by melting data frame
# remove block variable, don't need
rel.otu.meta <- rel.otu.meta[,-4]
head(rel.otu.meta[,1:10])
rel.melt <- melt(rel.otu.meta, id.cars=c("Treatment"))
head(rel.melt)

# add taxonomic data
head(rel.melt)
head(relabund.tax)
relabund.tax$otu <- row.names(relabund.tax)
rel.melt.tax <- merge(relabund.tax, rel.melt, by.x=4, by.y=5)
head(rel.melt.tax)

# calculate average relative abundance for each treatment
trt.summary <- ddply(rel.melt.tax, c("Treatment", "Phylum", "Class", "otu"), summarize, mean=mean(value), sd=sd(value), se=sd/sqrt(length(value)))
head(trt.summary)

# remove "c__" etc. in class names
trt.summary$Class <- gsub("c__", "", trt.summary$Class)
trt.summary$Class <- gsub("_cls_Incertae_sedis", "", trt.summary$Class)
trt.summary$Class <- gsub("\\?", "Unknown", trt.summary$Class)
trt.summary$Phylum <- gsub("p__", "", trt.summary$Phylum)
head(trt.summary)

# abbrev. phylum names
trt.summary$Phylum <- gsub("Ascomycota", "Asco.", trt.summary$Phylum)
trt.summary$Phylum <- gsub("Zygomycota", "Zygo.", trt.summary$Phylum)
trt.summary$Phylum <- gsub("Basidiomycota", "Basid.", trt.summary$Phylum)
trt.summary$Phylum <- gsub("Chytridiomycota", "Chytrid.", trt.summary$Phylum)

# combine phylum and class labels
trt.summary$phyclass <- paste(trt.summary$Phylum, trt.summary$Class)
head(trt.summary)

# create dev stage label
trt.summary$xlab1 <- trt.summary$Treatment
trt.summary$xlab1 <- gsub("TP1 .+", "V4", trt.summary$xlab1)
trt.summary$xlab1 <- gsub("TP2 .+", "V11", trt.summary$xlab1)
trt.summary$xlab1 <- gsub("TP3 .+", "R2", trt.summary$xlab1) 
trt.summary$xlab1 <- gsub("TP4 .+", "R5", trt.summary$xlab1)
head(trt.summary)
tail(trt.summary)

# create system label
trt.summary$xlab2 <- trt.summary$Treatment
trt.summary$xlab2 <- gsub("TP[0-9] Conv\\. .+", "C", trt.summary$xlab2)
trt.summary$xlab2 <- gsub("TP[0-9] Div\\. .+", "D", trt.summary$xlab2)
head(trt.summary)
tail(trt.summary)

# create x label
trt.summary$xlab <- paste(trt.summary$xlab1, trt.summary$xlab2)
head(trt.summary)
tail(trt.summary)

#change levels of x label
trt.summary$xlab <- factor(trt.summary$xlab, levels = c("V4 C", "V4 D", "V11 C", "V11 D", "R2 C", "R2 D", "R5 C", "R5 D"))
levels(trt.summary$xlab)

# create root compartment label
head(trt.summary)
trt.summary$compartment <- trt.summary$Treatment
trt.summary$compartment <- gsub(".+ bulk", "Bulk soil", trt.summary$compartment)
trt.summary$compartment <- gsub(".+ rhizosphere", "Rhizosphere soil", trt.summary$compartment)
trt.summary$compartment <- gsub(".+ rhizoplane", "Rhizoplane", trt.summary$compartment)
tail(trt.summary)

# change levels of root compartment
trt.summary$compartment <- factor(trt.summary$compartment, levels=c("Bulk soil", "Rhizosphere soil", "Rhizoplane"))
head(trt.summary)

#######
# Graph

#Palette, created using http://tools.medialab.sciences-po.fr/iwanthue/

palette <- c("#bd98da",
             "#b36ba7",
             "#6698e0",
             "#c74cba",
             "#757bb2",
             "#ca6d25",
             "#6d7fdd",
             "#ca6b5d",
             "#ac60dd",
             "#bc6871",
             "#5c59d4",
             "#a68657",
             "#8f4d99",
             "#b49339",
             "#da80a8",
             "#aab82f",
             "#dd80d3",
             "#5ba664",
             "#b77dde",
             "#4c916b",
             "#6d4ab0",
             "#43bea2",
             "#da4c23",
             "#4fa8d4",
             "#d43f69",
             "#5560aa",
             "#d9a02c",
             "#ce4295",
             "#56b642",
             "#9e466d",
             "#79973e",
             "#637df2",
             "#b17c41",
             "#d14b3f",
             "#3bbcc3")

list <- c("Zygo. Mortierellomycotina", "Asco. Sordariomycetes", "Asco. Dothideomycetes", "Asco. Unknown", "Basid. Agaricomycetes", "Basid. Tremellomycetes", "Asco. Leotiomycetes", "Asco. Eurotiomycetes", "Asco. Saccharomycetes", "Basid. Microbotryomycetes", "Basid. Unknown", "Zygo. Unknown", "Chytrid. Chytridiomycetes", "Zygo. Mucoromycotina", "Asco. Pezizomycotina")

relabund.its.graph <- ggplot(trt.summary, aes(x=xlab, y=mean, fill=reorder(phyclass, mean))) +
  geom_col(color="Black", width=0.7) +
  facet_wrap(~compartment) +
  labs(x="", y="Relative abundance") +
  scale_fill_manual(values=palette, breaks=list) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  guides(fill=guide_legend(title="Phylum Class")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        strip.text.x = element_text(size = 8, margin=margin(3,3,3,3)),
        legend.text=element_text(size=8),
        legend.title=element_text(size=8),
        panel.spacing = unit(1, "lines"),
        plot.margin = unit(c(0, 0, 0, 0.5), "cm"))
relabund.its.graph

save(relabund.its.graph, file="relabund.its.graph.RData")


#################
# Create figure 1

rm(list=ls())

load("relabund.16s.graph.RData")
load("relabund.its.graph.RData")

# extract legends and reformat
legend.16s <- get_legend(relabund.16s.graph + theme(legend.key.size = unit(0.25, "cm")))

legend.its <- get_legend(relabund.its.graph + theme(legend.key.size = unit(0.25, "cm")))

# remove legends from graphs
relabund.16s.graph <- relabund.16s.graph + theme(legend.position="none")
relabund.its.graph <- relabund.its.graph + theme(legend.position="none")

legenda <- plot_grid(legend.16s, NULL, ncol=1, rel_heights=c(1,0.5))
legendb <- plot_grid(legend.its, NULL, ncol=1, rel_heights=c(1,0.5))

# fig 1a
fig1a <- plot_grid(relabund.16s.graph, legenda, NULL, rel_widths=c(1,0.34,0.06), nrow=1)
fig1a

# fig 1b
fig1b <- plot_grid(relabund.its.graph, legendb, rel_widths=c(1,0.4), nrow=1)
fig1b

# figure 1
fig1 <- plot_grid(fig1a, fig1b, labels=c("A","B"), label_size=14, nrow=2)
fig1

ggsave("figure1.tiff", plot=fig1, scale=1, width=7, height=5, units="in")

