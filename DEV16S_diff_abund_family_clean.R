#########################################
#DEV 16S differential abundance family level
#Cassandra Wattenburger

#NOTES: 
#This code analyzes differential abundance at the family level using DESeq2
#This script uses files created in "DEV16S_deseq2_prep_clean.R"
#Compartment and soil refer to root proximity and system respectively

#Clear workspace and load necessary packages
rm(list=ls())

library("phyloseq")
library("DESeq2")
library("ggplot2")
library("plyr")

#set working directory to wherever works for you
setwd("Y:/Cassi/Sequencing Data/Argonne 2017/01 DEV Project/02 Github code")


############
#DESeq2 test

#Load in phyloseq object, file: DEV16S.deseq2.want.Rdata
load("DEV16S.deseq.want.Rdata")
deseq.want

#Aggregate OTUs to family level
deseq.fglom = tax_glom(deseq.want, taxrank="Family")
deseq.fglom

#create DESeq2 object from phyloseq object
deseq = phyloseq_to_deseq2(deseq.fglom, ~ group)
deseq

test = DESeq(deseq, test="Wald", fitType="parametric")
resultsNames(test)

#set alpha level to 0.01
#this limits family-wise error to 1%
alpha <- 0.01


#######################
#Rhizoplane Comparisons

############
#Timepoint 1

#Conv. bulk vs Conv. rhizoplane
#For contrasts, the treatment you want to act as the control should be the last listed
##here I set TP1Conv.bulk last because I want the bulk soil to be the baseline for comparison
result.1conv.brp = results(test, contrast=c("group", "TP1Conv.rhizoplane", "TP1Conv.bulk"))
mcols(result.1conv.brp, use.names=TRUE)
#Extract only significant results, at alpha level 0.01 (conservative)
sig.1conv.brp = result.1conv.brp[which(result.1conv.brp$padj < alpha),]
sig.1conv.brp = cbind(as(sig.1conv.brp, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.1conv.brp),], "matrix"))
head(sig.1conv.brp)
dim(sig.1conv.brp)
#166 diff abund taxa

#Div. bulk vs Div. rhizoplane
result.1div.brp = results(test, contrast=c("group", "TP1Div.rhizoplane", "TP1Div.bulk"))
mcols(result.1div.brp, use.names=TRUE)
#Extract only significant results
sig.1div.brp = result.1div.brp[which(result.1div.brp$padj < alpha),]
sig.1div.brp = cbind(as(sig.1div.brp, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.1div.brp),], "matrix"))
head(sig.1div.brp)
dim(sig.1div.brp)
#137 diff abund taxa

#Conv. rhizoplane vs Div. rhizoplane
#I chose Conv. rhizoplane to be the baseline for this comparison
result.1convdiv.rp = results(test, contrast=c("group", "TP1Div.rhizoplane", "TP1Conv.rhizoplane"))
mcols(result.1convdiv.rp, use.names=TRUE)
#extract only significant results
sig.1convdiv.rp = result.1convdiv.rp[which(result.1convdiv.rp$padj < alpha),]
sig.1convdiv.rp = cbind(as(sig.1convdiv.rp, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.1convdiv.rp),], "matrix"))
head(sig.1convdiv.rp)
dim(sig.1convdiv.rp)
#17 diff abund taxa 

#But, I'm interested in taxa that responded to the root, many of these differences probably also exist in the bulk soil, meaning the root wasn't the cause of the difference
#I want to remove those
#The code below retains only those OTUs that are significant in the rhizoplane vs rhizoplane comparison AND bulk vs rhizoplane comparison
tp1.diff1.rp <- sig.1convdiv.rp[(rownames(sig.1convdiv.rp) %in% rownames(sig.1conv.brp)),]
head(tp1.diff1.rp)
dim(tp1.diff1.rp)
#7 diff abund taxa between Conv. bulk and Conv. rhizoplane
tp1.diff2.rp <- sig.1convdiv.rp[(rownames(sig.1convdiv.rp) %in% rownames(sig.1div.brp)),]
head(tp1.diff2.rp)
dim(tp1.diff2.rp)
#9 diff abund taxa between Div. bulk and Div. rhizoplane
#Find OTUs that are not shared in common between diff1 and diff2
tp1.diff3.rp <- tp1.diff2.rp[!(rownames(tp1.diff2.rp) %in% rownames(tp1.diff1.rp)),]
head(tp1.diff3.rp)
dim(tp1.diff3.rp)
#4 unique OTUs
#Add unique OTUs to diff1
tp1.diff4.rp <- rbind(tp1.diff1.rp, tp1.diff3.rp)
head(tp1.diff4.rp)
dim(tp1.diff4.rp)
#11 diff abund taxa between rhizoplanes that are also diff abund from the bulk soils

############
#Timepoint 2

#Conv. bulk vs Conv. rhizoplane
result.2conv.brp = results(test, contrast=c("group", "TP2Conv.rhizoplane", "TP2Conv.bulk"))
mcols(result.2conv.brp, use.names=TRUE)
sig.2conv.brp = result.2conv.brp[which(result.2conv.brp$padj < alpha),]
sig.2conv.brp = cbind(as(sig.2conv.brp, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.2conv.brp),], "matrix"))
head(sig.2conv.brp)
dim(sig.2conv.brp)
#174 diff abund taxa

#Div. bulk vs Div. rhizoplane
result.2div.brp = results(test, contrast=c("group", "TP2Div.rhizoplane", "TP2Div.bulk"))
mcols(result.2div.brp, use.names=TRUE)
#Extract only significant results
sig.2div.brp = result.2div.brp[which(result.2div.brp$padj < alpha),]
sig.2div.brp = cbind(as(sig.2div.brp, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.2div.brp),], "matrix"))
head(sig.2div.brp)
dim(sig.2div.brp)
#123 diff abund taxa

#Conv. rhizoplane vs Div. rhizoplane
#I chose Conv. rhizoplane to be the baseline for this comparison
result.2convdiv.rp = results(test, contrast=c("group", "TP2Div.rhizoplane", "TP2Conv.rhizoplane"))
mcols(result.2convdiv.rp, use.names=TRUE)
#extract only significant results
sig.2convdiv.rp = result.2convdiv.rp[which(result.2convdiv.rp$padj < alpha),]
sig.2convdiv.rp = cbind(as(sig.2convdiv.rp, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.2convdiv.rp),], "matrix"))
head(sig.2convdiv.rp)
dim(sig.2convdiv.rp)
#46 diff abund taxa 

#But, I'm interested in taxa that responded to the root, many of these differences probably also exist in the bulk soil, meaning the root wasn't the cause of the difference
#I want to remove those
#The code below retains only those OTUs that are significant in the rhizoplane vs rhizoplane comparison AND bulk vs rhizoplane comparison
tp2.diff1.rp <- sig.2convdiv.rp[(rownames(sig.2convdiv.rp) %in% rownames(sig.2conv.brp)),]
head(tp2.diff1.rp)
dim(tp2.diff1.rp)
#26 diff abund taxa between Conv. bulk and Conv. rhizoplane
tp2.diff2.rp <- sig.2convdiv.rp[(rownames(sig.2convdiv.rp) %in% rownames(sig.2div.brp)),]
head(tp2.diff2.rp)
dim(tp2.diff2.rp)
#19 diff abund taxa between Div. bulk and Div. rhizoplane
#Find OTUs that are not shared in common between diff1 and diff2
tp2.diff3.rp <- tp2.diff2.rp[!(rownames(tp2.diff2.rp) %in% rownames(tp2.diff1.rp)),]
head(tp2.diff3.rp)
dim(tp2.diff3.rp)
#4 unique OTUs
#Add unique OTUs to diff1
tp2.diff4.rp <- rbind(tp2.diff1.rp, tp2.diff3.rp)
head(tp2.diff4.rp)
dim(tp2.diff4.rp)
#30 diff abund taxa between rhizoplanes that are also diff abund from the bulk soils

############
#Timepoint 3

#Conv. bulk vs Conv. rhizoplane
result.3conv.brp = results(test, contrast=c("group", "TP3Conv.rhizoplane", "TP3Conv.bulk"))
mcols(result.3conv.brp, use.names=TRUE)
#Extract only significant results
sig.3conv.brp = result.3conv.brp[which(result.3conv.brp$padj < alpha),]
sig.3conv.brp = cbind(as(sig.3conv.brp, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.3conv.brp),], "matrix"))
head(sig.3conv.brp)
dim(sig.3conv.brp)
#157 diff abund taxa

#Div. bulk vs Div. rhizoplane
result.3div.brp = results(test, contrast=c("group", "TP3Div.rhizoplane", "TP3Div.bulk"))
mcols(result.3div.brp, use.names=TRUE)
#Extract only significant results
sig.3div.brp = result.3div.brp[which(result.3div.brp$padj < alpha),]
sig.3div.brp = cbind(as(sig.3div.brp, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.3div.brp),], "matrix"))
head(sig.3div.brp)
dim(sig.3div.brp)
#151 diff abund taxa

#Conv. rhizoplane vs Div. rhizoplane
#I chose Conv. rhizoplane to be the baseline for this comparison
result.3convdiv.rp = results(test, contrast=c("group", "TP3Div.rhizoplane", "TP3Conv.rhizoplane"))
mcols(result.3convdiv.rp, use.names=TRUE)
#extract only significant results
sig.3convdiv.rp = result.3convdiv.rp[which(result.3convdiv.rp$padj < alpha),]
sig.3convdiv.rp = cbind(as(sig.3convdiv.rp, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.3convdiv.rp),], "matrix"))
head(sig.3convdiv.rp)
dim(sig.3convdiv.rp)
#15 diff abund taxa 

#But, I'm interested in taxa that responded to the root, many of these differences probably also exist in the bulk soil, meaning the root wasn't the cause of the difference
#I want to remove those
#The code below retains only those OTUs that are significant in the rhizoplane vs rhizoplane comparison AND bulk vs rhizoplane comparison
tp3.diff1.rp <- sig.3convdiv.rp[(rownames(sig.3convdiv.rp) %in% rownames(sig.3conv.brp)),]
head(tp3.diff1.rp)
dim(tp3.diff1.rp)
#6 diff abund taxa between Conv. bulk and Conv. rhizoplane
tp3.diff2.rp <- sig.3convdiv.rp[(rownames(sig.3convdiv.rp) %in% rownames(sig.3div.brp)),]
head(tp3.diff2.rp)
dim(tp3.diff2.rp)
#7 diff abund taxa between Div. bulk and Div. rhizoplane
#Find OTUs that are not shared in common between diff1 and diff2
tp3.diff3.rp <- tp3.diff2.rp[!(rownames(tp3.diff2.rp) %in% rownames(tp3.diff1.rp)),]
head(tp3.diff3.rp)
dim(tp3.diff3.rp)
#2 unique OTUs
#Add unique OTUs to diff1
tp3.diff4.rp <- rbind(tp3.diff1.rp, tp3.diff3.rp)
head(tp3.diff4.rp)
dim(tp3.diff4.rp)
#8 diff abund taxa between rhizoplanes that are also diff abund from the bulk soils

############
#Timepoint 4

#Conv. bulk vs Conv. rhizoplane
result.4conv.brp = results(test, contrast=c("group", "TP4Conv.rhizoplane", "TP4Conv.bulk"))
mcols(result.4conv.brp, use.names=TRUE)
#Extract only significant results
sig.4conv.brp = result.4conv.brp[which(result.4conv.brp$padj < alpha),]
sig.4conv.brp = cbind(as(sig.4conv.brp, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.4conv.brp),], "matrix"))
head(sig.4conv.brp)
dim(sig.4conv.brp)
#163 diff abund taxa

#Div. bulk vs Div. rhizoplane
result.4div.brp = results(test, contrast=c("group", "TP4Div.rhizoplane", "TP4Div.bulk"))
mcols(result.4div.brp, use.names=TRUE)
#Extract only significant results
sig.4div.brp = result.4div.brp[which(result.4div.brp$padj < alpha),]
sig.4div.brp = cbind(as(sig.4div.brp, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.4div.brp),], "matrix"))
head(sig.4div.brp)
dim(sig.4div.brp)
#153 diff abund taxa

#Conv. rhizoplane vs Div. rhizoplane
#I chose Conv. rhizoplane to be the baseline for this comparison
result.4convdiv.rp = results(test, contrast=c("group", "TP4Div.rhizoplane", "TP4Conv.rhizoplane"))
mcols(result.4convdiv.rp, use.names=TRUE)
#extract only significant results
sig.4convdiv.rp = result.4convdiv.rp[which(result.4convdiv.rp$padj < alpha),]
sig.4convdiv.rp = cbind(as(sig.4convdiv.rp, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.4convdiv.rp),], "matrix"))
head(sig.4convdiv.rp)
dim(sig.4convdiv.rp)
#10 diff abund taxa 

#But, I'm interested in taxa that responded to the root, many of these differences probably also exist in the bulk soil, meaning the root wasn't the cause of the difference
#I want to remove those
#The code below retains only those OTUs that are significant in the rhizoplane vs rhizoplane comparison AND bulk vs rhizoplane comparison
tp4.diff1.rp <- sig.4convdiv.rp[(rownames(sig.4convdiv.rp) %in% rownames(sig.4conv.brp)),]
head(tp4.diff1.rp)
dim(tp4.diff1.rp)
#3 diff abund taxa between Conv. bulk and Conv. rhizoplane
tp4.diff2.rp <- sig.4convdiv.rp[(rownames(sig.4convdiv.rp) %in% rownames(sig.4div.brp)),]
head(tp4.diff2.rp)
dim(tp4.diff2.rp)
#4 diff abund taxa between Div. bulk and Div. rhizoplane
#Find OTUs that are not shared in common between diff1 and diff2
tp4.diff3.rp <- tp4.diff2.rp[!(rownames(tp4.diff2.rp) %in% rownames(tp4.diff1.rp)),]
head(tp4.diff3.rp)
dim(tp4.diff3.rp)
#1 unique OTUs
#Add unique OTUs to diff1
tp4.diff4.rp <- rbind(tp4.diff1.rp, tp4.diff3.rp)
head(tp4.diff4.rp)
dim(tp4.diff4.rp)
#4 diff abund taxa between rhizoplanes that are also diff abund from the bulk soils


#######################
#Rhizosphere Comparisons

############
#Timepoint 1

#Conv. bulk vs Conv. rhizosphere
#For contrasts, the treatment you want to act as the control should be the last listed
##here I set TP1Conv.bulk last because I want the bulk soil to be the baseline for comparison
result.1conv.brs = results(test, contrast=c("group", "TP1Conv.rhizosphere", "TP1Conv.bulk"))
mcols(result.1conv.brs, use.names=TRUE)
#Extract only significant results, at alpha level 0.01 (conservative)
sig.1conv.brs = result.1conv.brs[which(result.1conv.brs$padj < alpha),]
sig.1conv.brs = cbind(as(sig.1conv.brs, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.1conv.brs),], "matrix"))
head(sig.1conv.brs)
dim(sig.1conv.brs)
#4 diff abund taxa

#Div. bulk vs Div. rhizosphere
result.1div.brs = results(test, contrast=c("group", "TP1Div.rhizosphere", "TP1Div.bulk"))
mcols(result.1div.brs, use.names=TRUE)
#Extract only significant results
sig.1div.brs = result.1div.brs[which(result.1div.brs$padj < alpha),]
sig.1div.brs = cbind(as(sig.1div.brs, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.1div.brs),], "matrix"))
head(sig.1div.brs)
dim(sig.1div.brs)
#1 diff abund taxa

#Conv. rhizosphere vs Div. rhizosphere
#I chose Conv. rhizosphere to be the baseline for this comparison
result.1convdiv.rs = results(test, contrast=c("group", "TP1Div.rhizosphere", "TP1Conv.rhizosphere"))
mcols(result.1convdiv.rs, use.names=TRUE)
#extract only significant results
sig.1convdiv.rs = result.1convdiv.rs[which(result.1convdiv.rs$padj < alpha),]
sig.1convdiv.rs = cbind(as(sig.1convdiv.rs, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.1convdiv.rs),], "matrix"))
head(sig.1convdiv.rs)
dim(sig.1convdiv.rs)
#33 diff abund taxa 

#But, I'm interested in taxa that responded to the root, many of these differences probably also exist in the bulk soil, meaning the root wasn't the cause of the difference
#I want to remove those
#The code below retains only those families that are significant in the rhizosphere vs rhizosphere comparison AND bulk vs rhizosphere comparison
tp1.diff1.rs <- sig.1convdiv.rs[(rownames(sig.1convdiv.rs) %in% rownames(sig.1conv.brs)),]
head(tp1.diff1.rs)
dim(tp1.diff1.rs)
#0 diff abund taxa between Conv. bulk and Conv. rhizosphere
tp1.diff2.rs <- sig.1convdiv.rs[(rownames(sig.1convdiv.rs) %in% rownames(sig.1div.brs)),]
head(tp1.diff2.rs)
dim(tp1.diff2.rs)
#1 diff abund taxa between Div. bulk and Div. rhizosphere
#Find OTUs that are not shared in common between diff1 and diff2
tp1.diff3.rs <- tp1.diff2.rs[!(rownames(tp1.diff2.rs) %in% rownames(tp1.diff1.rs)),]
head(tp1.diff3.rs)
dim(tp1.diff3.rs)
#1 unique OTUs
#Add unique OTUs to diff1
tp1.diff4.rs <- rbind(tp1.diff1.rs, tp1.diff3.rs)
head(tp1.diff4.rs)
dim(tp1.diff4.rs)
#1 diff abund taxa between rhizospheres that are also diff abund from the bulk soils

############
#Timepoint 2

#Conv. bulk vs Conv. rhizosphere
result.2conv.brs = results(test, contrast=c("group", "TP2Conv.rhizosphere", "TP2Conv.bulk"))
mcols(result.2conv.brs, use.names=TRUE)
sig.2conv.brs = result.2conv.brs[which(result.2conv.brs$padj < alpha),]
sig.2conv.brs = cbind(as(sig.2conv.brs, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.2conv.brs),], "matrix"))
head(sig.2conv.brs)
dim(sig.2conv.brs)
#26 diff abund taxa

#Div. bulk vs Div. rhizosphere
result.2div.brs = results(test, contrast=c("group", "TP2Div.rhizosphere", "TP2Div.bulk"))
mcols(result.2div.brs, use.names=TRUE)
#Extract only significant results
sig.2div.brs = result.2div.brs[which(result.2div.brs$padj < alpha),]
sig.2div.brs = cbind(as(sig.2div.brs, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.2div.brs),], "matrix"))
head(sig.2div.brs)
dim(sig.2div.brs)
#4 diff abund taxa

#Conv. rhizosphere vs Div. rhizosphere
#I chose Conv. rhizosphere to be the baseline for this comparison
result.2convdiv.rs = results(test, contrast=c("group", "TP2Div.rhizosphere", "TP2Conv.rhizosphere"))
mcols(result.2convdiv.rs, use.names=TRUE)
#extract only significant results
sig.2convdiv.rs = result.2convdiv.rs[which(result.2convdiv.rs$padj < alpha),]
sig.2convdiv.rs = cbind(as(sig.2convdiv.rs, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.2convdiv.rs),], "matrix"))
head(sig.2convdiv.rs)
dim(sig.2convdiv.rs)
#33 diff abund taxa 

#But, I'm interested in taxa that responded to the root, many of these differences probably also exist in the bulk soil, meaning the root wasn't the cause of the difference
#I want to remove those
#The code below retains only those families that are significant in the rhizosphere vs rhizosphere comparison AND bulk vs rhizosphere comparison
tp2.diff1.rs <- sig.2convdiv.rs[(rownames(sig.2convdiv.rs) %in% rownames(sig.2conv.brs)),]
head(tp2.diff1.rs)
dim(tp2.diff1.rs)
#1 diff abund taxa between Conv. bulk and Conv. rhizosphere
tp2.diff2.rs <- sig.2convdiv.rs[(rownames(sig.2convdiv.rs) %in% rownames(sig.2div.brs)),]
head(tp2.diff2.rs)
dim(tp2.diff2.rs)
#0 diff abund taxa between Div. bulk and Div. rhizosphere
#Find OTUs that are not shared in common between diff1 and diff2
tp2.diff3.rs <- tp2.diff2.rs[!(rownames(tp2.diff2.rs) %in% rownames(tp2.diff1.rs)),]
head(tp2.diff3.rs)
dim(tp2.diff3.rs)
#0 unique OTUs
#Add unique OTUs to diff1
tp2.diff4.rs <- rbind(tp2.diff1.rs, tp2.diff3.rs)
head(tp2.diff4.rs)
dim(tp2.diff4.rs)
#1 diff abund taxa between rhizospheres that are also diff abund from the bulk soils

############
#Timepoint 3

#Conv. bulk vs Conv. rhizosphere
result.3conv.brs = results(test, contrast=c("group", "TP3Conv.rhizosphere", "TP3Conv.bulk"))
mcols(result.3conv.brs, use.names=TRUE)
#Extract only significant results
sig.3conv.brs = result.3conv.brs[which(result.3conv.brs$padj < alpha),]
sig.3conv.brs = cbind(as(sig.3conv.brs, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.3conv.brs),], "matrix"))
head(sig.3conv.brs)
dim(sig.3conv.brs)
#15 diff abund taxa

#Div. bulk vs Div. rhizosphere
result.3div.brs = results(test, contrast=c("group", "TP3Div.rhizosphere", "TP3Div.bulk"))
mcols(result.3div.brs, use.names=TRUE)
#Extract only significant results
sig.3div.brs = result.3div.brs[which(result.3div.brs$padj < alpha),]
sig.3div.brs = cbind(as(sig.3div.brs, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.3div.brs),], "matrix"))
head(sig.3div.brs)
dim(sig.3div.brs)
#6 diff abund taxa

#Conv. rhizosphere vs Div. rhizosphere
#I chose Conv. rhizosphere to be the baseline for this comparison
result.3convdiv.rs = results(test, contrast=c("group", "TP3Div.rhizosphere", "TP3Conv.rhizosphere"))
mcols(result.3convdiv.rs, use.names=TRUE)
#extract only significant results
sig.3convdiv.rs = result.3convdiv.rs[which(result.3convdiv.rs$padj < alpha),]
sig.3convdiv.rs = cbind(as(sig.3convdiv.rs, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.3convdiv.rs),], "matrix"))
head(sig.3convdiv.rs)
dim(sig.3convdiv.rs)
#22 diff abund taxa 

#But, I'm interested in taxa that responded to the root, many of these differences probably also exist in the bulk soil, meaning the root wasn't the cause of the difference
#I want to remove those
#The code below retains only those families that are significant in the rhizosphere vs rhizosphere comparison AND bulk vs rhizosphere comparison
tp3.diff1.rs <- sig.3convdiv.rs[(rownames(sig.3convdiv.rs) %in% rownames(sig.3conv.brs)),]
head(tp3.diff1.rs)
dim(tp3.diff1.rs)
#1 diff abund taxa between Conv. bulk and Conv. rhizosphere
tp3.diff2.rs <- sig.3convdiv.rs[(rownames(sig.3convdiv.rs) %in% rownames(sig.3div.brs)),]
head(tp3.diff2.rs)
dim(tp3.diff2.rs)
#0 diff abund taxa between Div. bulk and Div. rhizosphere
#Find OTUs that are not shared in common between diff1 and diff2
tp3.diff3.rs <- tp3.diff2.rs[!(rownames(tp3.diff2.rs) %in% rownames(tp3.diff1.rs)),]
head(tp3.diff3.rs)
dim(tp3.diff3.rs)
#0 unique OTUs
#Add unique OTUs to diff1
tp3.diff4.rs <- rbind(tp3.diff1.rs, tp3.diff3.rs)
head(tp3.diff4.rs)
dim(tp3.diff4.rs)
#1 diff abund taxa between rhizospheres that are also diff abund from the bulk soils

############
#Timepoint 4

#Conv. bulk vs Conv. rhizosphere
result.4conv.brs = results(test, contrast=c("group", "TP4Conv.rhizosphere", "TP4Conv.bulk"))
mcols(result.4conv.brs, use.names=TRUE)
#Extract only significant results
sig.4conv.brs = result.4conv.brs[which(result.4conv.brs$padj < alpha),]
sig.4conv.brs = cbind(as(sig.4conv.brs, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.4conv.brs),], "matrix"))
head(sig.4conv.brs)
dim(sig.4conv.brs)
#9 diff abund taxa

#Div. bulk vs Div. rhizosphere
result.4div.brs = results(test, contrast=c("group", "TP4Div.rhizosphere", "TP4Div.bulk"))
mcols(result.4div.brs, use.names=TRUE)
#Extract only significant results
sig.4div.brs = result.4div.brs[which(result.4div.brs$padj < alpha),]
sig.4div.brs = cbind(as(sig.4div.brs, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.4div.brs),], "matrix"))
head(sig.4div.brs)
dim(sig.4div.brs)
#3 diff abund taxa

#Conv. rhizosphere vs Div. rhizosphere
#I chose Conv. rhizosphere to be the baseline for this comparison
result.4convdiv.rs = results(test, contrast=c("group", "TP4Div.rhizosphere", "TP4Conv.rhizosphere"))
mcols(result.4convdiv.rs, use.names=TRUE)
#extract only significant results
sig.4convdiv.rs = result.4convdiv.rs[which(result.4convdiv.rs$padj < alpha),]
sig.4convdiv.rs = cbind(as(sig.4convdiv.rs, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.4convdiv.rs),], "matrix"))
head(sig.4convdiv.rs)
dim(sig.4convdiv.rs)
#9 diff abund taxa 

#But, I'm interested in taxa that responded to the root, many of these differences probably also exist in the bulk soil, meaning the root wasn't the cause of the difference
#I want to remove those
#The code below retains only those families that are significant in the rhizosphere vs rhizosphere comparison AND bulk vs rhizosphere comparison
tp4.diff1.rs <- sig.4convdiv.rs[(rownames(sig.4convdiv.rs) %in% rownames(sig.4conv.brs)),]
head(tp4.diff1.rs)
dim(tp4.diff1.rs)
#0 diff abund taxa between Conv. bulk and Conv. rhizosphere
tp4.diff2.rs <- sig.4convdiv.rs[(rownames(sig.4convdiv.rs) %in% rownames(sig.4div.brs)),]
head(tp4.diff2.rs)
dim(tp4.diff2.rs)
#0 diff abund taxa between Div. bulk and Div. rhizosphere
#Find OTUs that are not shared in common between diff1 and diff2
tp4.diff3.rs <- tp4.diff2.rs[!(rownames(tp4.diff2.rs) %in% rownames(tp4.diff1.rs)),]
head(tp4.diff3.rs)
dim(tp4.diff3.rs)
#0 unique OTUs
#Add unique OTUs to diff1
tp4.diff4.rs <- rbind(tp4.diff1.rs, tp4.diff3.rs)
head(tp4.diff4.rs)
dim(tp4.diff4.rs)
#0 diff abund taxa between rhizospheres that are also diff abund from the bulk soils


#################
#Bulk comparisons

############
#Timepoint 1

#Conv. Bulk vs Div. Bulk
result.1conv.b = results(test, contrast=c("group", "TP1Div.bulk", "TP1Conv.bulk"))
mcols(result.1conv.b, use.names=TRUE)
#extract only significant results
sig.1conv.b = result.1conv.b[which(result.1conv.b$padj < alpha),]
sig.1conv.b = cbind(as(sig.1conv.b, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.1conv.b),], "matrix"))
head(sig.1conv.b)
dim(sig.1conv.b)
#31 diff abund taxa beteween bulk soils

############
#Timepoint 2

#Conv. Bulk vs Div. Bulk
result.2conv.b = results(test, contrast=c("group", "TP2Div.bulk", "TP2Conv.bulk"))
mcols(result.2conv.b, use.names=TRUE)
#extract only significant results
sig.2conv.b = result.2conv.b[which(result.2conv.b$padj < alpha),]
sig.2conv.b = cbind(as(sig.2conv.b, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.2conv.b),], "matrix"))
head(sig.2conv.b)
dim(sig.2conv.b)
#33 diff abund taxa beteween bulk soils

############
#Timepoint 3

#Conv. Bulk vs Div. Bulk
result.3conv.b = results(test, contrast=c("group", "TP3Div.bulk", "TP3Conv.bulk"))
mcols(result.3conv.b, use.names=TRUE)
#extract only significant results
sig.3conv.b = result.3conv.b[which(result.3conv.b$padj < alpha),]
sig.3conv.b = cbind(as(sig.3conv.b, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.3conv.b),], "matrix"))
head(sig.3conv.b)
dim(sig.3conv.b)
#33 diff abund taxa beteween bulk soils

############
#Timepoint 4

#Conv. Bulk vs Div. Bulk
result.4conv.b = results(test, contrast=c("group", "TP4Div.bulk", "TP4Conv.bulk"))
mcols(result.4conv.b, use.names=TRUE)
#extract only significant results
sig.4conv.b = result.4conv.b[which(result.4conv.b$padj < alpha),]
sig.4conv.b = cbind(as(sig.4conv.b, "data.frame"), as(tax_table(deseq.fglom)[rownames(sig.4conv.b),], "matrix"))
head(sig.4conv.b)
dim(sig.4conv.b)
#27 diff abund taxa beteween bulk soils


###################
#Create figures 6 A

###########################
#Figure 6 data frame create

#add timepoint metadata
tp1.diff4.rp$Timepoint <- rep("TP1", nrow(tp1.diff4.rp))
tp2.diff4.rp$Timepoint <- rep("TP2", nrow(tp2.diff4.rp))
tp3.diff4.rp$Timepoint <- rep("TP3", nrow(tp3.diff4.rp))
tp4.diff4.rp$Timepoint <- rep("TP4", nrow(tp4.diff4.rp))

#Create OTU column for each
#I've found that the rownames change when combined (because of duplicates), meaning the OTU identity is no longer reliable if you don't include this step
tp1.diff4.rp$OTU <- row.names(tp1.diff4.rp)
tp2.diff4.rp$OTU <- row.names(tp2.diff4.rp)
tp3.diff4.rp$OTU <- row.names(tp3.diff4.rp)
tp4.diff4.rp$OTU <- row.names(tp4.diff4.rp)

#Combine all rp vs rp data
rp.all <- rbind(tp1.diff4.rp, tp2.diff4.rp, tp3.diff4.rp, tp4.diff4.rp)
head(rp.all)

#Create a column to indicate if log2fold change is greater or less than 0
#color bars based on direction of change
rp.colors <- mutate(rp.all, color=log2FoldChange>0)
head(rp.colors)

#How many unknown families are there?
length(which(rp.colors$Family=="f__?"))
#11 unknown families, need to be separately labeled so they aren't conflated in the graph

#Make nice labels
rp.colors$Phylum <- gsub("p__", "", rp.colors$Phylum)
rp.colors$Class <- gsub("c__", "", rp.colors$Class)
rp.colors$Order <- gsub("o__", "", rp.colors$Order)
rp.colors$Family <- gsub("f__", "", rp.colors$Family)

#Make sure all families are graphed separately and not conflated
rp.alltax <- mutate(rp.colors, alltax = paste0(OTU, sep=",", Phylum, sep=",", Class, sep=",", Order, sep=",", Family))
head(rp.alltax)

#Unfortunately, I had to order these manually in excel to get the format of labelling that I wanted
y.taxlabs <- rev(c("Actinobacteria, Catenulisporaceae",
                   "Actinobacteria, 288-2",
                   "Armatimonadetes, Armatimonadales",
                   "Proteobacteria, Burkholderiaceae",
                   "Proteobacteria, Rickettsiales",
                   "Firmicutes, Erysipelotrichaceae",
                   "Proteobacteria, NKB5",
                   "Proteobacteria, Oligoflexales",
                   "Proteobacteria, Oxalobacteraceae",
                   "Actinobacteria, Rubrobacteriaceae",
                   "Proteobacteria, Bradyrhizobiaceae",
                   "Chloroflexi, C0119",
                   "Chloroflexi, Roseiflexaceae",
                   "Actinobacteria, Frankiaceae",
                   "Planctomycetes, WD2101 soil group",
                   "Proteobacteria, Sphingomonadaceae",
                   "Proteobacteria, mle1-27",
                   "Proteobacteria, Rhizobiales",
                   "Proteobacteria, Coxiellaceae",
                   "Actinobacteria, Solirubrobacteraceae",
                   "Thaumarchaeota, Soil Crenarchaeotic Group",
                   "Planctomycetes, Planctomycetaceae",
                   "Acidobacteria, Subgroup 6",
                   "Nitrospirae, Nitrospiraceae",
                   "Actinobacteria, Pseudonocardiaceae",
                   "Proteobacteria, GR-WP33-30",
                   "Proteobacteria, Oligoflexales",
                   "Acidobacteria, Subgroup 5",
                   "Verrucomicrobia, OPB35 soil group",
                   "Verrucomicrobia, Verrucomicrobiaceae",
                   "Chlorobi, OPB56",
                   "Proteobacteria, MNG7",
                   "Acidobacteria, Subgroup 11",
                   "Proteobacteria, TA18",
                   "Bacteroidetes, Saprospiraceae",
                   "Bacteroidetes, PHOS-HE51",
                   "Acidobacteria, Subgroup 18",
                   "Fibrobacteres, Fibrobacteraceae",
                   "Proteobacteria, Methylophilaceae",
                   "Actinobacteria, Glycomycetaceae",
                   "Acidobacteria, Subgroup 11",
                   "Bacteroidetes, Cyclobacteriaceae",
                   "Cyanobacteria, Vampirovibrionaceae",
                   "Proteobacteria, Pseudomonadales"))
#made sure these were in correct order

#graph figure 6
bar.rp <- ggplot(rp.alltax, aes(x=reorder(alltax, -log2FoldChange), y=log2FoldChange)) +
  geom_col(aes(fill=color)) +
  coord_flip() +
  labs(x="", y="log2Fold Change (relative to Conventional)", 
       title="") +
  scale_x_discrete(labels=y.taxlabs) +
  facet_grid(.~Timepoint) +
  geom_hline(yintercept=0, color="grey69", linetype=2) +
  theme_bw() +
  theme(axis.text.y = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill=FALSE)
bar.rp



