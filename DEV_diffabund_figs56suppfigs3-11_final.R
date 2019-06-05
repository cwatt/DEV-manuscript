
# Cassandra Wattenburger
# 01/25/18

# Cleaned up scripts to generate differential abundance data and figures 5, 6, and supp figures 3-11

############################################

library("phyloseq")
library("deseq2")
library("vegan")
library("ggplot2")
library("plyr")
library("ggplot2")
library("tidyr")
library("reshape2")
library("cowplot")

setwd("Y:/Cassi/Sequencing Data/Argonne 2017/01 DEV Project/03 Final Scripts and figures")


#########
#16s Data

rm(list=ls())

############
#DESeq2 test

##Load in phyloseq object for DESeq2
#Created in "DEV16S_deseq2_prep_clean.R"
load("DEV16S.deseq.want.Rdata")
deseq.want

#convert phyloseq object to deseq2 object
deseq = phyloseq_to_deseq2(deseq.want, ~ group)
deseq

test = DESeq(deseq, test="Wald", fitType="parametric")
#note: fitType='parametric', but the dispersion trend was not well captured by the
#function: y = a/x + b, and a local regression fit was automatically substituted.
#specify fitType='local' or 'mean' to avoid this message next time.

resultsNames(test)

#Set alpha level to 0.01
#This reduces family-wise error to 1%
alpha = 0.01


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
sig.1conv.brp = cbind(as(sig.1conv.brp, "data.frame"), as(tax_table(deseq.want)[rownames(sig.1conv.brp),], "matrix"))
head(sig.1conv.brp)
dim(sig.1conv.brp)
#705 diff abund taxa

#Div. bulk vs Div. rhizoplane
result.1div.brp = results(test, contrast=c("group", "TP1Div.rhizoplane", "TP1Div.bulk"))
mcols(result.1div.brp, use.names=TRUE)
#Extract only significant results
sig.1div.brp = result.1div.brp[which(result.1div.brp$padj < alpha),]
sig.1div.brp = cbind(as(sig.1div.brp, "data.frame"), as(tax_table(deseq.want)[rownames(sig.1div.brp),], "matrix"))
head(sig.1div.brp)
dim(sig.1div.brp)
#556 diff abund taxa

#Conv. rhizoplane vs Div. rhizoplane
#I chose Conv. rhizoplane to be the baseline for this comparison
result.1convdiv.rp = results(test, contrast=c("group", "TP1Div.rhizoplane", "TP1Conv.rhizoplane"))
mcols(result.1convdiv.rp, use.names=TRUE)
#extract only significant results
sig.1convdiv.rp = result.1convdiv.rp[which(result.1convdiv.rp$padj < alpha),]
sig.1convdiv.rp = cbind(as(sig.1convdiv.rp, "data.frame"), as(tax_table(deseq.want)[rownames(sig.1convdiv.rp),], "matrix"))
head(sig.1convdiv.rp)
dim(sig.1convdiv.rp)
#81 diff abund taxa 

#But, I'm interested in taxa that responded to the root, many of these differences probably also exist in the bulk soil, meaning the root wasn't the cause of the difference
#I want to remove those
#The code below retains only those OTUs that are significant in the rhizoplane vs rhizoplane comparison AND bulk vs rhizoplane comparison
tp1.diff1.rp <- sig.1convdiv.rp[(rownames(sig.1convdiv.rp) %in% rownames(sig.1conv.brp)),]
head(tp1.diff1.rp)
dim(tp1.diff1.rp)
#45 diff abund taxa between Conv. bulk and Conv. rhizoplane
tp1.diff2.rp <- sig.1convdiv.rp[(rownames(sig.1convdiv.rp) %in% rownames(sig.1div.brp)),]
head(tp1.diff2.rp)
dim(tp1.diff2.rp)
#36 diff abund taxa between Div. bulk and Div. rhizoplane
#Find OTUs that are not shared in common between diff1 and diff2
tp1.diff3.rp <- tp1.diff2.rp[!(rownames(tp1.diff2.rp) %in% rownames(tp1.diff1.rp)),]
head(tp1.diff3.rp)
dim(tp1.diff3.rp)
#8 unique OTUs
#Add unique OTUs to diff1
tp1.diff4.rp <- rbind(tp1.diff1.rp, tp1.diff3.rp)
head(tp1.diff4.rp)
dim(tp1.diff4.rp)
#53 diff abund taxa between rhizoplanes that are also diff abund from the bulk soils

############
#Timepoint 2

#Conv. bulk vs Conv. rhizoplane
result.2conv.brp = results(test, contrast=c("group", "TP2Conv.rhizoplane", "TP2Conv.bulk"))
mcols(result.2conv.brp, use.names=TRUE)
sig.2conv.brp = result.2conv.brp[which(result.2conv.brp$padj < alpha),]
sig.2conv.brp = cbind(as(sig.2conv.brp, "data.frame"), as(tax_table(deseq.want)[rownames(sig.2conv.brp),], "matrix"))
head(sig.2conv.brp)
dim(sig.2conv.brp)
#708 diff abund taxa

#Div. bulk vs Div. rhizoplane
result.2div.brp = results(test, contrast=c("group", "TP2Div.rhizoplane", "TP2Div.bulk"))
mcols(result.2div.brp, use.names=TRUE)
#Extract only significant results
sig.2div.brp = result.2div.brp[which(result.2div.brp$padj < alpha),]
sig.2div.brp = cbind(as(sig.2div.brp, "data.frame"), as(tax_table(deseq.want)[rownames(sig.2div.brp),], "matrix"))
head(sig.2div.brp)
dim(sig.2div.brp)
#578 diff abund taxa

#Conv. rhizoplane vs Div. rhizoplane
#I chose Conv. rhizoplane to be the baseline for this comparison
result.2convdiv.rp = results(test, contrast=c("group", "TP2Div.rhizoplane", "TP2Conv.rhizoplane"))
mcols(result.2convdiv.rp, use.names=TRUE)
#extract only significant results
sig.2convdiv.rp = result.2convdiv.rp[which(result.2convdiv.rp$padj < alpha),]
sig.2convdiv.rp = cbind(as(sig.2convdiv.rp, "data.frame"), as(tax_table(deseq.want)[rownames(sig.2convdiv.rp),], "matrix"))
head(sig.2convdiv.rp)
dim(sig.2convdiv.rp)
#125 diff abund taxa 

#But, I'm interested in taxa that responded to the root, many of these differences probably also exist in the bulk soil, meaning the root wasn't the cause of the difference
#I want to remove those
#The code below retains only those OTUs that are significant in the rhizoplane vs rhizoplane comparison AND bulk vs rhizoplane comparison
tp2.diff1.rp <- sig.2convdiv.rp[(rownames(sig.2convdiv.rp) %in% rownames(sig.2conv.brp)),]
head(tp2.diff1.rp)
dim(tp2.diff1.rp)
#69 diff abund taxa between Conv. bulk and Conv. rhizoplane
tp2.diff2.rp <- sig.2convdiv.rp[(rownames(sig.2convdiv.rp) %in% rownames(sig.2div.brp)),]
head(tp2.diff2.rp)
dim(tp2.diff2.rp)
#31 diff abund taxa between Div. bulk and Div. rhizoplane
#Find OTUs that are not shared in common between diff1 and diff2
tp2.diff3.rp <- tp2.diff2.rp[!(rownames(tp2.diff2.rp) %in% rownames(tp2.diff1.rp)),]
head(tp2.diff3.rp)
dim(tp2.diff3.rp)
#9 unique OTUs
#Add unique OTUs to diff1
tp2.diff4.rp <- rbind(tp2.diff1.rp, tp2.diff3.rp)
head(tp2.diff4.rp)
dim(tp2.diff4.rp)
#78 diff abund taxa between rhizoplanes that are also diff abund from the bulk soils

############
#Timepoint 3

#Conv. bulk vs Conv. rhizoplane
result.3conv.brp = results(test, contrast=c("group", "TP3Conv.rhizoplane", "TP3Conv.bulk"))
mcols(result.3conv.brp, use.names=TRUE)
#Extract only significant results
sig.3conv.brp = result.3conv.brp[which(result.3conv.brp$padj < alpha),]
sig.3conv.brp = cbind(as(sig.3conv.brp, "data.frame"), as(tax_table(deseq.want)[rownames(sig.3conv.brp),], "matrix"))
head(sig.3conv.brp)
dim(sig.3conv.brp)
#674 diff abund taxa

#Div. bulk vs Div. rhizoplane
result.3div.brp = results(test, contrast=c("group", "TP3Div.rhizoplane", "TP3Div.bulk"))
mcols(result.3div.brp, use.names=TRUE)
#Extract only significant results
sig.3div.brp = result.3div.brp[which(result.3div.brp$padj < alpha),]
sig.3div.brp = cbind(as(sig.3div.brp, "data.frame"), as(tax_table(deseq.want)[rownames(sig.3div.brp),], "matrix"))
head(sig.3div.brp)
dim(sig.3div.brp)
#685 diff abund taxa

#Conv. rhizoplane vs Div. rhizoplane
#I chose Conv. rhizoplane to be the baseline for this comparison
result.3convdiv.rp = results(test, contrast=c("group", "TP3Div.rhizoplane", "TP3Conv.rhizoplane"))
mcols(result.3convdiv.rp, use.names=TRUE)
#extract only significant results
sig.3convdiv.rp = result.3convdiv.rp[which(result.3convdiv.rp$padj < alpha),]
sig.3convdiv.rp = cbind(as(sig.3convdiv.rp, "data.frame"), as(tax_table(deseq.want)[rownames(sig.3convdiv.rp),], "matrix"))
head(sig.3convdiv.rp)
dim(sig.3convdiv.rp)
#61 diff abund taxa 

#But, I'm interested in taxa that responded to the root, many of these differences probably also exist in the bulk soil, meaning the root wasn't the cause of the difference
#I want to remove those
#The code below retains only those OTUs that are significant in the rhizoplane vs rhizoplane comparison AND bulk vs rhizoplane comparison
tp3.diff1.rp <- sig.3convdiv.rp[(rownames(sig.3convdiv.rp) %in% rownames(sig.3conv.brp)),]
head(tp3.diff1.rp)
dim(tp3.diff1.rp)
#26 diff abund taxa between Conv. bulk and Conv. rhizoplane
tp3.diff2.rp <- sig.3convdiv.rp[(rownames(sig.3convdiv.rp) %in% rownames(sig.3div.brp)),]
head(tp3.diff2.rp)
dim(tp3.diff2.rp)
#28 diff abund taxa between Div. bulk and Div. rhizoplane
#Find OTUs that are not shared in common between diff1 and diff2
tp3.diff3.rp <- tp3.diff2.rp[!(rownames(tp3.diff2.rp) %in% rownames(tp3.diff1.rp)),]
head(tp3.diff3.rp)
dim(tp3.diff3.rp)
#11 unique OTUs
#Add unique OTUs to diff1
tp3.diff4.rp <- rbind(tp3.diff1.rp, tp3.diff3.rp)
head(tp3.diff4.rp)
dim(tp3.diff4.rp)
#37 diff abund taxa between rhizoplanes that are also diff abund from the bulk soils

############
#Timepoint 4

#Conv. bulk vs Conv. rhizoplane
result.4conv.brp = results(test, contrast=c("group", "TP4Conv.rhizoplane", "TP4Conv.bulk"))
mcols(result.4conv.brp, use.names=TRUE)
#Extract only significant results
sig.4conv.brp = result.4conv.brp[which(result.4conv.brp$padj < alpha),]
sig.4conv.brp = cbind(as(sig.4conv.brp, "data.frame"), as(tax_table(deseq.want)[rownames(sig.4conv.brp),], "matrix"))
head(sig.4conv.brp)
dim(sig.4conv.brp)
#683 diff abund taxa

#Div. bulk vs Div. rhizoplane
result.4div.brp = results(test, contrast=c("group", "TP4Div.rhizoplane", "TP4Div.bulk"))
mcols(result.4div.brp, use.names=TRUE)
#Extract only significant results
sig.4div.brp = result.4div.brp[which(result.4div.brp$padj < alpha),]
sig.4div.brp = cbind(as(sig.4div.brp, "data.frame"), as(tax_table(deseq.want)[rownames(sig.4div.brp),], "matrix"))
head(sig.4div.brp)
dim(sig.4div.brp)
#664 diff abund taxa

#Conv. rhizoplane vs Div. rhizoplane
#I chose Conv. rhizoplane to be the baseline for this comparison
result.4convdiv.rp = results(test, contrast=c("group", "TP4Div.rhizoplane", "TP4Conv.rhizoplane"))
mcols(result.4convdiv.rp, use.names=TRUE)
#extract only significant results
sig.4convdiv.rp = result.4convdiv.rp[which(result.4convdiv.rp$padj < alpha),]
sig.4convdiv.rp = cbind(as(sig.4convdiv.rp, "data.frame"), as(tax_table(deseq.want)[rownames(sig.4convdiv.rp),], "matrix"))
head(sig.4convdiv.rp)
dim(sig.4convdiv.rp)
#39 diff abund taxa 

#But, I'm interested in taxa that responded to the root, many of these differences probably also exist in the bulk soil, meaning the root wasn't the cause of the difference
#I want to remove those
#The code below retains only those OTUs that are significant in the rhizoplane vs rhizoplane comparison AND bulk vs rhizoplane comparison
tp4.diff1.rp <- sig.4convdiv.rp[(rownames(sig.4convdiv.rp) %in% rownames(sig.4conv.brp)),]
head(tp4.diff1.rp)
dim(tp4.diff1.rp)
#15 diff abund taxa between Conv. bulk and Conv. rhizoplane
tp4.diff2.rp <- sig.4convdiv.rp[(rownames(sig.4convdiv.rp) %in% rownames(sig.4div.brp)),]
head(tp4.diff2.rp)
dim(tp4.diff2.rp)
#15 diff abund taxa between Div. bulk and Div. rhizoplane
#Find OTUs that are not shared in common between diff1 and diff2
tp4.diff3.rp <- tp4.diff2.rp[!(rownames(tp4.diff2.rp) %in% rownames(tp4.diff1.rp)),]
head(tp4.diff3.rp)
dim(tp4.diff3.rp)
#4 unique OTUs
#Add unique OTUs to diff1
tp4.diff4.rp <- rbind(tp4.diff1.rp, tp4.diff3.rp)
head(tp4.diff4.rp)
dim(tp4.diff4.rp)
#19 diff abund taxa between rhizoplanes that are also diff abund from the bulk soils

########################
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
sig.1conv.brs = cbind(as(sig.1conv.brs, "data.frame"), as(tax_table(deseq.want)[rownames(sig.1conv.brs),], "matrix"))
head(sig.1conv.brs)
dim(sig.1conv.brs)
#17 diff abund taxa

#Div. bulk vs Div. rhizosphere
result.1div.brs = results(test, contrast=c("group", "TP1Div.rhizosphere", "TP1Div.bulk"))
mcols(result.1div.brs, use.names=TRUE)
#Extract only significant results
sig.1div.brs = result.1div.brs[which(result.1div.brs$padj < alpha),]
sig.1div.brs = cbind(as(sig.1div.brs, "data.frame"), as(tax_table(deseq.want)[rownames(sig.1div.brs),], "matrix"))
head(sig.1div.brs)
dim(sig.1div.brs)
#11 diff abund taxa

#Conv. rhizosphere vs Div. rhizosphere
#I chose Conv. rhizosphere to be the baseline for this comparison
result.1convdiv.rs = results(test, contrast=c("group", "TP1Div.rhizosphere", "TP1Conv.rhizosphere"))
mcols(result.1convdiv.rs, use.names=TRUE)
#extract only significant results
sig.1convdiv.rs = result.1convdiv.rs[which(result.1convdiv.rs$padj < alpha),]
sig.1convdiv.rs = cbind(as(sig.1convdiv.rs, "data.frame"), as(tax_table(deseq.want)[rownames(sig.1convdiv.rs),], "matrix"))
head(sig.1convdiv.rs)
dim(sig.1convdiv.rs)
#194 diff abund taxa 

#But, I'm interested in taxa that responded to the root, many of these differences probably also exist in the bulk soil, meaning the root wasn't the cause of the difference
#I want to remove those
#The code below retains only those OTUs that are significant in the rhizosphere vs rhizosphere comparison AND bulk vs rhizosphere comparison
tp1.diff1.rs <- sig.1convdiv.rs[(rownames(sig.1convdiv.rs) %in% rownames(sig.1conv.brs)),]
head(tp1.diff1.rs)
dim(tp1.diff1.rs)
#5 diff abund taxa between Conv. bulk and Conv. rhizosphere
tp1.diff2.rs <- sig.1convdiv.rs[(rownames(sig.1convdiv.rs) %in% rownames(sig.1div.brs)),]
head(tp1.diff2.rs)
dim(tp1.diff2.rs)
#4 diff abund taxa between Div. bulk and Div. rhizosphere
#Find OTUs that are not shared in common between diff1 and diff2
tp1.diff3.rs <- tp1.diff2.rs[!(rownames(tp1.diff2.rs) %in% rownames(tp1.diff1.rs)),]
head(tp1.diff3.rs)
dim(tp1.diff3.rs)
#4 unique OTUs
#Add unique OTUs to diff1
tp1.diff4.rs <- rbind(tp1.diff1.rs, tp1.diff3.rs)
head(tp1.diff4.rs)
dim(tp1.diff4.rs)
#9 diff abund taxa between rhizospheres that are also diff abund from the bulk soils

############
#Timepoint 2

#Conv. bulk vs Conv. rhizosphere
result.2conv.brs = results(test, contrast=c("group", "TP2Conv.rhizosphere", "TP2Conv.bulk"))
mcols(result.2conv.brs, use.names=TRUE)
sig.2conv.brs = result.2conv.brs[which(result.2conv.brs$padj < alpha),]
sig.2conv.brs = cbind(as(sig.2conv.brs, "data.frame"), as(tax_table(deseq.want)[rownames(sig.2conv.brs),], "matrix"))
head(sig.2conv.brs)
dim(sig.2conv.brs)
#80 diff abund taxa

#Div. bulk vs Div. rhizosphere
result.2div.brs = results(test, contrast=c("group", "TP2Div.rhizosphere", "TP2Div.bulk"))
mcols(result.2div.brs, use.names=TRUE)
#Extract only significant results
sig.2div.brs = result.2div.brs[which(result.2div.brs$padj < alpha),]
sig.2div.brs = cbind(as(sig.2div.brs, "data.frame"), as(tax_table(deseq.want)[rownames(sig.2div.brs),], "matrix"))
head(sig.2div.brs)
dim(sig.2div.brs)
#26 diff abund taxa

#Conv. rhizosphere vs Div. rhizosphere
#I chose Conv. rhizosphere to be the baseline for this comparison
result.2convdiv.rs = results(test, contrast=c("group", "TP2Div.rhizosphere", "TP2Conv.rhizosphere"))
mcols(result.2convdiv.rs, use.names=TRUE)
#extract only significant results
sig.2convdiv.rs = result.2convdiv.rs[which(result.2convdiv.rs$padj < alpha),]
sig.2convdiv.rs = cbind(as(sig.2convdiv.rs, "data.frame"), as(tax_table(deseq.want)[rownames(sig.2convdiv.rs),], "matrix"))
head(sig.2convdiv.rs)
dim(sig.2convdiv.rs)
#162 diff abund taxa 

#But, I'm interested in taxa that responded to the root, many of these differences probably also exist in the bulk soil, meaning the root wasn't the cause of the difference
#I want to remove those
#The code below retains only those OTUs that are significant in the rhizosphere vs rhizosphere comparison AND bulk vs rhizosphere comparison
tp2.diff1.rs <- sig.2convdiv.rs[(rownames(sig.2convdiv.rs) %in% rownames(sig.2conv.brs)),]
head(tp2.diff1.rs)
dim(tp2.diff1.rs)
#16 diff abund taxa between Conv. bulk and Conv. rhizosphere
tp2.diff2.rs <- sig.2convdiv.rs[(rownames(sig.2convdiv.rs) %in% rownames(sig.2div.brs)),]
head(tp2.diff2.rs)
dim(tp2.diff2.rs)
#2 diff abund taxa between Div. bulk and Div. rhizosphere
#Find OTUs that are not shared in common between diff1 and diff2
tp2.diff3.rs <- tp2.diff2.rs[!(rownames(tp2.diff2.rs) %in% rownames(tp2.diff1.rs)),]
head(tp2.diff3.rs)
dim(tp2.diff3.rs)
#1 unique OTUs
#Add unique OTUs to diff1
tp2.diff4.rs <- rbind(tp2.diff1.rs, tp2.diff3.rs)
head(tp2.diff4.rs)
dim(tp2.diff4.rs)
#17 diff abund taxa between rhizospheres that are also diff abund from the bulk soils

############
#Timepoint 3

#Conv. bulk vs Conv. rhizosphere
result.3conv.brs = results(test, contrast=c("group", "TP3Conv.rhizosphere", "TP3Conv.bulk"))
mcols(result.3conv.brs, use.names=TRUE)
#Extract only significant results
sig.3conv.brs = result.3conv.brs[which(result.3conv.brs$padj < alpha),]
sig.3conv.brs = cbind(as(sig.3conv.brs, "data.frame"), as(tax_table(deseq.want)[rownames(sig.3conv.brs),], "matrix"))
head(sig.3conv.brs)
dim(sig.3conv.brs)
#56 diff abund taxa

#Div. bulk vs Div. rhizosphere
result.3div.brs = results(test, contrast=c("group", "TP3Div.rhizosphere", "TP3Div.bulk"))
mcols(result.3div.brs, use.names=TRUE)
#Extract only significant results
sig.3div.brs = result.3div.brs[which(result.3div.brs$padj < alpha),]
sig.3div.brs = cbind(as(sig.3div.brs, "data.frame"), as(tax_table(deseq.want)[rownames(sig.3div.brs),], "matrix"))
head(sig.3div.brs)
dim(sig.3div.brs)
#30 diff abund taxa

#Conv. rhizosphere vs Div. rhizosphere
#I chose Conv. rhizosphere to be the baseline for this comparison
result.3convdiv.rs = results(test, contrast=c("group", "TP3Div.rhizosphere", "TP3Conv.rhizosphere"))
mcols(result.3convdiv.rs, use.names=TRUE)
#extract only significant results
sig.3convdiv.rs = result.3convdiv.rs[which(result.3convdiv.rs$padj < alpha),]
sig.3convdiv.rs = cbind(as(sig.3convdiv.rs, "data.frame"), as(tax_table(deseq.want)[rownames(sig.3convdiv.rs),], "matrix"))
head(sig.3convdiv.rs)
dim(sig.3convdiv.rs)
#72 diff abund taxa 

#But, I'm interested in taxa that responded to the root, many of these differences probably also exist in the bulk soil, meaning the root wasn't the cause of the difference
#I want to remove those
#The code below retains only those OTUs that are significant in the rhizosphere vs rhizosphere comparison AND bulk vs rhizosphere comparison
tp3.diff1.rs <- sig.3convdiv.rs[(rownames(sig.3convdiv.rs) %in% rownames(sig.3conv.brs)),]
head(tp3.diff1.rs)
dim(tp3.diff1.rs)
#3 diff abund taxa between Conv. bulk and Conv. rhizosphere
tp3.diff2.rs <- sig.3convdiv.rs[(rownames(sig.3convdiv.rs) %in% rownames(sig.3div.brs)),]
head(tp3.diff2.rs)
dim(tp3.diff2.rs)
#2 diff abund taxa between Div. bulk and Div. rhizosphere
#Find OTUs that are not shared in common between diff1 and diff2
tp3.diff3.rs <- tp3.diff2.rs[!(rownames(tp3.diff2.rs) %in% rownames(tp3.diff1.rs)),]
head(tp3.diff3.rs)
dim(tp3.diff3.rs)
#1 unique OTUs
#Add unique OTUs to diff1
tp3.diff4.rs <- rbind(tp3.diff1.rs, tp3.diff3.rs)
head(tp3.diff4.rs)
dim(tp3.diff4.rs)
#4 diff abund taxa between rhizospheres that are also diff abund from the bulk soils

############
#Timepoint 4

#Conv. bulk vs Conv. rhizosphere
result.4conv.brs = results(test, contrast=c("group", "TP4Conv.rhizosphere", "TP4Conv.bulk"))
mcols(result.4conv.brs, use.names=TRUE)
#Extract only significant results
sig.4conv.brs = result.4conv.brs[which(result.4conv.brs$padj < alpha),]
sig.4conv.brs = cbind(as(sig.4conv.brs, "data.frame"), as(tax_table(deseq.want)[rownames(sig.4conv.brs),], "matrix"))
head(sig.4conv.brs)
dim(sig.4conv.brs)
#50 diff abund taxa

#Div. bulk vs Div. rhizosphere
result.4div.brs = results(test, contrast=c("group", "TP4Div.rhizosphere", "TP4Div.bulk"))
mcols(result.4div.brs, use.names=TRUE)
#Extract only significant results
sig.4div.brs = result.4div.brs[which(result.4div.brs$padj < alpha),]
sig.4div.brs = cbind(as(sig.4div.brs, "data.frame"), as(tax_table(deseq.want)[rownames(sig.4div.brs),], "matrix"))
head(sig.4div.brs)
dim(sig.4div.brs)
#29 diff abund taxa

#Conv. rhizosphere vs Div. rhizosphere
#I chose Conv. rhizosphere to be the baseline for this comparison
result.4convdiv.rs = results(test, contrast=c("group", "TP4Div.rhizosphere", "TP4Conv.rhizosphere"))
mcols(result.4convdiv.rs, use.names=TRUE)
#extract only significant results
sig.4convdiv.rs = result.4convdiv.rs[which(result.4convdiv.rs$padj < alpha),]
sig.4convdiv.rs = cbind(as(sig.4convdiv.rs, "data.frame"), as(tax_table(deseq.want)[rownames(sig.4convdiv.rs),], "matrix"))
head(sig.4convdiv.rs)
dim(sig.4convdiv.rs)
#38 diff abund taxa 

#But, I'm interested in taxa that responded to the root, many of these differences probably also exist in the bulk soil, meaning the root wasn't the cause of the difference
#I want to remove those
#The code below retains only those OTUs that are significant in the rhizosphere vs rhizosphere comparison AND bulk vs rhizosphere comparison
tp4.diff1.rs <- sig.4convdiv.rs[(rownames(sig.4convdiv.rs) %in% rownames(sig.4conv.brs)),]
head(tp4.diff1.rs)
dim(tp4.diff1.rs)
#2 diff abund taxa between Conv. bulk and Conv. rhizosphere
tp4.diff2.rs <- sig.4convdiv.rs[(rownames(sig.4convdiv.rs) %in% rownames(sig.4div.brs)),]
head(tp4.diff2.rs)
dim(tp4.diff2.rs)
#1 diff abund taxa between Div. bulk and Div. rhizosphere
#Find OTUs that are not shared in common between diff1 and diff2
tp4.diff3.rs <- tp4.diff2.rs[!(rownames(tp4.diff2.rs) %in% rownames(tp4.diff1.rs)),]
head(tp4.diff3.rs)
dim(tp4.diff3.rs)
#1 unique OTUs
#Add unique OTUs to diff1
tp4.diff4.rs <- rbind(tp4.diff1.rs, tp4.diff3.rs)
head(tp4.diff4.rs)
dim(tp4.diff4.rs)
#3 diff abund taxa between rhizospheres that are also diff abund from the bulk soils

#################
#Bulk comparisons

############
#Timepoint 1

#Conv. Bulk vs Div. Bulk
result.1conv.b = results(test, contrast=c("group", "TP1Div.bulk", "TP1Conv.bulk"))
mcols(result.1conv.b, use.names=TRUE)
#extract only significant results
sig.1conv.b = result.1conv.b[which(result.1conv.b$padj < alpha),]
sig.1conv.b = cbind(as(sig.1conv.b, "data.frame"), as(tax_table(deseq.want)[rownames(sig.1conv.b),], "matrix"))
head(sig.1conv.b)
dim(sig.1conv.b)
#121 diff abund taxa beteween bulk soils

############
#Timepoint 2

#Conv. Bulk vs Div. Bulk
result.2conv.b = results(test, contrast=c("group", "TP2Div.bulk", "TP2Conv.bulk"))
mcols(result.2conv.b, use.names=TRUE)
#extract only significant results
sig.2conv.b = result.2conv.b[which(result.2conv.b$padj < alpha),]
sig.2conv.b = cbind(as(sig.2conv.b, "data.frame"), as(tax_table(deseq.want)[rownames(sig.2conv.b),], "matrix"))
head(sig.2conv.b)
dim(sig.2conv.b)
#125 diff abund taxa beteween bulk soils

############
#Timepoint 3

#Conv. Bulk vs Div. Bulk
result.3conv.b = results(test, contrast=c("group", "TP3Div.bulk", "TP3Conv.bulk"))
mcols(result.3conv.b, use.names=TRUE)
#extract only significant results
sig.3conv.b = result.3conv.b[which(result.3conv.b$padj < alpha),]
sig.3conv.b = cbind(as(sig.3conv.b, "data.frame"), as(tax_table(deseq.want)[rownames(sig.3conv.b),], "matrix"))
head(sig.3conv.b)
dim(sig.3conv.b)
#164 diff abund taxa beteween bulk soils

############
#Timepoint 4

#Conv. Bulk vs Div. Bulk
result.4conv.b = results(test, contrast=c("group", "TP4Div.bulk", "TP4Conv.bulk"))
mcols(result.4conv.b, use.names=TRUE)
#extract only significant results
sig.4conv.b = result.4conv.b[which(result.4conv.b$padj < alpha),]
sig.4conv.b = cbind(as(sig.4conv.b, "data.frame"), as(tax_table(deseq.want)[rownames(sig.4conv.b),], "matrix"))
head(sig.4conv.b)
dim(sig.4conv.b)
#107 diff abund taxa beteween bulk soils

################
# Create figures

# Figure 6

###########################################################################################################################
load("DEV16S.tp1.conv.brs.RData")
load("DEV16S.tp1.div.brs.RData")
load("DEV16S.tp2.conv.brs.RData")
load("DEV16S.tp2.div.brs.RData")
load("DEV16S.tp3.conv.brs.RData")
load("DEV16S.tp3.div.brs.RData")
load("DEV16S.tp4.conv.brs.RData")
load("DEV16S.tp4.div.brs.RData")
###########################################################################################################################

#Add timepoint metadata
sig.1conv.brs$Timepoint <- rep("V4", nrow(sig.1conv.brs))
sig.1div.brs$Timepoint <- rep("V4", nrow(sig.1div.brs))
sig.2conv.brs$Timepoint <- rep("V11", nrow(sig.2conv.brs))
sig.2div.brs$Timepoint <- rep("V11", nrow(sig.2div.brs))
sig.3conv.brs$Timepoint <- rep("R2", nrow(sig.3conv.brs))
sig.3div.brs$Timepoint <- rep("R2", nrow(sig.3div.brs))
sig.4conv.brs$Timepoint <- rep("R5", nrow(sig.4conv.brs))
sig.4div.brs$Timepoint <- rep("R5", nrow(sig.4div.brs))

#Add soil metadata
sig.1conv.brs$Soil <- rep("Conventional", nrow(sig.1conv.brs))
sig.1div.brs$Soil <- rep("Diversified", nrow(sig.1div.brs))
sig.2conv.brs$Soil <- rep("Conventional", nrow(sig.2conv.brs))
sig.2div.brs$Soil <- rep("Diversified", nrow(sig.2div.brs))
sig.3conv.brs$Soil <- rep("Conventional", nrow(sig.3conv.brs))
sig.3div.brs$Soil <- rep("Diversified", nrow(sig.3div.brs))
sig.4conv.brs$Soil <- rep("Conventional", nrow(sig.4conv.brs))
sig.4div.brs$Soil <- rep("Diversified", nrow(sig.4div.brs))

#combine all bulk vs rs data
brs.all <- rbind(sig.1conv.brs, 
                 sig.1div.brs, 
                 sig.2conv.brs, 
                 sig.2div.brs, 
                 sig.3conv.brs, 
                 sig.3div.brs, 
                 sig.4conv.brs, 
                 sig.4div.brs)
head(brs.all)

#relevel stages to correct order
brs.all$Timepoint <- factor(brs.all$Timepoint, c("V4", "V11", "R2", "R5"))

#remove p__
brs.all$Phylum <- gsub("p__?", "", brs.all$Phylum)

scale_fill_discrete <- function(palname = "Set1", ...) {scale_fill_brewer(pallette = palname, ...)}
x.brs = tapply(brs.all$log2FoldChange, brs.all$Phylum, function(x) max(x))
x.brs = sort(x.brs, TRUE)
brs.all$Phylum = factor(as.character(brs.all$Phylum), levels=names(x.brs))

# Number of OTUs, see above data
label.brs <- data.frame(Phylum = 13, log2FoldChange = 6, 
                        label = c("Total OTUs = 17",
                                  "Total OTUs = 11",
                                  "Total OTUs = 80",
                                  "Total OTUs = 26",
                                  "Total OTUs = 56",
                                  "Total OTUs = 30",
                                  "Total OTUs = 50",
                                  "Total OTUs = 29"),
                        Timepoint = c("V4", "V4", "V11", "V11", "R2", "R2", "R5", "R5"),
                        Soil = c("Conventional", "Diversified", "Conventional", "Diversified", "Conventional", "Diversified", "Conventional", "Diversified"))
head(label.brs)

# graph
fig6 <- ggplot(brs.all, aes(x=Phylum, y=log2FoldChange, color= Phylum)) + 
  geom_jitter(size=0.75, width=0.1) + 
  facet_grid(Timepoint~Soil) +
  geom_hline(yintercept=0, color="grey69", linetype=2) +
  geom_text(data = label.brs, aes(Phylum, log2FoldChange, label=label), inherit.aes=FALSE, size=2.5) +
  labs(title = "",
       y="Log2 fold change (compared to bulk soil)",
       x="") +
  theme_bw() +
  theme(axis.text.x = element_text(size=8, angle = 45, hjust=1),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        legend.position="none",
        strip.text.x = element_text(size = 8, margin=margin(3,3,3,3)),
        strip.text.y = element_text(size = 8, margin=margin(3,3,3,3)),
        plot.margin = unit(c(0, 0, 0, 1), "cm"))
fig6

ggsave("figure6.tiff", plot=fig6, scale=1, width=5, height=5, units="in")


# Supp figure 3

#########################################################################################################################
load("DEV16S.tp1.convdiv.b.RData")
load("DEV16S.tp2.convdiv.b.RData")
load("DEV16S.tp3.convdiv.b.RData")
load("DEV16S.tp4.convdiv.b.RData")
#########################################################################################################################

#Add timepoint metadata
sig.1conv.b$Timepoint <- rep("V4", nrow(sig.1conv.b))
sig.2conv.b$Timepoint <- rep("V11", nrow(sig.2conv.b))
sig.3conv.b$Timepoint <- rep("R2", nrow(sig.3conv.b))
sig.4conv.b$Timepoint <- rep("R5", nrow(sig.4conv.b))

#combine
b.all <- rbind(sig.1conv.b, sig.2conv.b, sig.3conv.b, sig.4conv.b)
head(b.all)

#relevel timepoints to correct order
b.all$Timepoint <- factor(b.all$Timepoint, c("V4", "V11", "R2", "R5"))

#remove p__
b.all$Phylum <- gsub("p__?", "", b.all$Phylum)

#Point graph
scale_fill_discrete <- function(palname = "Set1", ...) {scale_fill_brewer(pallette = palname, ...)}
x.b = tapply(b.all$log2FoldChange, b.all$Phylum, function(x) max(x))
x.b = sort(x.b, TRUE)
b.all$Phylum = factor(as.character(b.all$Phylum), levels=names(x.b))

#Labels
label.b <- data.frame(Phylum = 14, log2FoldChange = 4.5, 
                      label = c("Total OTUs = 121",
                                "Total OTUs = 125",
                                "Total OTUs = 164",
                                "Total OTUs = 107"),
                      Timepoint = c("V4", "V11", "R2", "R5"))
head(label.b)

# graph
suppfig3 <- ggplot(b.all, aes(x=Phylum, y=log2FoldChange, color=Phylum)) + 
  geom_jitter(size=0.5, width=0.1) + 
  facet_grid(Timepoint~.) +
  geom_hline(yintercept=0, color="grey69", linetype=2) +
  geom_text(data = label.b, aes(Phylum, log2FoldChange, label=label), inherit.aes=FALSE, size=2.5) +
  labs(title="", 
       y="Log2 fold change (compared to conv. bulk soil)",
       x="") +
  theme_bw() +
  theme(axis.text.x = element_text(size=8, angle = 45, hjust=1),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        legend.position="none",
        strip.text.x = element_text(size = 8, margin=margin(3,3,3,3)),
        strip.text.y = element_text(size = 8, margin=margin(3,3,3,3)),
        plot.margin = unit(c(0, 0, 0, 1), "cm"))
suppfig3

ggsave("suppfig3.tiff", plot=suppfig3, scale=1, width=7, height=4, units="in")


# Supp fig 5

########################################################################################################################
load("DEV16S.tp1.convdiv.rs.RData")
load("DEV16S.tp2.convdiv.rs.RData")
load("DEV16S.tp3.convdiv.rs.RData")
load("DEV16S.tp4.convdiv.rs.RData")
#########################################################################################################################

#add timepoint metadata
tp1.diff4.rs$Timepoint <- rep("V4", nrow(tp1.diff4.rs))
tp2.diff4.rs$Timepoint <- rep("V11", nrow(tp2.diff4.rs))
tp3.diff4.rs$Timepoint <- rep("R2", nrow(tp3.diff4.rs))
tp4.diff4.rs$Timepoint <- rep("R5", nrow(tp4.diff4.rs))

#Combine all rs vs rs data
rs.all <- rbind(tp1.diff4.rs, tp2.diff4.rs, tp3.diff4.rs, tp4.diff4.rs)
head(rs.all)

#relevel stages to correct order
rs.all$Timepoint <- factor(rs.all$Timepoint, c("V4", "V11", "R2", "R5"))

#remove p__
rs.all$Phylum <- gsub("p__?", "", rs.all$Phylum)

scale_fill_discrete <- function(palname = "Set1", ...) {scale_fill_brewer(pallette = palname, ...)}
x.rs = tapply(rs.all$log2FoldChange, rs.all$Phylum, function(x) max(x))
x.rs = sort(x.rs, TRUE)
rs.all$Phylum = factor(as.character(rs.all$Phylum), levels=names(x.rs))

#Labels
label.rs <- data.frame(Phylum = 4.5, log2FoldChange = 5, 
                       label = c("Total OTUs = 9",
                                 "Total OTUs = 17",
                                 "Total OTUs = 4",
                                 "Total OTUs = 3"),
                       Timepoint = c("V4", "V11", "R2", "R5"))
head(label.rs)

# graph
suppfig5 <- ggplot(rs.all, aes(x=Phylum, y=log2FoldChange, color=Phylum)) + 
  geom_jitter(size=0.5, width=0.1) + 
  facet_grid(Timepoint~.) +
  geom_hline(yintercept=0, color="grey69", linetype=2) +
  geom_text(data = label.rs, aes(Phylum, log2FoldChange, label=label), inherit.aes=FALSE, size=2.5) +
  labs(title="", 
       y="Log2 fold change (compared to conv. rhizopshere soil)",
       x="") +
  theme_bw() +
  theme(axis.text.x = element_text(size=8, angle = 45, hjust=1),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        legend.position="none",
        strip.text.x = element_text(size = 8, margin=margin(3,3,3,3)),
        strip.text.y = element_text(size = 8, margin=margin(3,3,3,3)),
        plot.margin = unit(c(0, 0, 0, 1), "cm"))
suppfig5

ggsave("suppfig5.tiff", plot=suppfig5, scale=1, width=3.25, height=4, units="in")


# Supp fig 8

#######################################################################################################################
load("DEV16S.tp1.conv.brp.RData")
load("DEV16S.tp1.div.brp.RData")
load("DEV16S.tp2.conv.brp.RData")
load("DEV16S.tp2.div.brp.RData")
load("DEV16S.tp3.conv.brp.RData")
load("DEV16S.tp3.div.brp.RData")
load("DEV16S.tp4.conv.brp.RData")
load("DEV16S.tp4.div.brp.RData")
########################################################################################################################

#Add timepoint metadata
sig.1conv.brp$Timepoint <- rep("V4", nrow(sig.1conv.brp))
sig.1div.brp$Timepoint <- rep("V4", nrow(sig.1div.brp))
sig.2conv.brp$Timepoint <- rep("V11", nrow(sig.2conv.brp))
sig.2div.brp$Timepoint <- rep("V11", nrow(sig.2div.brp))
sig.3conv.brp$Timepoint <- rep("R2", nrow(sig.3conv.brp))
sig.3div.brp$Timepoint <- rep("R2", nrow(sig.3div.brp))
sig.4conv.brp$Timepoint <- rep("R5", nrow(sig.4conv.brp))
sig.4div.brp$Timepoint <- rep("R5", nrow(sig.4div.brp))

#Add soil metadata
sig.1conv.brp$Soil <- rep("Conventional", nrow(sig.1conv.brp))
sig.1div.brp$Soil <- rep("Diversified", nrow(sig.1div.brp))
sig.2conv.brp$Soil <- rep("Conventional", nrow(sig.2conv.brp))
sig.2div.brp$Soil <- rep("Diversified", nrow(sig.2div.brp))
sig.3conv.brp$Soil <- rep("Conventional", nrow(sig.3conv.brp))
sig.3div.brp$Soil <- rep("Diversified", nrow(sig.3div.brp))
sig.4conv.brp$Soil <- rep("Conventional", nrow(sig.4conv.brp))
sig.4div.brp$Soil <- rep("Diversified", nrow(sig.4div.brp))

#combine all bulk vs rp data
brp.all <- rbind(sig.1conv.brp, 
                 sig.1div.brp, 
                 sig.2conv.brp, 
                 sig.2div.brp, 
                 sig.3conv.brp, 
                 sig.3div.brp, 
                 sig.4conv.brp, 
                 sig.4div.brp)
head(brp.all)

#relevel stages to correct order
brp.all$Timepoint <- factor(brp.all$Timepoint, c("V4", "V11", "R2", "R5"))

#remove p__
brp.all$Phylum <- gsub("p__?", "", brp.all$Phylum)

scale_fill_discrete <- function(palname = "Set1", ...) {scale_fill_brewer(pallette = palname, ...)}
x.brp = tapply(brp.all$log2FoldChange, brp.all$Phylum, function(x) max(x))
x.brp = sort(x.brp, TRUE)
brp.all$Phylum = factor(as.character(brp.all$Phylum), levels=names(x.brp))

#For labeling
label.brp <- data.frame(Phylum = 21, log2FoldChange = 9.5, 
                        label = c("Total OTUs = 705",
                                  "Total OTUs = 556",
                                  "Total OTUs = 708",
                                  "Total OTUs = 578",
                                  "Total OTUs = 674",
                                  "Total OTUs = 685",
                                  "Total OTUs = 683",
                                  "Total OTUs = 664"),
                        Timepoint = c("V4", "V4", "V11", "V11", "R2", "R2", "R5", "R5"),
                        Soil = c("Conventional", "Diversified", "Conventional", "Diversified", "Conventional", "Diversified", "Conventional", "Diversified"))
head(label.brp)

# graph
suppfig8 <- ggplot(brp.all, aes(x=Phylum, y=log2FoldChange, color=Phylum)) + 
  geom_jitter(size=0.5, width=0.1) + 
  facet_grid(Timepoint~Soil) +
  geom_hline(yintercept=0, color="grey69", linetype=2) +
  geom_text(data = label.brp, aes(Phylum, log2FoldChange, label=label), inherit.aes=FALSE, size=2.5) +
  labs(title="", 
       y="Log2 fold change (compared to bulk soil)",
       x="") +
  theme_bw() +
  theme(axis.text.x = element_text(size=8, angle = 45, hjust=1),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        legend.position="none",
        strip.text.x = element_text(size = 8, margin=margin(3,3,3,3)),
        strip.text.y = element_text(size = 8, margin=margin(3,3,3,3)),
        plot.margin = unit(c(0, 0, 0, 1), "cm"))
suppfig8

ggsave("suppfig8.tiff", plot=suppfig8, scale=1, width=7, height=4, units="in")


# Supp fig 11

############################################################################################################################

load("DEV16S.tp1.convdiv.rp.RData")
load("DEV16S.tp2.convdiv.rp.RData")
load("DEV16S.tp3.convdiv.rp.RData")
load("DEV16S.tp4.convdiv.rp.RData")

############################################################################################################################

#add timepoint metadata
tp1.diff4.rp$Timepoint <- rep("V4", nrow(tp1.diff4.rp))
tp2.diff4.rp$Timepoint <- rep("V11", nrow(tp2.diff4.rp))
tp3.diff4.rp$Timepoint <- rep("R2", nrow(tp3.diff4.rp))
tp4.diff4.rp$Timepoint <- rep("R5", nrow(tp4.diff4.rp))

#Combine all rp vs rp data
rp.all <- rbind(tp1.diff4.rp, tp2.diff4.rp, tp3.diff4.rp, tp4.diff4.rp)
head(rp.all)

#relevel stages to correct order
rp.all$Timepoint <- factor(rp.all$Timepoint, c("V4", "V11", "R2", "R5"))

#remove p__
rp.all$Phylum <- gsub("p__?", "", rp.all$Phylum)

scale_fill_discrete <- function(palname = "Set1", ...) {scale_fill_brewer(pallette = palname, ...)}
x.rp = tapply(rp.all$log2FoldChange, rp.all$Phylum, function(x) max(x))
x.rp = sort(x.rp, TRUE)
rp.all$Phylum = factor(as.character(rp.all$Phylum), levels=names(x.rp))

#Labels
label.rp <- data.frame(Phylum = 10, log2FoldChange = 4.5, 
                       label = c("Total OTUs = 53",
                                 "Total OTUs = 78",
                                 "Total OTUs = 37",
                                 "Total OTUs = 19"),
                       Timepoint = c("V4", "V11", "R2", "R5"))
head(label.rp)

# graph
suppfig11 <- ggplot(rp.all, aes(x=Phylum, y=log2FoldChange, color=Phylum)) + 
  geom_jitter(size=0.5, width=0.1) + 
  facet_grid(Timepoint~.) +
  geom_hline(yintercept=0, color="grey69", linetype=2) +
  geom_text(data = label.rp, aes(Phylum, log2FoldChange, label=label), inherit.aes=FALSE, size=2.5) +
  labs(title="", 
       y="Log2 fold change (compared to conv. rhizoplane)",
       x="") +
  theme_bw() +
  theme(axis.text.x = element_text(size=8, angle = 45, hjust=1),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        legend.position="none",
        strip.text.x = element_text(size = 8, margin=margin(3,3,3,3)),
        strip.text.y = element_text(size = 8, margin=margin(3,3,3,3)),
        plot.margin = unit(c(0, 0, 0, 1), "cm"))
suppfig11

ggsave("suppfig11.tiff", plot=suppfig11, scale=1, width=3.25, height=4, units="in")


#########
#ITS Data

rm(list=ls())


############
#DESeq2 test

##Load in phyloseq object for DESeq2
#Created in "DEV16S_deseq2_prep_clean.R"
load("DEVITS.deseq.want.Rdata")
deseq.want

#convert phyloseq object to deseq2 object
deseq = phyloseq_to_deseq2(deseq.want, ~ group)
deseq

test = DESeq(deseq, test="Wald", fitType="parametric")

resultsNames(test)

#Set alpha level to 0.01
#This reduces family-wise error to 1%
alpha = 0.01


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
sig.1conv.brp = cbind(as(sig.1conv.brp, "data.frame"), as(tax_table(deseq.want)[rownames(sig.1conv.brp),], "matrix"))
head(sig.1conv.brp)
dim(sig.1conv.brp)
#69 diff abund taxa

#Div. bulk vs Div. rhizoplane
result.1div.brp = results(test, contrast=c("group", "TP1Div.rhizoplane", "TP1Div.bulk"))
mcols(result.1div.brp, use.names=TRUE)
#Extract only significant results
sig.1div.brp = result.1div.brp[which(result.1div.brp$padj < alpha),]
sig.1div.brp = cbind(as(sig.1div.brp, "data.frame"), as(tax_table(deseq.want)[rownames(sig.1div.brp),], "matrix"))
head(sig.1div.brp)
dim(sig.1div.brp)
#60 diff abund taxa

#Conv. rhizoplane vs Div. rhizoplane
#I chose Conv. rhizoplane to be the baseline for this comparison
result.1convdiv.rp = results(test, contrast=c("group", "TP1Div.rhizoplane", "TP1Conv.rhizoplane"))
mcols(result.1convdiv.rp, use.names=TRUE)
#extract only significant results
sig.1convdiv.rp = result.1convdiv.rp[which(result.1convdiv.rp$padj < alpha),]
sig.1convdiv.rp = cbind(as(sig.1convdiv.rp, "data.frame"), as(tax_table(deseq.want)[rownames(sig.1convdiv.rp),], "matrix"))
head(sig.1convdiv.rp)
dim(sig.1convdiv.rp)
#35 diff abund taxa 

#But, I'm interested in taxa that responded to the root, many of these differences probably also exist in the bulk soil, meaning the root wasn't the cause of the difference
#I want to remove those
#The code below retains only those OTUs that are significant in the rhizoplane vs rhizoplane comparison AND bulk vs rhizoplane comparison
tp1.diff1.rp <- sig.1convdiv.rp[(rownames(sig.1convdiv.rp) %in% rownames(sig.1conv.brp)),]
head(tp1.diff1.rp)
dim(tp1.diff1.rp)
#9 diff abund taxa between Conv. bulk and Conv. rhizoplane
tp1.diff2.rp <- sig.1convdiv.rp[(rownames(sig.1convdiv.rp) %in% rownames(sig.1div.brp)),]
head(tp1.diff2.rp)
dim(tp1.diff2.rp)
#10 diff abund taxa between Div. bulk and Div. rhizoplane
#Find OTUs that are not shared in common between diff1 and diff2
tp1.diff3.rp <- tp1.diff2.rp[!(rownames(tp1.diff2.rp) %in% rownames(tp1.diff1.rp)),]
head(tp1.diff3.rp)
dim(tp1.diff3.rp)
#5 unique OTUs
#Add unique OTUs to diff1
tp1.diff4.rp <- rbind(tp1.diff1.rp, tp1.diff3.rp)
head(tp1.diff4.rp)
dim(tp1.diff4.rp)
#14 diff abund taxa between rhizoplanes that are also diff abund from the bulk soils

############
#Timepoint 2

#Conv. bulk vs Conv. rhizoplane
result.2conv.brp = results(test, contrast=c("group", "TP2Conv.rhizoplane", "TP2Conv.bulk"))
mcols(result.2conv.brp, use.names=TRUE)
sig.2conv.brp = result.2conv.brp[which(result.2conv.brp$padj < alpha),]
sig.2conv.brp = cbind(as(sig.2conv.brp, "data.frame"), as(tax_table(deseq.want)[rownames(sig.2conv.brp),], "matrix"))
head(sig.2conv.brp)
dim(sig.2conv.brp)
#51 diff abund taxa

#Div. bulk vs Div. rhizoplane
result.2div.brp = results(test, contrast=c("group", "TP2Div.rhizoplane", "TP2Div.bulk"))
mcols(result.2div.brp, use.names=TRUE)
#Extract only significant results
sig.2div.brp = result.2div.brp[which(result.2div.brp$padj < alpha),]
sig.2div.brp = cbind(as(sig.2div.brp, "data.frame"), as(tax_table(deseq.want)[rownames(sig.2div.brp),], "matrix"))
head(sig.2div.brp)
dim(sig.2div.brp)
#58 diff abund taxa

#Conv. rhizoplane vs Div. rhizoplane
#I chose Conv. rhizoplane to be the baseline for this comparison
result.2convdiv.rp = results(test, contrast=c("group", "TP2Div.rhizoplane", "TP2Conv.rhizoplane"))
mcols(result.2convdiv.rp, use.names=TRUE)
#extract only significant results
sig.2convdiv.rp = result.2convdiv.rp[which(result.2convdiv.rp$padj < alpha),]
sig.2convdiv.rp = cbind(as(sig.2convdiv.rp, "data.frame"), as(tax_table(deseq.want)[rownames(sig.2convdiv.rp),], "matrix"))
head(sig.2convdiv.rp)
dim(sig.2convdiv.rp)
#38 diff abund taxa 

#But, I'm interested in taxa that responded to the root, many of these differences probably also exist in the bulk soil, meaning the root wasn't the cause of the difference
#I want to remove those
#The code below retains only those OTUs that are significant in the rhizoplane vs rhizoplane comparison AND bulk vs rhizoplane comparison
tp2.diff1.rp <- sig.2convdiv.rp[(rownames(sig.2convdiv.rp) %in% rownames(sig.2conv.brp)),]
head(tp2.diff1.rp)
dim(tp2.diff1.rp)
#69 diff abund taxa between Conv. bulk and Conv. rhizoplane
tp2.diff2.rp <- sig.2convdiv.rp[(rownames(sig.2convdiv.rp) %in% rownames(sig.2div.brp)),]
head(tp2.diff2.rp)
dim(tp2.diff2.rp)
#7 diff abund taxa between Div. bulk and Div. rhizoplane
#Find OTUs that are not shared in common between diff1 and diff2
tp2.diff3.rp <- tp2.diff2.rp[!(rownames(tp2.diff2.rp) %in% rownames(tp2.diff1.rp)),]
head(tp2.diff3.rp)
dim(tp2.diff3.rp)
#5 unique OTUs
#Add unique OTUs to diff1
tp2.diff4.rp <- rbind(tp2.diff1.rp, tp2.diff3.rp)
head(tp2.diff4.rp)
dim(tp2.diff4.rp)
#8 diff abund taxa between rhizoplanes that are also diff abund from the bulk soils

############
#Timepoint 3

#Conv. bulk vs Conv. rhizoplane
result.3conv.brp = results(test, contrast=c("group", "TP3Conv.rhizoplane", "TP3Conv.bulk"))
mcols(result.3conv.brp, use.names=TRUE)
#Extract only significant results
sig.3conv.brp = result.3conv.brp[which(result.3conv.brp$padj < alpha),]
sig.3conv.brp = cbind(as(sig.3conv.brp, "data.frame"), as(tax_table(deseq.want)[rownames(sig.3conv.brp),], "matrix"))
head(sig.3conv.brp)
dim(sig.3conv.brp)
#71 diff abund taxa

#Div. bulk vs Div. rhizoplane
result.3div.brp = results(test, contrast=c("group", "TP3Div.rhizoplane", "TP3Div.bulk"))
mcols(result.3div.brp, use.names=TRUE)
#Extract only significant results
sig.3div.brp = result.3div.brp[which(result.3div.brp$padj < alpha),]
sig.3div.brp = cbind(as(sig.3div.brp, "data.frame"), as(tax_table(deseq.want)[rownames(sig.3div.brp),], "matrix"))
head(sig.3div.brp)
dim(sig.3div.brp)
#56 diff abund taxa

#Conv. rhizoplane vs Div. rhizoplane
#I chose Conv. rhizoplane to be the baseline for this comparison
result.3convdiv.rp = results(test, contrast=c("group", "TP3Div.rhizoplane", "TP3Conv.rhizoplane"))
mcols(result.3convdiv.rp, use.names=TRUE)
#extract only significant results
sig.3convdiv.rp = result.3convdiv.rp[which(result.3convdiv.rp$padj < alpha),]
sig.3convdiv.rp = cbind(as(sig.3convdiv.rp, "data.frame"), as(tax_table(deseq.want)[rownames(sig.3convdiv.rp),], "matrix"))
head(sig.3convdiv.rp)
dim(sig.3convdiv.rp)
#25 diff abund taxa 

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
#3 diff abund taxa between Div. bulk and Div. rhizoplane
#Find OTUs that are not shared in common between diff1 and diff2
tp3.diff3.rp <- tp3.diff2.rp[!(rownames(tp3.diff2.rp) %in% rownames(tp3.diff1.rp)),]
head(tp3.diff3.rp)
dim(tp3.diff3.rp)
#1 unique OTUs
#Add unique OTUs to diff1
tp3.diff4.rp <- rbind(tp3.diff1.rp, tp3.diff3.rp)
head(tp3.diff4.rp)
dim(tp3.diff4.rp)
#7 diff abund taxa between rhizoplanes that are also diff abund from the bulk soils

############
#Timepoint 4

#Conv. bulk vs Conv. rhizoplane
result.4conv.brp = results(test, contrast=c("group", "TP4Conv.rhizoplane", "TP4Conv.bulk"))
mcols(result.4conv.brp, use.names=TRUE)
#Extract only significant results
sig.4conv.brp = result.4conv.brp[which(result.4conv.brp$padj < alpha),]
sig.4conv.brp = cbind(as(sig.4conv.brp, "data.frame"), as(tax_table(deseq.want)[rownames(sig.4conv.brp),], "matrix"))
head(sig.4conv.brp)
dim(sig.4conv.brp)
#49 diff abund taxa

#Div. bulk vs Div. rhizoplane
result.4div.brp = results(test, contrast=c("group", "TP4Div.rhizoplane", "TP4Div.bulk"))
mcols(result.4div.brp, use.names=TRUE)
#Extract only significant results
sig.4div.brp = result.4div.brp[which(result.4div.brp$padj < alpha),]
sig.4div.brp = cbind(as(sig.4div.brp, "data.frame"), as(tax_table(deseq.want)[rownames(sig.4div.brp),], "matrix"))
head(sig.4div.brp)
dim(sig.4div.brp)
#53 diff abund taxa

#Conv. rhizoplane vs Div. rhizoplane
#I chose Conv. rhizoplane to be the baseline for this comparison
result.4convdiv.rp = results(test, contrast=c("group", "TP4Div.rhizoplane", "TP4Conv.rhizoplane"))
mcols(result.4convdiv.rp, use.names=TRUE)
#extract only significant results
sig.4convdiv.rp = result.4convdiv.rp[which(result.4convdiv.rp$padj < alpha),]
sig.4convdiv.rp = cbind(as(sig.4convdiv.rp, "data.frame"), as(tax_table(deseq.want)[rownames(sig.4convdiv.rp),], "matrix"))
head(sig.4convdiv.rp)
dim(sig.4convdiv.rp)
#21 diff abund taxa 

#But, I'm interested in taxa that responded to the root, many of these differences probably also exist in the bulk soil, meaning the root wasn't the cause of the difference
#I want to remove those
#The code below retains only those OTUs that are significant in the rhizoplane vs rhizoplane comparison AND bulk vs rhizoplane comparison
tp4.diff1.rp <- sig.4convdiv.rp[(rownames(sig.4convdiv.rp) %in% rownames(sig.4conv.brp)),]
head(tp4.diff1.rp)
dim(tp4.diff1.rp)
#4 diff abund taxa between Conv. bulk and Conv. rhizoplane
tp4.diff2.rp <- sig.4convdiv.rp[(rownames(sig.4convdiv.rp) %in% rownames(sig.4div.brp)),]
head(tp4.diff2.rp)
dim(tp4.diff2.rp)
#5 diff abund taxa between Div. bulk and Div. rhizoplane
#Find OTUs that are not shared in common between diff1 and diff2
tp4.diff3.rp <- tp4.diff2.rp[!(rownames(tp4.diff2.rp) %in% rownames(tp4.diff1.rp)),]
head(tp4.diff3.rp)
dim(tp4.diff3.rp)
#3 unique OTUs
#Add unique OTUs to diff1
tp4.diff4.rp <- rbind(tp4.diff1.rp, tp4.diff3.rp)
head(tp4.diff4.rp)
dim(tp4.diff4.rp)
#7 diff abund taxa between rhizoplanes that are also diff abund from the bulk soils

########################
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
sig.1conv.brs = cbind(as(sig.1conv.brs, "data.frame"), as(tax_table(deseq.want)[rownames(sig.1conv.brs),], "matrix"))
head(sig.1conv.brs)
dim(sig.1conv.brs)
#15 diff abund taxa

#Div. bulk vs Div. rhizosphere
result.1div.brs = results(test, contrast=c("group", "TP1Div.rhizosphere", "TP1Div.bulk"))
mcols(result.1div.brs, use.names=TRUE)
#Extract only significant results
sig.1div.brs = result.1div.brs[which(result.1div.brs$padj < alpha),]
sig.1div.brs = cbind(as(sig.1div.brs, "data.frame"), as(tax_table(deseq.want)[rownames(sig.1div.brs),], "matrix"))
head(sig.1div.brs)
dim(sig.1div.brs)
#7 diff abund taxa

#Conv. rhizosphere vs Div. rhizosphere
#I chose Conv. rhizosphere to be the baseline for this comparison
result.1convdiv.rs = results(test, contrast=c("group", "TP1Div.rhizosphere", "TP1Conv.rhizosphere"))
mcols(result.1convdiv.rs, use.names=TRUE)
#extract only significant results
sig.1convdiv.rs = result.1convdiv.rs[which(result.1convdiv.rs$padj < alpha),]
sig.1convdiv.rs = cbind(as(sig.1convdiv.rs, "data.frame"), as(tax_table(deseq.want)[rownames(sig.1convdiv.rs),], "matrix"))
head(sig.1convdiv.rs)
dim(sig.1convdiv.rs)
#70 diff abund taxa 

#But, I'm interested in taxa that responded to the root, many of these differences probably also exist in the bulk soil, meaning the root wasn't the cause of the difference
#I want to remove those
#The code below retains only those OTUs that are significant in the rhizosphere vs rhizosphere comparison AND bulk vs rhizosphere comparison
tp1.diff1.rs <- sig.1convdiv.rs[(rownames(sig.1convdiv.rs) %in% rownames(sig.1conv.brs)),]
head(tp1.diff1.rs)
dim(tp1.diff1.rs)
#7 diff abund taxa between Conv. bulk and Conv. rhizosphere
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
#8 diff abund taxa between rhizospheres that are also diff abund from the bulk soils

############
#Timepoint 2

#Conv. bulk vs Conv. rhizosphere
result.2conv.brs = results(test, contrast=c("group", "TP2Conv.rhizosphere", "TP2Conv.bulk"))
mcols(result.2conv.brs, use.names=TRUE)
sig.2conv.brs = result.2conv.brs[which(result.2conv.brs$padj < alpha),]
sig.2conv.brs = cbind(as(sig.2conv.brs, "data.frame"), as(tax_table(deseq.want)[rownames(sig.2conv.brs),], "matrix"))
head(sig.2conv.brs)
dim(sig.2conv.brs)
#1 diff abund taxa

#Div. bulk vs Div. rhizosphere
result.2div.brs = results(test, contrast=c("group", "TP2Div.rhizosphere", "TP2Div.bulk"))
mcols(result.2div.brs, use.names=TRUE)
#Extract only significant results
sig.2div.brs = result.2div.brs[which(result.2div.brs$padj < alpha),]
sig.2div.brs = cbind(as(sig.2div.brs, "data.frame"), as(tax_table(deseq.want)[rownames(sig.2div.brs),], "matrix"))
head(sig.2div.brs)
dim(sig.2div.brs)
#5 diff abund taxa

#Conv. rhizosphere vs Div. rhizosphere
#I chose Conv. rhizosphere to be the baseline for this comparison
result.2convdiv.rs = results(test, contrast=c("group", "TP2Div.rhizosphere", "TP2Conv.rhizosphere"))
mcols(result.2convdiv.rs, use.names=TRUE)
#extract only significant results
sig.2convdiv.rs = result.2convdiv.rs[which(result.2convdiv.rs$padj < alpha),]
sig.2convdiv.rs = cbind(as(sig.2convdiv.rs, "data.frame"), as(tax_table(deseq.want)[rownames(sig.2convdiv.rs),], "matrix"))
head(sig.2convdiv.rs)
dim(sig.2convdiv.rs)
#73 diff abund taxa 

#But, I'm interested in taxa that responded to the root, many of these differences probably also exist in the bulk soil, meaning the root wasn't the cause of the difference
#I want to remove those
#The code below retains only those OTUs that are significant in the rhizosphere vs rhizosphere comparison AND bulk vs rhizosphere comparison
tp2.diff1.rs <- sig.2convdiv.rs[(rownames(sig.2convdiv.rs) %in% rownames(sig.2conv.brs)),]
head(tp2.diff1.rs)
dim(tp2.diff1.rs)
#1 diff abund taxa between Conv. bulk and Conv. rhizosphere
tp2.diff2.rs <- sig.2convdiv.rs[(rownames(sig.2convdiv.rs) %in% rownames(sig.2div.brs)),]
head(tp2.diff2.rs)
dim(tp2.diff2.rs)
#1 diff abund taxa between Div. bulk and Div. rhizosphere
#Find OTUs that are not shared in common between diff1 and diff2
tp2.diff3.rs <- tp2.diff2.rs[!(rownames(tp2.diff2.rs) %in% rownames(tp2.diff1.rs)),]
head(tp2.diff3.rs)
dim(tp2.diff3.rs)
#1 unique OTUs
#Add unique OTUs to diff1
tp2.diff4.rs <- rbind(tp2.diff1.rs, tp2.diff3.rs)
head(tp2.diff4.rs)
dim(tp2.diff4.rs)
#2 diff abund taxa between rhizospheres that are also diff abund from the bulk soils

############
#Timepoint 3

#Conv. bulk vs Conv. rhizosphere
result.3conv.brs = results(test, contrast=c("group", "TP3Conv.rhizosphere", "TP3Conv.bulk"))
mcols(result.3conv.brs, use.names=TRUE)
#Extract only significant results
sig.3conv.brs = result.3conv.brs[which(result.3conv.brs$padj < alpha),]
sig.3conv.brs = cbind(as(sig.3conv.brs, "data.frame"), as(tax_table(deseq.want)[rownames(sig.3conv.brs),], "matrix"))
head(sig.3conv.brs)
dim(sig.3conv.brs)
#7 diff abund taxa

#Div. bulk vs Div. rhizosphere
result.3div.brs = results(test, contrast=c("group", "TP3Div.rhizosphere", "TP3Div.bulk"))
mcols(result.3div.brs, use.names=TRUE)
#Extract only significant results
sig.3div.brs = result.3div.brs[which(result.3div.brs$padj < alpha),]
sig.3div.brs = cbind(as(sig.3div.brs, "data.frame"), as(tax_table(deseq.want)[rownames(sig.3div.brs),], "matrix"))
head(sig.3div.brs)
dim(sig.3div.brs)
#11 diff abund taxa

#Conv. rhizosphere vs Div. rhizosphere
#I chose Conv. rhizosphere to be the baseline for this comparison
result.3convdiv.rs = results(test, contrast=c("group", "TP3Div.rhizosphere", "TP3Conv.rhizosphere"))
mcols(result.3convdiv.rs, use.names=TRUE)
#extract only significant results
sig.3convdiv.rs = result.3convdiv.rs[which(result.3convdiv.rs$padj < alpha),]
sig.3convdiv.rs = cbind(as(sig.3convdiv.rs, "data.frame"), as(tax_table(deseq.want)[rownames(sig.3convdiv.rs),], "matrix"))
head(sig.3convdiv.rs)
dim(sig.3convdiv.rs)
#41 diff abund taxa 

#But, I'm interested in taxa that responded to the root, many of these differences probably also exist in the bulk soil, meaning the root wasn't the cause of the difference
#I want to remove those
#The code below retains only those OTUs that are significant in the rhizosphere vs rhizosphere comparison AND bulk vs rhizosphere comparison
tp3.diff1.rs <- sig.3convdiv.rs[(rownames(sig.3convdiv.rs) %in% rownames(sig.3conv.brs)),]
head(tp3.diff1.rs)
dim(tp3.diff1.rs)
#3 diff abund taxa between Conv. bulk and Conv. rhizosphere
tp3.diff2.rs <- sig.3convdiv.rs[(rownames(sig.3convdiv.rs) %in% rownames(sig.3div.brs)),]
head(tp3.diff2.rs)
dim(tp3.diff2.rs)
#1 diff abund taxa between Div. bulk and Div. rhizosphere
#Find OTUs that are not shared in common between diff1 and diff2
tp3.diff3.rs <- tp3.diff2.rs[!(rownames(tp3.diff2.rs) %in% rownames(tp3.diff1.rs)),]
head(tp3.diff3.rs)
dim(tp3.diff3.rs)
#1 unique OTUs
#Add unique OTUs to diff1
tp3.diff4.rs <- rbind(tp3.diff1.rs, tp3.diff3.rs)
head(tp3.diff4.rs)
dim(tp3.diff4.rs)
#4 diff abund taxa between rhizospheres that are also diff abund from the bulk soils

############
#Timepoint 4

#Conv. bulk vs Conv. rhizosphere
result.4conv.brs = results(test, contrast=c("group", "TP4Conv.rhizosphere", "TP4Conv.bulk"))
mcols(result.4conv.brs, use.names=TRUE)
#Extract only significant results
sig.4conv.brs = result.4conv.brs[which(result.4conv.brs$padj < alpha),]
sig.4conv.brs = cbind(as(sig.4conv.brs, "data.frame"), as(tax_table(deseq.want)[rownames(sig.4conv.brs),], "matrix"))
head(sig.4conv.brs)
dim(sig.4conv.brs)
#7 diff abund taxa

#Div. bulk vs Div. rhizosphere
result.4div.brs = results(test, contrast=c("group", "TP4Div.rhizosphere", "TP4Div.bulk"))
mcols(result.4div.brs, use.names=TRUE)
#Extract only significant results
sig.4div.brs = result.4div.brs[which(result.4div.brs$padj < alpha),]
sig.4div.brs = cbind(as(sig.4div.brs, "data.frame"), as(tax_table(deseq.want)[rownames(sig.4div.brs),], "matrix"))
head(sig.4div.brs)
dim(sig.4div.brs)
#15 diff abund taxa

#Conv. rhizosphere vs Div. rhizosphere
#I chose Conv. rhizosphere to be the baseline for this comparison
result.4convdiv.rs = results(test, contrast=c("group", "TP4Div.rhizosphere", "TP4Conv.rhizosphere"))
mcols(result.4convdiv.rs, use.names=TRUE)
#extract only significant results
sig.4convdiv.rs = result.4convdiv.rs[which(result.4convdiv.rs$padj < alpha),]
sig.4convdiv.rs = cbind(as(sig.4convdiv.rs, "data.frame"), as(tax_table(deseq.want)[rownames(sig.4convdiv.rs),], "matrix"))
head(sig.4convdiv.rs)
dim(sig.4convdiv.rs)
#34 diff abund taxa 

#But, I'm interested in taxa that responded to the root, many of these differences probably also exist in the bulk soil, meaning the root wasn't the cause of the difference
#I want to remove those
#The code below retains only those OTUs that are significant in the rhizosphere vs rhizosphere comparison AND bulk vs rhizosphere comparison
tp4.diff1.rs <- sig.4convdiv.rs[(rownames(sig.4convdiv.rs) %in% rownames(sig.4conv.brs)),]
head(tp4.diff1.rs)
dim(tp4.diff1.rs)
#3 diff abund taxa between Conv. bulk and Conv. rhizosphere
tp4.diff2.rs <- sig.4convdiv.rs[(rownames(sig.4convdiv.rs) %in% rownames(sig.4div.brs)),]
head(tp4.diff2.rs)
dim(tp4.diff2.rs)
#2 diff abund taxa between Div. bulk and Div. rhizosphere
#Find OTUs that are not shared in common between diff1 and diff2
tp4.diff3.rs <- tp4.diff2.rs[!(rownames(tp4.diff2.rs) %in% rownames(tp4.diff1.rs)),]
head(tp4.diff3.rs)
dim(tp4.diff3.rs)
#1 unique OTUs
#Add unique OTUs to diff1
tp4.diff4.rs <- rbind(tp4.diff1.rs, tp4.diff3.rs)
head(tp4.diff4.rs)
dim(tp4.diff4.rs)
#4 diff abund taxa between rhizospheres that are also diff abund from the bulk soils

#################
#Bulk comparisons

############
#Timepoint 1

#Conv. Bulk vs Div. Bulk
result.1conv.b = results(test, contrast=c("group", "TP1Div.bulk", "TP1Conv.bulk"))
mcols(result.1conv.b, use.names=TRUE)
#extract only significant results
sig.1conv.b = result.1conv.b[which(result.1conv.b$padj < alpha),]
sig.1conv.b = cbind(as(sig.1conv.b, "data.frame"), as(tax_table(deseq.want)[rownames(sig.1conv.b),], "matrix"))
head(sig.1conv.b)
dim(sig.1conv.b)
#66 diff abund taxa beteween bulk soils

############
#Timepoint 2

#Conv. Bulk vs Div. Bulk
result.2conv.b = results(test, contrast=c("group", "TP2Div.bulk", "TP2Conv.bulk"))
mcols(result.2conv.b, use.names=TRUE)
#extract only significant results
sig.2conv.b = result.2conv.b[which(result.2conv.b$padj < alpha),]
sig.2conv.b = cbind(as(sig.2conv.b, "data.frame"), as(tax_table(deseq.want)[rownames(sig.2conv.b),], "matrix"))
head(sig.2conv.b)
dim(sig.2conv.b)
#47 diff abund taxa beteween bulk soils

############
#Timepoint 3

#Conv. Bulk vs Div. Bulk
result.3conv.b = results(test, contrast=c("group", "TP3Div.bulk", "TP3Conv.bulk"))
mcols(result.3conv.b, use.names=TRUE)
#extract only significant results
sig.3conv.b = result.3conv.b[which(result.3conv.b$padj < alpha),]
sig.3conv.b = cbind(as(sig.3conv.b, "data.frame"), as(tax_table(deseq.want)[rownames(sig.3conv.b),], "matrix"))
head(sig.3conv.b)
dim(sig.3conv.b)
#62 diff abund taxa beteween bulk soils

############
#Timepoint 4

#Conv. Bulk vs Div. Bulk
result.4conv.b = results(test, contrast=c("group", "TP4Div.bulk", "TP4Conv.bulk"))
mcols(result.4conv.b, use.names=TRUE)
#extract only significant results
sig.4conv.b = result.4conv.b[which(result.4conv.b$padj < alpha),]
sig.4conv.b = cbind(as(sig.4conv.b, "data.frame"), as(tax_table(deseq.want)[rownames(sig.4conv.b),], "matrix"))
head(sig.4conv.b)
dim(sig.4conv.b)
#49 diff abund taxa beteween bulk soils


#####################
# Create figures

# Supp fig 4

#######################################################################################################################
load("DEVITS.tp1.convdiv.b.RData")
load("DEVITS.tp2.convdiv.b.RData")
load("DEVITS.tp3.convdiv.b.RData")
load("DEVITS.tp4.convdiv.b.RData")
#######################################################################################################################

#Add timepoint metadata
sig.1conv.b$Timepoint <- rep("V4", nrow(sig.1conv.b))
sig.2conv.b$Timepoint <- rep("V11", nrow(sig.2conv.b))
sig.3conv.b$Timepoint <- rep("R2", nrow(sig.3conv.b))
sig.4conv.b$Timepoint <- rep("R5", nrow(sig.4conv.b))

#combine
b.all <- rbind(sig.1conv.b, sig.2conv.b, sig.3conv.b, sig.4conv.b)
head(b.all)

#relevel timepoints to correct order
b.all$Timepoint <- factor(b.all$Timepoint, c("V4", "V11", "R2", "R5"))

#remove p__
b.all$Phylum <- gsub("p__?", "", b.all$Phylum)

#Point graph
scale_fill_discrete <- function(palname = "Set1", ...) {scale_fill_brewer(pallette = palname, ...)}
x.b = tapply(b.all$log2FoldChange, b.all$Phylum, function(x) max(x))
x.b = sort(x.b, TRUE)
b.all$Phylum = factor(as.character(b.all$Phylum), levels=names(x.b))

#Labels
label.b <- data.frame(Phylum = 4.5, log2FoldChange = 1, 
                      label = c("Total OTUs = 66",
                                "Total OTUs = 47",
                                "Total OTUs = 62",
                                "Total OTUs = 49"),
                      Timepoint = c("V4", "V11", "R2", "R5"))
head(label.b)

# graph
suppfig4 <- ggplot(b.all, aes(x=Phylum, y=log2FoldChange, color=Phylum)) + 
  geom_jitter(size=0.5, width=0.1) + 
  facet_grid(Timepoint~.) +
  geom_hline(yintercept=0, color="grey69", linetype=2) +
  geom_text(data = label.b, aes(Phylum, log2FoldChange, label=label), inherit.aes=FALSE, size=2.5) +
  labs(title="", 
       y="Log2 fold change (compared to conv. bulk soil)",
       x="") +
  theme_bw() +
  theme(axis.text.x = element_text(size=8, angle = 45, hjust=1),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        legend.position="none",
        strip.text.x = element_text(size = 8, margin=margin(3,3,3,3)),
        strip.text.y = element_text(size = 8, margin=margin(3,3,3,3)),
        plot.margin = unit(c(0, 0, 0, 1), "cm"))
suppfig4

ggsave("suppfig4.tiff", plot=suppfig4, scale=1, width=3.25, height=4, units="in")


# Supp fig 6
##########################################################################################################################33
load("DEVITS.tp1.convdiv.rs.RData")
load("DEVITS.tp2.convdiv.rs.RData")
load("DEVITS.tp3.convdiv.rs.RData")
load("DEVITS.tp4.convdiv.rs.RData")
###########################################################################################################################

#add timepoint metadata
tp1.diff4.rs$Timepoint <- rep("V4", nrow(tp1.diff4.rs))
tp2.diff4.rs$Timepoint <- rep("V11", nrow(tp2.diff4.rs))
tp3.diff4.rs$Timepoint <- rep("R2", nrow(tp3.diff4.rs))
tp4.diff4.rs$Timepoint <- rep("R5", nrow(tp4.diff4.rs))

#Combine all rs vs rs data
rs.all <- rbind(tp1.diff4.rs, tp2.diff4.rs, tp3.diff4.rs, tp4.diff4.rs)
head(rs.all)

#relevel stages to correct order
rs.all$Timepoint <- factor(rs.all$Timepoint, c("V4", "V11", "R2", "R5"))

#remove p__
rs.all$Phylum <- gsub("p__?", "", rs.all$Phylum)

scale_fill_discrete <- function(palname = "Set1", ...) {scale_fill_brewer(pallette = palname, ...)}
x.rs = tapply(rs.all$log2FoldChange, rs.all$Phylum, function(x) max(x))
x.rs = sort(x.rs, TRUE)
rs.all$Phylum = factor(as.character(rs.all$Phylum), levels=names(x.rs))

#Labels
label.rs <- data.frame(Phylum = 3, log2FoldChange = 4.5, 
                       label = c("Total OTUs = 8",
                                 "Total OTUs = 2",
                                 "Total OTUs = 4",
                                 "Total OTUs = 4"),
                       Timepoint = c("V4", "V11", "R2", "R5"))
head(label.rs)

#graph
suppfig6 <- ggplot(rs.all, aes(x=Phylum, y=log2FoldChange, color=Phylum)) + 
  geom_jitter(size=0.5, width=0.1) + 
  facet_grid(Timepoint~.) +
  geom_hline(yintercept=0, color="grey69", linetype=2) +
  geom_text(data = label.rs, aes(Phylum, log2FoldChange, label=label), inherit.aes=FALSE, size=2.5) +
  labs(title="", 
       y="Log2 fold change (compared to conv. rhizosphere soil)",
       x="") +
  theme_bw() +
  theme(axis.text.x = element_text(size=8, angle = 45, hjust=1),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        legend.position="none",
        strip.text.x = element_text(size = 8, margin=margin(3,3,3,3)),
        strip.text.y = element_text(size = 8, margin=margin(3,3,3,3)),
        plot.margin = unit(c(0, 0, 0, 1), "cm"))
suppfig6

ggsave("suppfig6.tiff", plot=suppfig6, scale=1, width=3.25, height=4, units="in")


# Supp fig 7

########################################################################################################################
load("DEVITS.tp1.conv.brs.RData")
load("DEVITS.tp1.div.brs.RData")
load("DEVITS.tp2.conv.brs.RData")
load("DEVITS.tp2.div.brs.RData")
load("DEVITS.tp3.conv.brs.RData")
load("DEVITS.tp3.div.brs.RData")
load("DEVITS.tp4.conv.brs.RData")
load("DEVITS.tp4.div.brs.RData")
########################################################################################################################

#Add timepoint metadata
sig.1conv.brs$Timepoint <- rep("V4", nrow(sig.1conv.brs))
sig.1div.brs$Timepoint <- rep("V4", nrow(sig.1div.brs))
sig.2conv.brs$Timepoint <- rep("V11", nrow(sig.2conv.brs))
sig.2div.brs$Timepoint <- rep("V11", nrow(sig.2div.brs))
sig.3conv.brs$Timepoint <- rep("R2", nrow(sig.3conv.brs))
sig.3div.brs$Timepoint <- rep("R2", nrow(sig.3div.brs))
sig.4conv.brs$Timepoint <- rep("R5", nrow(sig.4conv.brs))
sig.4div.brs$Timepoint <- rep("R5", nrow(sig.4div.brs))

#Add soil metadata
sig.1conv.brs$Soil <- rep("Conventional", nrow(sig.1conv.brs))
sig.1div.brs$Soil <- rep("Diversified", nrow(sig.1div.brs))
sig.2conv.brs$Soil <- rep("Conventional", nrow(sig.2conv.brs))
sig.2div.brs$Soil <- rep("Diversified", nrow(sig.2div.brs))
sig.3conv.brs$Soil <- rep("Conventional", nrow(sig.3conv.brs))
sig.3div.brs$Soil <- rep("Diversified", nrow(sig.3div.brs))
sig.4conv.brs$Soil <- rep("Conventional", nrow(sig.4conv.brs))
sig.4div.brs$Soil <- rep("Diversified", nrow(sig.4div.brs))

#combine all bulk vs rs data
brs.all <- rbind(sig.1conv.brs, 
                 sig.1div.brs, 
                 sig.2conv.brs, 
                 sig.2div.brs, 
                 sig.3conv.brs, 
                 sig.3div.brs, 
                 sig.4conv.brs, 
                 sig.4div.brs)
head(brs.all)

#relevel stages to correct order
brs.all$Timepoint <- factor(brs.all$Timepoint, c("V4", "V11", "R2", "R5"))

#remove p__
brs.all$Phylum <- gsub("p__?", "", brs.all$Phylum)

scale_fill_discrete <- function(palname = "Set1", ...) {scale_fill_brewer(pallette = palname, ...)}
x.brs = tapply(brs.all$log2FoldChange, brs.all$Phylum, function(x) max(x))
x.brs = sort(x.brs, TRUE)
brs.all$Phylum = factor(as.character(brs.all$Phylum), levels=names(x.brs))

#For labeling
label.brs <- data.frame(Phylum = 3.5, log2FoldChange = 6.75, 
                        label = c("Total OTUs = 15",
                                  "Total OTUs = 7",
                                  "Total OTUs = 1",
                                  "Total OTUs = 5",
                                  "Total OTUs = 7",
                                  "Total OTUs = 11",
                                  "Total OTUs = 7",
                                  "Total OTUs = 15"),
                        Timepoint = c("V4", "V4", "V11", "V11", "R2", "R2", "R5", "R5"),
                        Soil = c("Conventional", "Diversified", "Conventional", "Diversified", "Conventional", "Diversified", "Conventional", "Diversified"))
head(label.brs)

# graph
suppfig7 <- ggplot(brs.all, aes(x=Phylum, y=log2FoldChange, color=Phylum)) + 
  geom_jitter(size=0.5, width=0.1) + 
  facet_grid(Timepoint~Soil) +
  geom_hline(yintercept=0, color="grey69", linetype=2) +
  geom_text(data = label.brs, aes(Phylum, log2FoldChange, label=label), inherit.aes=FALSE, size=2.5) +
  labs(title="", 
       y="Log2 fold change (compared to bulk soil)",
       x="") +
  theme_bw() +
  theme(axis.text.x = element_text(size=8, angle = 45, hjust=1),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        legend.position="none",
        strip.text.x = element_text(size = 8, margin=margin(3,3,3,3)),
        strip.text.y = element_text(size = 8, margin=margin(3,3,3,3)),
        plot.margin = unit(c(0, 0, 0, 1), "cm"))
suppfig7

ggsave("suppfig7.tiff", plot=suppfig7, scale=1, width=5, height=4, units="in")


# Supp fig 9

######################################################################################################################
load("DEVITS.tp1.conv.brp.RData")
load("DEVITS.tp1.div.brp.RData")
load("DEVITS.tp2.conv.brp.RData")
load("DEVITS.tp2.div.brp.RData")
load("DEVITS.tp3.conv.brp.RData")
load("DEVITS.tp3.div.brp.RData")
load("DEVITS.tp4.conv.brp.RData")
load("DEVITS.tp4.div.brp.RData")
#####################################################################################################################

#Add timepoint metadata
sig.1conv.brp$Timepoint <- rep("V4", nrow(sig.1conv.brp))
sig.1div.brp$Timepoint <- rep("V4", nrow(sig.1div.brp))
sig.2conv.brp$Timepoint <- rep("V11", nrow(sig.2conv.brp))
sig.2div.brp$Timepoint <- rep("V11", nrow(sig.2div.brp))
sig.3conv.brp$Timepoint <- rep("R2", nrow(sig.3conv.brp))
sig.3div.brp$Timepoint <- rep("R2", nrow(sig.3div.brp))
sig.4conv.brp$Timepoint <- rep("R5", nrow(sig.4conv.brp))
sig.4div.brp$Timepoint <- rep("R5", nrow(sig.4div.brp))

#Add soil metadata
sig.1conv.brp$Soil <- rep("Conventional", nrow(sig.1conv.brp))
sig.1div.brp$Soil <- rep("Diversified", nrow(sig.1div.brp))
sig.2conv.brp$Soil <- rep("Conventional", nrow(sig.2conv.brp))
sig.2div.brp$Soil <- rep("Diversified", nrow(sig.2div.brp))
sig.3conv.brp$Soil <- rep("Conventional", nrow(sig.3conv.brp))
sig.3div.brp$Soil <- rep("Diversified", nrow(sig.3div.brp))
sig.4conv.brp$Soil <- rep("Conventional", nrow(sig.4conv.brp))
sig.4div.brp$Soil <- rep("Diversified", nrow(sig.4div.brp))

#combine all bulk vs rp data
brp.all <- rbind(sig.1conv.brp, 
                 sig.1div.brp, 
                 sig.2conv.brp, 
                 sig.2div.brp, 
                 sig.3conv.brp, 
                 sig.3div.brp, 
                 sig.4conv.brp, 
                 sig.4div.brp)
head(brp.all)


#relevel stages to correct order
brp.all$Timepoint <- factor(brp.all$Timepoint, c("V4", "V11", "R2", "R5"))

#remove p__
brp.all$Phylum <- gsub("p__?", "", brp.all$Phylum)

scale_fill_discrete <- function(palname = "Set1", ...) {scale_fill_brewer(pallette = palname, ...)}
x.brp = tapply(brp.all$log2FoldChange, brp.all$Phylum, function(x) max(x))
x.brp = sort(x.brp, TRUE)
brp.all$Phylum = factor(as.character(brp.all$Phylum), levels=names(x.brp))
#For labeling
label.brp <- data.frame(Phylum = 3.5, log2FoldChange = 9, 
                        label = c("Total OTUs = 69",
                                  "Total OTUs = 60",
                                  "Total OTUs = 51",
                                  "Total OTUs = 58",
                                  "Total OTUs = 71",
                                  "Total OTUs = 56",
                                  "Total OTUs = 49",
                                  "Total OTUs = 53"),
                        Timepoint = c("V4", "V4", "V11", "V11", "R2", "R2", "R5", "R5"),
                        Soil = c("Conventional", "Diversified", "Conventional", "Diversified", "Conventional", "Diversified", "Conventional", "Diversified"))
head(label.brp)

# graph
suppfig9 <- ggplot(brp.all, aes(x=Phylum, y=log2FoldChange, color=Phylum)) + 
  geom_jitter(size=0.5, width=0.1) + 
  facet_grid(Timepoint~Soil) +
  geom_hline(yintercept=0, color="grey69", linetype=2) +
  geom_text(data = label.brp, aes(Phylum, log2FoldChange, label=label), inherit.aes=FALSE, size=2.5) +
  labs(title="", 
       y="Log2 fold change (compared to bulk soil)",
       x="") +
  theme_bw() +
  theme(axis.text.x = element_text(size=8, angle = 45, hjust=1),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        legend.position="none",
        strip.text.x = element_text(size = 8, margin=margin(3,3,3,3)),
        strip.text.y = element_text(size = 8, margin=margin(3,3,3,3)),
        plot.margin = unit(c(0, 0, 0, 1), "cm"))
suppfig9

ggsave("suppfig9.tiff", plot=suppfig9, scale=1, width=5, height=4, units="in")


# Supp fig 10

###########################################################################################################################
load("DEVITS.tp1.convdiv.rp.RData")
load("DEVITS.tp2.convdiv.rp.RData")
load("DEVITS.tp3.convdiv.rp.RData")
load("DEVITS.tp4.convdiv.rp.RData")
##########################################################################################################################

#add timepoint metadata
tp1.diff4.rp$Timepoint <- rep("V4", nrow(tp1.diff4.rp))
tp2.diff4.rp$Timepoint <- rep("V11", nrow(tp2.diff4.rp))
tp3.diff4.rp$Timepoint <- rep("R2", nrow(tp3.diff4.rp))
tp4.diff4.rp$Timepoint <- rep("R5", nrow(tp4.diff4.rp))

#Combine all rp vs rp data
rp.all <- rbind(tp1.diff4.rp, tp2.diff4.rp, tp3.diff4.rp, tp4.diff4.rp)
head(rp.all)

#relevel stages to correct order
rp.all$Timepoint <- factor(rp.all$Timepoint, c("V4", "V11", "R2", "R5"))

#remove p__
rp.all$Phylum <- gsub("p__?", "", rp.all$Phylum)

scale_fill_discrete <- function(palname = "Set1", ...) {scale_fill_brewer(pallette = palname, ...)}
x.rp = tapply(rp.all$log2FoldChange, rp.all$Phylum, function(x) max(x))
x.rp = sort(x.rp, TRUE)
rp.all$Phylum = factor(as.character(rp.all$Phylum), levels=names(x.rp))
#Labels
label.rp <- data.frame(Phylum = 2, log2FoldChange = 5, 
                       label = c("Total OTUs = 14",
                                 "Total OTUs = 8",
                                 "Total OTUs = 7",
                                 "Total OTUs = 7"),
                       Timepoint = c("V4", "V11", "R2", "R5"))
head(label.rp)

#graph
suppfig10 <- ggplot(rp.all, aes(x=Phylum, y=log2FoldChange, color=Phylum)) + 
  geom_jitter(size=0.5, width=0.1) + 
  facet_grid(Timepoint~.) +
  geom_hline(yintercept=0, color="grey69", linetype=2) +
  geom_text(data = label.rp, aes(Phylum, log2FoldChange, label=label), inherit.aes=FALSE, size=2.5) +
  labs(title="", 
       y="Log2 fold change (compared to conv. rhizoplane)",
       x="") +
  theme_bw() +
  theme(axis.text.x = element_text(size=8, angle = 45, hjust=1),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        legend.position="none",
        strip.text.x = element_text(size = 8, margin=margin(3,3,3,3)),
        strip.text.y = element_text(size = 8, margin=margin(3,3,3,3)),
        plot.margin = unit(c(0, 0, 0, 1), "cm"))
suppfig10

ggsave("suppfig10.tiff", plot=suppfig10, scale=1, width=3.25, height=4, units="in")



# Figure 5

# metadata
comparison <- c(rep("Conv. vs Div. bulk", 4), rep("Conv. vs Div. rhizosphere", 4), rep("Conv. vs Div. rhizoplane", 4))
days <- as.numeric(rep(c(35,55,82,119), 3))

# numbers of differentially abundant prokaryotic OTUs (see above script)
value <- as.numeric(c(121,125,164,107,9,17,4,3,53,78,37,19))

# create dataframe
fig5.data16s <- cbind(comparison, days, value)
fig5.data16s <- as.data.frame(fig5.data16s)
head(fig5.data16s)

# add 16s label
fig5.data16s$micro <- rep("Prokaryotic", nrow(fig5.data16s))
head(fig5.data16s)

# number of diff abund fungal OTUs (see above script)
value <- as.numeric(rep(c(66,47,62,49,8,2,4,4,14,8,7,7)))

# create dataframe
fig5.dataits <- cbind(comparison, days, value)
fig5.dataits <- as.data.frame(fig5.dataits)
head(fig5.dataits)

# add 16s label
fig5.dataits$micro <- rep("Fungal", nrow(fig5.dataits))
head(fig5.dataits)

# combine data
fig5.data <- rbind(fig5.data16s, fig5.dataits)
head(fig5.data)
tail(fig5.data)

# fix numeric values
fig5.data$value <- as.matrix(fig5.data$value)
fig5.data$value <- as.numeric(fig5.data$value)
fig5.data$days <- as.matrix(fig5.data$days)
fig5.data$days <- as.numeric(fig5.data$days)

# relevel factors
fig5.data$micro <- factor(fig5.data$micro, c("Prokaryotic", "Fungal"))

# graph
fig5 <- ggplot(fig5.data, aes(x=days, y=value, shape=comparison)) +
  geom_point(size=2) +
  geom_line(aes(group=comparison)) +
  facet_wrap(~micro) +
  labs(title="", y="Number of differentially abundant OTUs", x="Days after planting") +
  theme_bw() +
  theme(axis.title=element_text(size=8),
        axis.text=element_text(size=8),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        panel.background = element_rect(colour="black", size=0.5),
        strip.text.x = element_text(size = 8, margin=margin(3,3,3,3)),
        legend.position="none",
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  guides(fill=FALSE)
fig5

ggsave("figure5.tiff", plot=fig5, scale=1, width=5, height=3, units="in")
