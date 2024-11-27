rm(list=ls())
setwd("")
library(TCGAmutations)
tmp=as.data.frame(tcga_available())
tmp
PRAD<-TCGAmutations::tcga_load(study = "PRAD")
PRAD

library(maftools)
library(dplyr)
library(ggplot2)
p30<-oncoplot(maf = PRAD,
         top = 30,          fontSize = 0.6,          showTumorSampleBarcodes = F) 
p30

gene<-read.csv("res_cox_sig.csv",header = T,row.names = 1)
head(gene)
gene_name<-colnames(gene)
oncoplot(maf = PRAD,
         genes=gene_name )








