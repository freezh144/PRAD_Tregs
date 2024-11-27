rm(list = ls())
setwd("")
library(psych)
library(pheatmap)
library(reshape2)
risk<-read.csv("risk_score.csv",header=TRUE,row.names = 1)
head(risk)
dim(risk)

M<-read.csv("survival_exp1.csv",header=TRUE,row.names = 1)
head(M)[1:6]
#M1<-t(M)
M1<-as.data.frame(M)
head(M1)[1:6]
M1<-M1[order(rownames(M1)),]
dim(M1)
head(M1)[1:6]

cor = corr.test(M1, risk, method="spearman", adjust="none")
cmt <-cor$r
head(cmt)
pmt <- cor$p
head(pmt)
write.csv(cmt,"cmt.csv")
write.csv(pmt,"pmt.csv")

cmt<-as.data.frame(cmt)
cmt1<-subset(cmt,cmt$risk_score>0.3 |  cmt$risk_score<(-0.3))
dim(cmt1)
head(cmt1)
write.csv(cmt1,"cmt1.csv")

Correlation<-cmt1
Correlation$Correlation<-ifelse(cmt1$risk_score>0.3,"Positive_correlation","Negative_correlation" )
Correlation<-Correlation[,-1]
head(Correlation)
names(Correlation)<-row.names(cmt1)
write.csv(Correlation,"Correlation.csv" )

genes<-row.names(cmt1)
head(genes)
M4<-M1[,genes]
head(M4)[,1:6]
dim(M4)
write.csv(M4,"exp_selcect.csv")


rm(list = ls())
setwd("")
exp_genes<-read.csv("exp_selcect.csv",header=TRUE,row.names = 1)
head(exp_genes)[,1:6]
genes<-colnames(exp_genes)
head(genes)
genes<-as.data.frame(genes)
colnames(genes)[1]<-"symbol"
head(genes)
library(dplyr)
library(clusterProfiler)
library(topGO)
s2e <- bitr(unique(genes$symbol), fromType = "SYMBOL",  
            toType = c( "ENTREZID"),
            OrgDb = org.Hs.eg.db)
head(s2e)
head(genes)
deg <- inner_join(genes,s2e,by=c("symbol"="SYMBOL"))
head(deg)
dim(deg)
gene_ENTREZID<-deg$ENTREZID

library(ggplot2)
library(ggpubr)
#GO
ego<-enrichGO(
  gene = gene_ENTREZID,     keyType = "ENTREZID",   
  OrgDb = org.Hs.eg.db,     ont = "all",           
  pAdjustMethod = "BH",    pvalueCutoff = 0.5,      readable = TRUE
)
head(ego)
pdf("GO_bar.pdf",width = 10,height = 9)
b1<-barplot(ego,
            drop=TRUE,
            showCategory = 10,   
            split="ONTOLOGY")+
  facet_grid(ONTOLOGY~.,scales = "free")
b1
dev.off()

##KEGG
ekegg<-enrichKEGG(gene = gene_ENTREZID,
                  keyType = "kegg", 
                  organism  = "hsa",
                  pvalueCutoff = 0.5,
                  pAdjustMethod = "BH")
head(ekegg)
pdf(file = "KEGG_bar.pdf",width = 9,height = 8)
b4<-barplot(ekegg,showCategory = 10)
b4
dev.off()


rm(list = ls())
exp<-read.csv("exp_selcect.csv",header=TRUE,row.names = 1)
head(exp)[1:6]
exp1<-t(exp)
exp1[1:6,1:6]

risk<-read.csv("Clinical_km2.csv",header=TRUE,row.names = 1)
head(risk)
risk1<-risk[,1]
head(risk1)
names(risk1)<-row.names(risk)
head(risk1)

#limma 
options(stringsAsFactors = F)  
library(limma)  
design <- model.matrix(~0+factor(risk1))
colnames(design)=levels(factor(risk1))
rownames(design)=colnames(exp1)
design  

contrast.matrix<-makeContrasts(paste0(c("high","low"),collapse = "-"),levels = design)
contrast.matrix

fit <- lmFit(exp1,design)
fit<- contrasts.fit(fit, contrast.matrix) 
fit <- eBayes(fit)  ## default no trend !!!
tempOutput = topTable(fit, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
head(nrDEG)

library(dplyr)
deg = nrDEG

library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
#gsea
geneList=deg$logFC
names(geneList)=deg$ENTREZID
geneList=sort(geneList,decreasing = T)
geneList[1:6]

#GSEA——KEGG
KEGG_gseresult <- gseKEGG(geneList, nPerm = 1000, minGSSize = 10, #organism     = 'hsa',
                          maxGSSize = 1000, pvalueCutoff=1)

#GSEA——GO
Go_gseresult <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", 
                      ont="all", nPerm = 1000, minGSSize = 10, 
                      maxGSSize = 1000, pvalueCutoff=1)
Go_gseresult
#保存文件
write.table (Go_gseresult, file ="Go_gseresult.csv", sep =",", row.names =TRUE)
write.table (KEGG_gseresult, file ="KEGG_gseresult.csv", sep =",", row.names =TRUE)







