rm(list = ls())
#CIBERSORT 
setwd("")
library(rtracklayer)

exp<-read.csv("survival_exp1.csv",header=TRUE,row.names = 1)  
dim(exp)
head(exp)[,1:6]
write.table(exp, "expression.txt", sep="\t")

source("CIBERSORT.R")
# Define LM22 file
LM22.file <- "LM22.txt"
exp.file <- "expression.txt"

TME.results = CIBERSORT(LM22.file, exp.file, perm = 1000, QN = TRUE)

# output CIBERSORT results
write.table(TME.results, "TME.results.output.txt", 
            sep = "\t", row.names = T, col.names = T, quote = F)

save(TME.results,file = "ciber_CHOL.Rdata")
load("ciber_CHOL.Rdata")
TME.results[1:4,1:4]
dim(TME.results)
View(TME.results)
re <- TME.results[,-(23:25)]  

library(pheatmap)
k <- apply(re,2,function(x) {sum(x == 0) < nrow(TME.results)/2})
table(k)
re2 <- as.data.frame(t(re[,k]))
head(re2)
write.csv(re2,"res.csv")
re2<-read.csv("res.csv",header = T,row.names = 1)

dim(re2)
#an = data.frame(group = B1$group,
              #  row.names = colnames(exp))
Risk.score<-read.csv("Risk_score_km.csv",header=TRUE,row.names = 1)  
head(Risk.score)

p1<-pheatmap(re2,scale = "row",
             show_colnames = F,cluster_cols = F,
             annotation_col = Risk.score,treeheight_row = 0,
             color = colorRampPalette(c("#01A187", "#EFE9D9", "#DC0100"))(50))
p1

library(magrittr)
library(tidyr)
library(dplyr)
library(data.table)
library(tibble)
library(ggsci)
library(ggplot2)
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
dat <- re %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)
head(dat)
p2<-ggplot(dat,aes(Sample,Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(22))
p2

library(stringr )
head(dat)
Risk.score<-read.csv("Risk_score_km.csv",header=TRUE)  
head(Risk.score)
dat1<-cbind(dat,Risk.score)
head(dat1)
dat2<-dat1[,-4]
head(dat2)

library(ggpubr)
p4<-ggplot(dat2,aes(Cell_type,Proportion,fill = Risk_score)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_pubr() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "bottom") + coord_flip()+
  #theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = mypalette(22)[c(6,1)])+ 
  stat_compare_means(aes(group = Risk_score,label = ..p.signif..),method = "kruskal.test")
p4
ggview::ggview(p4,width = 10, height = 8)




