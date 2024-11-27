rm(list = ls())
setwd("")
library(pheatmap)  
library(ggplot2)
cos_exp<-read.csv("res_cox_sig.csv",header = T,row.names = 1)
head(cos_exp)
cos_exp<-t(cos_exp)
head(cos_exp)[,1:6]

col<-read.csv("Clinical_km2.csv",header = T,row.names = 1)
head(col)

p1<-pheatmap(cos_exp,annotation_col = col,
             scale = "row",cluster_cols =F, #fontsize_row = 8,
             show_colnames = F, clustering_method="ward.D2",
             #clustering_distance_cols = "correlation",
             treeheight_row = 0,treeheight_col = 0,#cutree_cols = 0,
             color=colorRampPalette(rev(c("#FF5126","white","#26A7FF")))(102)
             )
p1
print(p1,newpage=FALSE)
ggsave("heatmap.pdf",width =6,height =6)
ggsave("heatmap.tiff",width =6,height =6)




