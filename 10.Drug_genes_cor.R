rm(list = ls())
setwd("")
library(psych)
library(pheatmap)
library(ggplot2)
library(reshape2)
gene<-read.csv("select.csv",header=TRUE,row.names = 1)
head(gene)[,1:5]
drug<-read.csv("DrugPredictions_select.csv",header=TRUE,row.names = 1)
head(drug)[,1:6]
colnames(drug)

cor = corr.test(gene, drug, method="spearman", adjust="none")
cmt <-cor$r
head(cmt)
pmt <- cor$p
head(pmt)
write.csv(cmt,"cmt.csv")
write.csv(pmt,"pmt.csv")

#cmt[abs(cmt)<0.3]<-0
if (!is.null(pmt)){
  
  ssmt <- pmt< 0.001
  
  pmt[ssmt] <-'***'
  smt <- pmt >0.001& pmt <0.01
  
  pmt[smt] <- '**'
  
  smt <- pmt >0.01& pmt <0.05
  
  pmt[smt] <- '*'
  
  pmt[!ssmt&!smt]<- ''
  
} else {
  
  pmt <- F
  
}
mycol<-colorRampPalette(c("#107AB0","#D0FEFE","#FD5956"))(800)
mycol1<-colorRampPalette(c("#107AB0","#D0FEFE","#FD5956"))(800)

p1<-pheatmap(t(cmt),scale = "none", cluster_row = T, cluster_col = T, 
             treeheight_row = 0,treeheight_col = 0,
         border=NA, #fontsize_row=8, fontsize_col=8,
         display_numbers = t(pmt), 
         fontsize_number = 10, number_color = "black", 
         cellwidth = 13, 
         cellheight =13,color=mycol)
p1$gtable
library(ggplotify)
p1<-as.ggplot(p1)
p1
ggview::ggview(p1,width = 6,height = 6)
ggsave("HEATPLOT.pdf",width = 7,height = 6)


library(reshape2)
cmt1<-melt(cmt)
colnames(cmt1)<-c("genes","drug","cor")
head(cmt1)
write.csv(cmt1,"cmt_correlation.csv")

pmt1<-melt(pmt)   #p_value
colnames(pmt1)<-c("genes","drug","p_value")
head(pmt1)
gd<-cbind(cmt1,pmt1$p_value)
colnames(gd)<-c("genes","drug","cor","p_value")
head(gd)
write.csv(gd,"gd.csv")

library(ggplot2)
library(ggpubr)
head(gd)
gd$Cor<-abs(gd$cor)
gd$type<-ifelse(gd$cor>0,"Postive","Negative")
p2<-ggplot(gd,aes(x=genes,y=drug))+
  geom_point(aes(color=cor,size=Cor,shape=type))+
  scale_shape_manual(values = c(15,16))+
  scale_color_gradient2(low = "#3B4992FF",mid="#E0DEC5",high = "#BB0021FF")+
  geom_text(label= round(gd$cor,2))+xlab("")+ylab("")+
  theme_light()+theme(legend.position = "bottom")#+guides(shape=TRUE)
p2
ggview::ggview(p2,width = 8,height =6)
ggsave("PLOT.pdf",width = 8,height = 6)






