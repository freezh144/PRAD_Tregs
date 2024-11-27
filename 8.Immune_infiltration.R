rm(list = ls())
library(utils)
library(estimate)
setwd("")
dat<-read.csv("survival_exp1.csv",header = T,row.names = 1)
library(estimate)
estimate <- function(dat){
  filterCommonGenes(input.f=input.f,
                    output.f=output.f ,
                    id="GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds=output.ds,
                platform="illumina")   ## 
  scores=read.table(output.ds,skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  return(scores)
}
pro=''
scores=estimate(dat,pro)

head(scores)
#Tumour purity=cos (0.6049872018+0.0001467884 × ESTIMATE score)

TumorPurity = cos(0.6049872018+0.0001467884 * scores[,3])   
head(TumorPurity)

scores_ALL<-cbind(scores,TumorPurity)
head(scores_ALL)
write.csv(scores_ALL,"scores_estimate.csv")

library(ggplot2)
library(ggpubr)
library(ggsci)
head(scores_ALL)
scores_ALL<-as.data.frame(scores_ALL)
scores_p1<-ggplot(scores_ALL)+
  geom_point(aes(x=StromalScore,y=TumorPurity),color="#78C25E")+
  theme_pubr() #+ylim(0.40,1)
scores_p1

scores_p2<-ggplot(scores_ALL)+
  geom_point(aes(x=ImmuneScore,y=TumorPurity),color="#DE9F54")+
  theme_pubr() #+ylim(0.40,1)
scores_p2

library(patchwork)
library(ggview)
scores_p<-scores_p1+scores_p2+plot_annotation(tag_levels = 'A')
scores_p
ggview(scores_p,width = 10,height = 6)
ggsave("scores_p.pdf",width = 10,height = 6)

head(scores_ALL)
Risk_score<-read.csv("Risk_score_km.csv",header = T,row.names = 1)
head(Risk_score)
scores_ALL_risk<-cbind(Risk_score,scores_ALL )
head(scores_ALL_risk)
write.csv(scores_ALL_risk,"scores_estimate1.csv")

scores_estimate1<-read.csv("scores_estimate1.csv",header = T,row.names = 1)
head(scores_estimate1)
library(ggsci)
library(ggplot2)
library(ggpubr)
p1<-ggplot(scores_estimate1)+geom_violin(aes(x=Risk_score, y=ImmuneScore,fill=Risk_score),alpha=0.6)+
  theme_pubr(legend = "bottom")+ggsci::scale_fill_nejm(alpha = 0.7)+
  geom_boxplot(aes(x=Risk_score, y=ImmuneScore),width=0.15,outlier.colour = "NA")+
  stat_compare_means(aes(x=Risk_score,y=ImmuneScore,group=Risk_score)
                     ,hide.ns = TRUE)
p1

head(scores_estimate1)
p2<-ggplot(scores_estimate1)+geom_violin(aes(x=Risk_score, y=TumorPurity,fill=Risk_score),alpha=0.6)+
  theme_pubr(legend = "bottom")+ ggsci::scale_fill_nejm(alpha = 0.7)+
  geom_boxplot(aes(x=Risk_score, y=TumorPurity),width=0.15,outlier.colour = "NA")+
  stat_compare_means(aes(x=Risk_score,y=TumorPurity,group=Risk_score)
                     ,hide.ns = TRUE)
p2


p<-p1+p2+plot_annotation(tag_levels = 'A')
p
ggview(p,width = 10,height = 6)


p_zh<-scores_p/p+plot_annotation(tag_levels = 'A')
p_zh
ggview(p_zh,width = 12,height = 12)


