rm(list=ls())
setwd("")
multi_cox_res_sig_gene<-read.csv("res_cox_sig.csv",header = T,row.names = 1)
head(multi_cox_res_sig_gene)[,1:6]
colnames(multi_cox_res_sig_gene)
multi_cox_res_sig_gene<-t(multi_cox_res_sig_gene)
head(multi_cox_res_sig_gene)[,1:6]
dat<-multi_cox_res_sig_gene

set.seed(1234)
#consensus clustering
d=dat
d = sweep(d,1, apply(d,1,median,na.rm=T))
library(ConsensusClusterPlus)
results = ConsensusClusterPlus(as.matrix(d),maxK=6,reps=1000,pItem=0.8,
                               pFeature=1,
                               title="output",
                               clusterAlg="pam",distance="pearson",
                               seed=1234,plot="pdf",writeTable=F)

#group  
cluster<-read.csv("result_Cluster_to_2.csv",row.names = 1,header = T)
head(cluster)
colnames(cluster)<-"Cluster"
cluster$Cluster<-gsub("1","Cluster_1",cluster$Cluster)  #
cluster$Cluster<-gsub("2","Cluster_2",cluster$Cluster)
cluster$Cluster<-gsub("3","Cluster_3",cluster$Cluster)
head(cluster)
table(cluster$Cluster)
dim(cluster)

multi_expr_cli<-read.csv("select.csv",header = T,row.names = 1)
head(multi_expr_cli)[,1:6]  
colnames(multi_expr_cli)[1:2]<-c("os_status","os_followup" )
head(multi_expr_cli)[,1:6]
dim(multi_expr_cli)
multi_expr_cli_cluster<-cbind(multi_expr_cli[,1:2],cluster)
head(multi_expr_cli_cluster)


library(ggplot2)
library(survminer)
cluster.cox<-coxph(Surv(os_followup,os_status)~Cluster,data = multi_expr_cli_cluster)
x<-summary(cluster.cox)
x
fit_cluster<-survfit(Surv(os_followup,os_status)~Cluster,data=multi_expr_cli_cluster)
fit_cluster
p_cluster<-ggsurvplot(fit_cluster,linetype = "strata",
                      pval = TRUE,pval.coord = c(0, 0.7),
                      risk.table = T,conf.int = F,ylim=c(0.6,1),
                      palette =  c("#CA645B","#877BAF" ),xlab="Time (year)" ,
                      legend.labs=c("Cluster_1","Cluster_2")  #,"Cluster_3"
)
print(p_cluster)
library(patchwork)
p<-p_cluster$plot/p_cluster$table+
  plot_layout(nrow = 2,#
              heights = c(3.5, 1))#
p
ggview::ggview(p,width = 6,height =7)




