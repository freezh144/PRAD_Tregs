rm(list = ls())
setwd("")
#Lasso
multi_expr_cli<-read.csv("select.csv",header = T,row.names = 1)
head(multi_expr_cli)[,1:6]
colnames(multi_expr_cli)[1:2]<-c("os_status","os_followup" )
dim(multi_expr_cli)
set.seed(123456)
library(glmnet)
library(rms)
library(survivalROC)
library(plotROC)
library(survival)   library(survminer)
x<-as.matrix(multi_expr_cli[,3:ncol(multi_expr_cli)])
head(x)
head(multi_expr_cli)[,1:6]
y<-multi_expr_cli[,c("os_followup","os_status")]
head(y)

#Lasso regression
fit<-glmnet(x,y,family="cox")
plot(fit,label=TRUE)   
lasso_fit<-cv.glmnet(x,y,family="cox",type.measure = "deviance", nfolds = 10)
plot(lasso_fit,xvar = 'lambda',label = T) 

head(multi_expr_cli)[,1:6]
multi_expr_cli_sig<-multi_expr_cli[,sig_gene_multi_cox]
multi_expr_cli_sig<-cbind(multi_expr_cli[,1:2],multi_expr_cli_sig)
head(multi_expr_cli_sig)[,1:6]
dim(multi_expr_cli_sig)
res_sig<-multi_expr_cli_sig
dim(res_sig)

#cox regression
multi_res.cox<-coxph(Surv(os_followup,os_status)~.,data = res_sig)
x<-summary(multi_res.cox)
pvalue<-signif(as.matrix(x$coefficients)[,5],2)
head(pvalue)
HR<-signif(as.matrix(x$coefficients)[,2],2)
head(HR)
low<-signif(x$conf.int[,3],2)
high<-signif(x$conf.int[,4],2)
multi_res<-data.frame(p.value=pvalue,
                      HR=paste(HR,"(",low,"-",high,")",sep = ""),
                      stringsAsFactors = F)
head(multi_res)
dim(multi_res)
head(multi_res_sig)
dim(multi_res_sig)
multi_res_sig_gene<-row.names(multi_res_sig) 
multi_res_sig_gene
dim(multi_res_sig)
res_cox_sig<-multi_expr_cli[,multi_res_sig_gene]
head(res_cox_sig)[,1:6]
dim(res_cox_sig)
summary(res_sig$os_followup)

head(res_sig)[,1:6]
dim(res_sig)
multi_res.cox<-coxph(Surv(os_followup,os_status)~.,data = res_sig)

#risk score
coefficients<-multi_res.cox$coefficients   
coefficients[1:6]
length(coefficients)
res_sig_gene<-t(res_sig[,-1:-2])   
head(res_sig_gene)[,1:6]
dim(res_sig_gene)
head(res_sig_gene)[,1:6]
multi_res_sig_gene
coefficients1<-coefficients[multi_res_sig_gene]
coefficients1
res_sig_gene1<-res_sig_gene[multi_res_sig_gene,]

head(risk)[,1:6]
risk_score<-apply(risk, 2, sum)
risk_score[1:6]
risk_score_median<-median(risk_score)
risk_score<-as.data.frame(risk_score)
risk_score$high_low<-ifelse(risk_score$risk_score>risk_score_median,"high","low")  
head(risk_score)
risk_score_1<-cbind(res_sig[,1:2],risk_score)
head(risk_score_1)

library(ggplot2)
library(ggsci)
library(ggpubr)
risk_score_boxplot<-ggplot(risk_score_1,aes(x=high_low,y=risk_score,fill=high_low))+
  geom_boxplot(alpha=0.7)+theme_pubr(legend = "top")+
  scale_fill_manual(values = c("#FF1080","#0237E2" ))+
  #scale_fill_uchicago(alpha = 0.8,name="Strata")+
  xlab("")+ylab("Risk score")+labs(fill = "Strata")
risk_score_boxplot
ggview::ggview(risk_score_boxplot,width = 6,height = 6)
ggsave("risk_score_boxplot.pdf",width = 6,height = 6)

risk_score.cox<-coxph(Surv(os_followup,os_status)~risk_score,data = risk_score_1)
x<-summary(risk_score.cox)
x

fit_high_low<-survfit(Surv(os_followup,os_status)~high_low,data=risk_score_1)
p_high_low<-ggsurvplot(fit_high_low,linetype = "strata", ylim=c(0.35,1 ),
                       pval = TRUE,pval.coord = c(0, 0.4), #set pvalue position
                       risk.table = F,conf.int = F,
                       palette =  c("#FF1080","#0237E2" ),xlab="Time (year)",
                       legend.labs=c("high","low") #,legend = c(0.8, 0.2)
                       )
print(p_high_low)
ggview::ggview(p_high_low$plot,width = 6,height = 6)
