rm(list = ls())
setwd("")
#Import expression matrix
expr_cli<-read.csv("survival_exp.csv",header = T,row.names = 1)
head(expr_cli)[,1:6]
dim(expr_cli)

#survival analysis
library(survival)   #Load Package
library(survminer)
covariates<-names(expr_cli)[3:ncol(expr_cli)]   
univ_formulas<-apply(covariates,  function(x) 
  as.formula(paste("Surv(os_followup,os_status)~",x)) )     
head(univ_formulas)
#loop analysis
univ_models<-apply(univ_formulas,function(x){coxph(x,data = expr_cli)})   
head(univ_models)

#Extract HR and p-value
univ_results<-apply(univ_models, function(x) {
  x<-summary(x)
  p.value<-signif(x$wald["pvalue"],digits = 2)  
  HR<-signif(x$coef[2],digits = 2)   
  res<-c(p.value,HR)
  names(res)<-c("p.value","HR (95% CI for HR)")
  return(res)
})
head(univ_results)
res<-t(as.data.frame(univ_results,check.names=FALSE))
res<-as.data.frame(res)
head(res)

res1<-res[order(res$p.value,decreasing = T),]
head(res1)
dim(res)
write.csv(res,"univariate_cox_result.csv")








