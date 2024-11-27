rm(list=ls())
#install.packages("xlsx")
library("xlsx")
library(survival)   
library(survminer)
setwd("")
expr_cli<-read.csv("Clinical.csv",header = T,row.names = 1)
head(expr_cli)

covariates<-names(expr_cli)
univ_formulas<-apply(covariates,  function(x) 
  as.formula(paste("Surv(Time,Status)~",x)) )    
head(univ_formulas)
univ_models<-lapply(univ_formulas,function(x){coxph(x,data = expr_cli)})   
head(univ_models)

#HR,p-value
univ_results<-apply(univ_models, function(x) {
  x<-summary(x)
  p.value<-signif(x$wald["pvalue"],digits = 2)  #
  HR<-signif(x$coef[2],digits = 2)   #
  res<-c(p.value,HR)
  names(res)<-c("p.value","HR (95% CI for HR)")
  return(res)
})
head(univ_results)
res<-t(as.data.frame(univ_results,check.names=FALSE))
res<-as.data.frame(res)
head(res)

head(expr_cli)[,1:6]
res_sig<-expr_cli
multi_res.cox<-coxph(Surv(Time,Status)~.,data = res_sig)
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
#rownames(multi_res)<-colnames(res_sig[,3:ncol(res_sig)])
head(multi_res)
dim(multi_res)
write.csv(multi_res,"multi_cox_result.csv")



#KM_plot
km<-read.csv("Clinical_km.csv",header = T,row.names = 1)
head(km)
genes<-names(km)[3:ncol(km)]
head(genes)
for (i in genes) {
  print(i)
  fit<-survfit(Surv(Time,Status)~km[,i],data=km)
  p<-ggsurvplot(fit,linetype = "strata",
                pval = TRUE,
                palette = "Dark2",
                #legend.labs=c(paste0(i,"=H"),paste0(i,"=L"))
  )
  pdf(paste0(i,"_surv.pdf"),width=5,height=5)
  print(p,newpage=FALSE)
  dev.off()
}
