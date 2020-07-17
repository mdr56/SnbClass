# Project

This contains the package and code scripts for the following paper,
Rahman, T., Huang, H. E., Li, Y., Tai, A. S., Hsieh, W. P., & Tseng, G. (2019). A sparse negative binomial classifier with covariate adjustment for RNA-seq data. BioRxiv, 636340.

### Installing

For installation, we recommend to unzip the tar.gz file first and then use devtools::install() to install the package, which can make sure to install all the depends. 

##
## The following shows how to get result for one run of simulation.

##Simulate the data 

```
library(SnbClass)
sim<-Simulat_data(N_sample=1200,pareto_shape=8,pareto_scale=10,lfc_mean=0.25,prob_inform = 0.1250,cov_mean =0.25,N_train=200)

## For this version the folds for cross-validation is found from "PoiClaClu" package. Next version will automatically determine the folds.
plda_list<-Classify.cv(sim$train_data$dat,y=sim$train_data$ytr, rhos = seq(0,6,len=31), beta = 1, nfolds = 10, type="quantile",
                       folds = NULL, transform=TRUE)
plda_folds<-plda_list$folds
```
##Running the snbGLM cross-validations, full model and prediction for testing data

```

result<-cross_val_N(dat=sim$train_data$dat,y=sim$train_data$ytr,lambda1=0,ncpu = 10,flds=plda_folds)

cross_result<-data.frame(do.call(rbind,lapply(result,function(x) do.call(cbind, x[1:3]))))

opt_lambda<-cross_result$lambda[which.min(cross_result$err)]
Model<-new_sBETA_hatN(datTr = sim$train_data$dat,y=sim$train_data$ytr,lambda1 = seq(0,10,1))
ypred<-Pred_function_N(model=Model,xte=sim$test_data$dat)
test_acc<-mean(ypred$ypred==sim$test_data$yte)



```
##SnbGLM_sC


```
result_glmsc<-cross_val_NC(dat=sim$train_data$dat,y=sim$train_data$ytr,clinical=sim$train_data$clinical,lambda1=seq(0,5,1),lambda2=seq(0,5,1),ncpu = 10,flds=plda_folds)
cross_result_glmsc<-data.frame(do.call(rbind,lapply(result_glmsc,function(x) do.call(cbind, x[1:4]))))
opt_lambda1<-cross_result_glmsc$lambda1[which.min(cross_result_glmsc$err)]
opt_lambda2<-cross_result_glmsc$lambda2[which.min(cross_result_glmsc$err)]

model_glmsc<-new_sBETA_hatNC(dat=sim$train_data$dat,y=sim$train_data$ytr,lambda1=0,lambda2=0,clinical=sim$train_data$clinical)

pred_glmsc<-Pred_function_NC(model=model_glmsc,xte=sim$test_data)

test_acc_glmsc<-mean(pred_glmsc$ypred==sim$test_data$yte)


```
