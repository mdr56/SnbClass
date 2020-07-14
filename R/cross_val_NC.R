#' @title Cross-validation for sNBLDA_glmsc model
#'
#' @description Uses cross-validation to select the best tuning paramter combination
#'
#' @param dat Data matrix of N by P where N is the number of subjects and P is the number of features
#'
#' @param y class label (numeric)
#' @param nfold Number of folds to be considered
#' @param lambda1 Tuning parameter vector for genes
#' @param lambda2 Tuning parameter vector for covariates
#' @param clinical matrix of clinical information
#' @param ncpu Number of CPUS to be parallel
#' @param flds list of fld sets
#' @export cross_val_NC
#'
cross_val_NC<-function(dat,y,clinical,nfold=10,lambda1=0,lambda2=0,ncpu=4,flds){
  lambda<-expand.grid(lambda1,lambda2)
  lambda_ind<-c(1:dim(lambda)[1])
  #  flds <- createFolds(y, k = nfold, list = TRUE, returnTrain = FALSE)

  wrapper_function<-function(lambda_ind){
    err<-inform_genes<-c()

    for(i in 1:nfold){

      daTr<-dat[-flds[[i]],]
      ytr<-y[-flds[[i]]]
      clinical_tr<-clinical[-flds[[i]],]
      daTe<-dat[flds[[i]],]
      true_labe<-y[flds[[i]]]
      clinical_test<-clinical[flds[[i]],]
      model<-new_sBETA_hatNC(dat=daTr,y=ytr,clinical=clinical_tr,lambda1=lambda[lambda_ind,1],lambda2=lambda[lambda_ind,2])
      cross_xte<-list('dat'=daTe,'clinical'=clinical_test)
      ypred<-Pred_function_NC(model=model,xte=cross_xte)$ypred
      err[i]<-1-mean(ypred==true_labe)
      inform_genes[i]<-sum(apply(model[[1]],2,var)>0)
      print(i)
    }
    list('lambda1'=lambda[lambda_ind,1],'lambda2'=lambda[lambda_ind,2],'err'=mean(err),'inform'=mean(inform_genes),'model'=model)

  }
  sfInit(parallel=T,cpus=ncpu)

  sfExport("dat")
  sfLibrary("SnbClass")
  sfExport('clinical')
  sfExport('lambda')
  sfExport('y')
  sfExport('flds')

  result<-sfLapply(lambda_ind,wrapper_function)
  sfStop()
  return(result)

}

