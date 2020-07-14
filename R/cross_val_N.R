#' @title Cross-validation for sNBLDA_glm model
#'
#' @description Uses cross-validation to select the best tuning paramter
#'
#' @param dat Data matrix of N by P where N is the number of subjects and P is the number of features
#'
#' @param y class label (numeric)
#' @param nfold Number of folds to be considered
#' @param lambda1 Tuning parameter vector
#' @param ncpu Number of CPUS to be parallel
#' @param flds list of fld sets
#' @export cross_val_N

cross_val_N<-function(dat,y,nfold=10,lambda1=0,ncpu=4,flds){

  #  flds <- createFolds(y, k = nfold, list = TRUE, returnTrain = FALSE)

  wrapper_function<-function(lambda1){
    err<-inform_genes<-c()

    for(i in 1:nfold){

      daTr<-dat[-flds[[i]],]
      ytr<-y[-flds[[i]]]
      daTe<-dat[flds[[i]],]
      true_labe<-y[flds[[i]]]

      model<-new_sBETA_hatN(datTr=daTr,y=ytr,lambda1=lambda1)
      ypred<-Pred_function_N(model=model,xte=daTe)$ypred
      err[i]<-1-mean(ypred==true_labe)
      inform_genes[i]<-sum(apply(model[[1]],2,var)>0)
      print(i)
    }
    list('lambda'=lambda1,'err'=mean(err),'inform'=mean(inform_genes),'model'=model)

  }
  sfInit(parallel=T,cpus=ncpu)

  sfExport("dat")
  sfSource("function_codes.R")
  #sfExport('lambd')
  sfExport('y')
  sfExport('flds')

  result<-sfLapply(lambda1,wrapper_function)
  sfStop()
  return(result)

}


cross_val_N<-function(dat,y,nfold=10,lambda1=0,ncpu=4,flds){

  #  flds <- createFolds(y, k = nfold, list = TRUE, returnTrain = FALSE)

  wrapper_function<-function(lambda1){
    err<-inform_genes<-c()

    for(i in 1:nfold){

      daTr<-dat[-flds[[i]],]
      ytr<-y[-flds[[i]]]
      daTe<-dat[flds[[i]],]
      true_labe<-y[flds[[i]]]

      model<-new_sBETA_hatN(datTr=daTr,y=ytr,lambda1=lambda1)
      ypred<-Pred_function_N(model=model,xte=daTe)$ypred
      err[i]<-1-mean(ypred==true_labe)
      inform_genes[i]<-sum(apply(model[[1]],2,var)>0)
      print(i)
    }
    list('lambda'=lambda1,'err'=mean(err),'inform'=mean(inform_genes),'model'=model)

  }
  sfInit(parallel=T,cpus=ncpu)

  sfExport("dat")
  sfSource("function_codes.R")
  #sfExport('lambd')
  sfExport('y')
  sfExport('flds')

  result<-sfLapply(lambda1,wrapper_function)
  sfStop()
  return(result)

}

