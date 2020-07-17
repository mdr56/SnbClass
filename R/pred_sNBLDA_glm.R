#' @title prediction for a new observation using sNBLDA_glm model
#'
#' @description Predicts the class for a new observation
#'
#' @param xte Data matrix of N by P where N is the number of subjects and P is the number of features
#'
#' @param model sNBLDA_glm model built using training data
#' @export Pred_function_N
Pred_function_N<-function(model,xte){
  OUTPUT<-model
  sbethat=as.matrix(OUTPUT[[1]])
  K<-dim(sbethat)[1]
  chec=apply(sbethat,2,var)
  selA=c(which(chec>0))

  sbethatc=as.matrix(sbethat[,selA])
  phic=OUTPUT$phi[selA]
  pri<-OUTPUT$pri

  Tesampfac<- NullModelTest1(ref=OUTPUT$ref,xte=xte,type='quantile')$sizeste
  datTe1=as.matrix(xte[,selA])
 # print(Tesampfac)
  fitlab=c()
  scores<-matrix(,nrow=length(Tesampfac),ncol=K)
  for(l in 1:length(Tesampfac)){
    pk=c()
    for(k in 1:K){
      pp=t(as.matrix(datTe1[l,]))%*%(as.matrix(sbethatc[k,]-log(brob(log(Tesampfac[l])+sbethatc[k,])+phic)))-sum(phic*log(brob(log(Tesampfac[l])+sbethatc[k,])+phic))
      +log(pri[k])
      pk[k]=as.numeric(pp)
    }
    scores[l,]<-pk
    fitlab[l]=which(pk==max(pk))
  }
  pred_result<-list('ypred'=fitlab,'scores'=scores)
  return(pred_result)
}

