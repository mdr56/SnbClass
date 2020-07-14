#' @title prediction for a new observation using sNBLDA_glmsc model
#'
#' @description Predicts the class for a new observation
#'
#' @param xte Data matrix of N by P where N is the number of subjects and P is the number of features
#'
#' @param model sNBLDA_glm model built using training data
#' @export Pred_function_N
Pred_function_NC<-function(model,xte){
  OUTPUT<-model

  sbethat=as.matrix(OUTPUT[[1]])
  salphat=as.matrix(OUTPUT[[4]])
  chec=apply(sbethat,2,var)
  selD=c(which(chec>0))
  sbethatc=as.matrix(sbethat[,selD])
  salphatc=as.matrix(salphat[,selD])
  phic=OUTPUT$phi[selD]
  pri<-OUTPUT$pri
  K<-dim(sbethat)[1]
  Tesampfac<- NullModelTest1(ref=OUTPUT$reference,xte=xte$dat,type='quantile')$sizeste
  if(length(OUTPUT$ind_nc)>0){
    xte$dat<-xte$dat[,-OUTPUT$ind_nc]
  }
  datTe1=as.matrix(xte$dat[,selD])

  Covte=xte$clinical #salphatc
  scores<-matrix(,nrow=length(Tesampfac),ncol=K)
  fitlab=c()
  for(l in 1:length(Tesampfac)){
    pk=c()
    for(k in 1:K){
      pp=t(as.matrix(datTe1[l,]))%*%(as.matrix(sbethatc[k,]+(Covte%*%salphatc)[l,]-log(brob(log(Tesampfac[l])+sbethatc[k,]+(Covte%*%salphatc)[l,])+phic)))-sum(phic*log(brob(log(Tesampfac[l])+sbethatc[k,]+(Covte%*%salphatc)[l,])+phic))
      +log(pri[k])
      pk[k]=as.numeric(pp)
    }
    scores[l,]<-pk
    fitlab[l]=which(pk==max(pk))
  }

  pred_result<-list('ypred'=fitlab,'scores'=scores)
  return(pred_result)
}
