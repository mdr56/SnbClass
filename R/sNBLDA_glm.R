#' @title sNBLDA_glm model
#'
#' @description Builds the predictive model
#'
#' @param dat Data matrix of N by P where N is the number of subjects and P is the number of features
#'
#' @param y class label (numeric)
#' @param lambda1 Tuning parameter vector
#' @export new_sBETA_hatN
new_sBETA_hatN=function(datTr,y,lambda1){
  library(MASS)

  design=cbind(1,y)
  phi=1/estimateDisp(t(datTr),design)$tagwise.dispersion

  featurlen=ncol(datTr)
  K=length(levels(factor(y)))
  clas_by_labe=matrix(0,1,ncol(datTr))
  Nk=c()
  for(k in 1:K){
    c_b_l=as.matrix(datTr[y==(k),])
    Nk[k]=nrow(c_b_l)
    clas_by_labe=rbind(clas_by_labe,c_b_l)
  }
  clas_by_labe=clas_by_labe[-1,]

  norm_factors=NullModel1(datTr,type="quantile",stat='mean')$sizes
  ref<-NullModel1(datTr,type="quantile",stat='mean')$reference
  norm_factors1=NullModel1(datTr,type="quantile",stat='mean')$sizes

  repeat{
    inibet=matrix(0,K,featurlen)
    for(j in 1:featurlen){
      resp=log(c(datTr[,j]+1))/norm_factors
      intec=matrix(0,nrow(datTr),K)
      for(i in 1:nrow(datTr)){
        intec[i,y[i]]=1
      }
      des=intec
      est=ginv(t(des)%*%des)%*%t(des)%*%resp
      inibet[,j]=est

    }
    if((sum(inibet)!="NaN")&(sum(inibet)!="NA")){break}
  }
  #
  iteration <- 0
  precision <- 1
  norm_factors=NullModel1(clas_by_labe,type="quantile",stat='mean')$sizes
  s=matrix(norm_factors,sum(Nk),featurlen)
  s=log(s)
  y_clea=as.matrix(clas_by_labe[,1:featurlen])
  PHI_exte=t(matrix(phi,featurlen,sum(Nk)))
  pt=proc.time()
  PP=c()
  repeat{
    bet=matrix(0,1,featurlen)
    caculat=matrix(0,K,sum(Nk))
    for(k in 1:K){
      if(k==1){
        caculat[k,]=rep(c(1,0),c(Nk[k],sum(Nk)-Nk[k]))
      } else{
        befsum=sum(Nk[1:(k-1)])
        caculat[k,]=rep(c(0,1,0),c(befsum,Nk[k],sum(Nk)-befsum-Nk[k]))
      }
      bet=rbind(bet,t(matrix(c(inibet[k,]),featurlen,Nk[k])))
    }
    bet=bet[-1,]
    lnmu=s+bet
    mu=exp(lnmu)
    w=(PHI_exte)/((PHI_exte/mu)+1)
    tau=(lnmu)+(y_clea/mu)-1
    bs=(caculat%*%(w*(tau-s)))/(caculat%*%w)
    #bs=(caculat%*%(w*(bet+(y_clea/mu)-1)))/(caculat%*%w)
    bbar=colMeans(bs) #1 x j
    bm=(caculat%*%w)
    bet_pen=(rbind(bs,bs-((K-1)*lambda1/K),bs+((K-1)*lambda1/K))) #K*3 x j

    new_bet=matrix(0,K,featurlen)
    for(k in 1:K){
      pen1=(bet_pen[k,]>bbar)*(bet_pen[K+k,]>bbar)
      pen2=(bet_pen[k,]<bbar)*(bet_pen[2*K+k,]<bbar)
      pen3=1-(pen1+pen2)
      new_bet[k,]=bet_pen[K+k,]*pen1+bet_pen[2*K+k,]*pen2+bbar*pen3
    }
    precision <- sqrt((mean((inibet-new_bet)^2)))
    PP=c(PP,precision)
    max_diff<-max(abs(inibet-new_bet))
    inibet=new_bet
    if(max_diff<0.001){break}
    iteration <- iteration + 1
    if(iteration==15000){break}
  }
  proc.time()-pt
  pri<-Nk/sum(Nk)

  OUTPU=list('inibet'=inibet,'bbar'=bbar,'pri'=pri,'norm_factors'=norm_factors1,'phi'=phi,'ref'=ref,'iteration'=iteration)

  return(OUTPU)
}
