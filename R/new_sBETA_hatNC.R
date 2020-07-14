#' @title sNBLDA_GLMsC
#'
#' @description Fit the predictive model for sNBLDA_GLMsC
#'
#' @param dat Data matrix of N by P where N is the number of subjects and P is the number of features
#'
#' @param y class label (numeric)
#' @param lambda1 Tuning parameter vector for genes
#' @param lambda2 Tuning parameter vector for covariates
#' @param clinical matrix of clinical information
#' @export new_sBETA_hatNC
new_sBETA_hatNC=function(dat,y,lambda1,lambda2,clinical){

  designC=cbind(1,clinical,y)
  phi=1/estimateDisp(t(dat),designC)$tagwise.dispersion

  covlen<-dim(clinical)[2]
  norm_factors=NullModel1(dat,type="quantile",stat='mean')$sizes
  ref<-NullModel1(dat,type="quantile",stat='mean')$reference


  datTr=cbind(dat,clinical,y,norm_factors)
  library(MASS)
  featurlen=ncol(datTr)-(covlen+2)
  K=length(levels(factor(datTr[,featurlen+(covlen+1)])))
  clas_by_labe=matrix(0,1,ncol(datTr))
  Nk=c()
  for(k in 1:K){
    c_b_l=as.matrix(datTr[datTr[,(featurlen+(covlen+1))]==(k),])
    Nk[k]=nrow(c_b_l)
    clas_by_labe=rbind(clas_by_labe,c_b_l)
  }
  clas_by_labe=clas_by_labe[-1,]

  repeat{
    inibet=matrix(0,K,featurlen)
    inialp=matrix(0,covlen,featurlen)
    for(j in 1:featurlen){
      resp=log(c(datTr[,j]+1)/c(datTr[,featurlen+(covlen+2)]))
      intec=matrix(0,nrow(datTr),K)
      for(i in 1:nrow(datTr)){
        intec[i,c(datTr[i,featurlen+(covlen+1)])]=1
      }
      des=cbind(intec,datTr[,c((featurlen+1):(featurlen+covlen))])
      est=ginv(t(des)%*%des)%*%t(des)%*%resp
      inibet[,j]=est[c(1:K)]
      inialp[,j]=est[-c(1:K)]
    }
    if((sum(inibet)+sum(inialp)!="NaN")&(sum(inibet)+sum(inialp)!="NA")){break}
  }
  #inibet<-init_bet
  iteration <- 0
  precision <- 1

  s=matrix(c(clas_by_labe[,(featurlen+covlen+2)]),sum(Nk),featurlen)
  s=log(s)
  y_clea=as.matrix(clas_by_labe[,1:featurlen])
  PHI_exte=t(matrix(rep(phi,sum(Nk)),featurlen,sum(Nk)))
  pt=proc.time()
  PP=c()
  repeat{
    #
    inibet[which(abs(inibet)>=100)]=sign(inibet[which(abs(inibet)>=100)])*5
    inialp[which(abs(inialp)>=100)]=sign(inialp[which(abs(inialp)>=100)])*sample(seq(-1,1,0.1),1)
    #inibet[which(abs(inibet)>=100)]=10
    #inialp[which(abs(inialp)>=100)]=0

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
    cova=as.matrix(clas_by_labe[,(featurlen+1):(featurlen+covlen)])%*%inialp
    lnmu=s+bet+cova

    mu=exp(lnmu)
    w=(PHI_exte)/((PHI_exte/mu)+1)
    tau=(lnmu)+(y_clea/mu)-1
    bs=(caculat%*%(w*(tau-s-cova)))/(caculat%*%w)
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
    cova_clea=as.matrix(clas_by_labe[,(featurlen+1):(featurlen+covlen)])
    new_alp=matrix(0,covlen,featurlen)
    for(q in 1:covlen){
      alp_pen=matrix(0,3,featurlen)
      amo=colSums(apply(w,2,function(x){x*(cova_clea[,q]^2)}))
      w_pen<-matrix(1,nrow=dim(w)[1],ncol=dim(w)[2])
      amo1=colSums(apply(w_pen,2,function(x){x*(cova_clea[,q]^2)}))
      aso=colSums((apply(w,2,function(x){x*cova_clea[,q]}))*(tau-s-bet-as.matrix(cova_clea[,-q])%*%inialp[-q,]))
      #aso=colSums((apply(w,2,function(x){x*cova_clea[,q]}))*(as.matrix(cova_clea[,q])%*%inialp[q,])+(y_clea/mu)-1)
      alp_pen[1,]=aso/amo
      alp_pen[2,]=(aso/amo)-lambda2/amo1
      alp_pen[3,]=(aso/amo)+lambda2/amo1
      pen1=(alp_pen[1,]>0)*(alp_pen[2,]>0)
      pen2=(alp_pen[1,]<0)*(alp_pen[3,]<0)
      pen3=1-(pen1+pen2)
      new_alp[q,]=alp_pen[2,]*pen1+alp_pen[3,]*pen2+rep(0,featurlen)*pen3
    }

    precision <- sqrt((max(mean((inialp-new_alp)^2),mean((inibet-new_bet)^2))))
    PP=c(PP,precision)
    max_diff<-max(max(abs(inialp-new_alp)),max(abs(inibet-new_bet)))



    print(max_diff)
    ind_nc<-union(unique(which(abs(inialp-new_alp)>0.01,arr.ind=T)[,2]),unique(which(abs(inibet-new_bet)>0.01,arr.ind=T)[,2]))
    print(ind_nc)
    if(max_diff<0.01){break}
    print(max_diff)
    print(iteration)
    iteration <- iteration + 1
    if(iteration==80){break}
    inialp=new_alp
    inibet=new_bet

  }
  proc.time()-pt
  if(length(ind_nc)>0){
    inibet<-inibet[,-ind_nc]
    inialp<-inialp[,-ind_nc]
    phi<-phi[-ind_nc]
  }
  pri<-Nk/sum(Nk)
  OUTPU=list('beta'=inibet,'pri'=pri,'bbar'=bbar,'alpha'=inialp,'ind_nc'=ind_nc,'iteration'=iteration,'norm'=norm_factors,'phi'=phi,'reference'=ref)

  return(OUTPU)
}

