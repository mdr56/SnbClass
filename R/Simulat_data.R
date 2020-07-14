

#' @title Simulate count data
#'
#' @description Simulate RNA-seq data with clinical variables
#'
#' @param N_sample Total number of samples to simulate
#' @param pareto_shape Value of the shape paramter
#' @param pareto_scale Value of the scale paramter
#' @param lfc_mean logfold change
#' @param cov_mean Mean effect size of the covariates
#' @param prob_inform Number of genes affected by clinical variables
#' @param N_train Number of training samples
#' @export Simulat_data




Simulat_data<-function(N_sample=1100,pareto_shape=8,pareto_scale=10,lfc_mean=0.5,cov_mean=0,prob_inform=0,N_train=100){

  n=N_sample
  phi=1/c(rpareto(1000,shape=pareto_shape,scale=pareto_scale)) #1 0.5 0.25
  rho=lfc_mean
  eta=cov_mean
  s=c(runif(n, min = 0.75, max = 1.25))
  b=c(abs(log(rowMeans(dat_mice_ctrl[[1]])[sample(1:10280,300,replace = F)]*(1)))) #baseline
  b=rbind(b,b,b)
  delset=cbind(c(-1,0,1),c(1,0,1),c(-1,0,-1))
  Del=c(rtruncnorm(ncol(b),a=rho/2,b=Inf,mean=rho, sd = 0.1))
  del=apply(delset[,sample(1:ncol(delset),ncol(b),replace=T)],2
            ,function(x){x[sample(1:nrow(delset),nrow(delset),replace=F)]})
  B=b+del*matrix(rep(Del,nrow(del)),nrow(del),ncol(del),byro=T)

  sex=c(scale(sample(1:3,prob=c(0.334,0.333,0.333),size=n,replace = T)))#c(scale(rbinom(n,1,0.5)))
  age=c(scale(rgamma(n,scale=5, shape=10)))
  eps=matrix(c(rtruncnorm(2*1000,a=eta/2,b=Inf,mean=eta, sd = 0.1)*sample(c(-1,1),size=2000,replace=T)),2,1000)
  gamset=cbind(c(1,1),c(1,0),c(0,1),c(0,0))

  # gam=apply(gamset[,sample(1:6,ncol(eps),replace=T,prob=c(rep((1-prop_noninform)/5,5),prob_noninform)],2
  #           ,function(x){x[sample(1:nrow(gamset),nrow(gamset),replace=F)]})

  gam<-gamset[,sample(1:4,1000,replace=T,c(rep((prob_inform)/3,3),1-prob_inform))]
  A=cbind(gam*eps)
  covar=cbind(sex,age)%*%A

  datori=matrix(0,n,ncol(b))
  y<-sample(c(1,2,3),size=n,replace=T)

  for(i in 1:n){
    for(j in 1:ncol(b)){
      lnMu=log(s[i])+B[y[i],j] +covar[i,j]
      datori[i,j]=rnbinom(1,size=phi[j], mu=exp(lnMu))
    }
  }


  # y=c(rep(1,(n/3)),rep(2,(n/3)),rep(3,(n/3)))
  nois =matrix(0,n,700)
  bno=c(abs(log(rowMeans(dat_mice_ctrl[[1]])[sample(1:10280,700,replace = F)])))
  for(i in 1:n){
    for(j in 1:700){
      lnnoismu=log(s[i])+bno[j] +covar[i,(ncol(b)+j)]
      nois[i,j]=rnbinom(1,size=phi[(ncol(b)+j)], mu=exp(lnnoismu))
    }
  }
  dat=cbind(datori,nois)
  clinical<-cbind(sex,age)
  #s=NullModel(x,type="quantile")$sizes*n
  #s=c(calcNormFactors(t(x),method="TMM"))
  #dat=cbind(x,sex,age,y,s)
  Bnois<-rbind(bno,bno,bno)
  B_true<-cbind(B,Bnois)
  #design=cbind(1,y)
  #phi=1/estimateDisp(t(x),design)$tagwise.dispersion

  #  designC=cbind(1,sex,age,y)
  # phiC=1/estimateDisp(t(x),designC)$tagwise.dispersion

  #seperate
  ord=sample(1:nrow(dat),nrow(dat),replace=FALSE) #MESS
  dat=dat[ord,]
  y=y[ord]
  clinical=clinical[ord,]
  datTr=dat[c(1:N_train),]
  datTe=dat[-c(1:N_train),]
  ytr<-y[c(1:N_train)]
  yte<-y[-c(1:N_train)]
  clinical_tr<-clinical[c(1:N_train),]
  clinical_test<-clinical[-c(1:N_train),]
  trualp=matrix(1,2,1000)
  trualp[which(A==0)]=0
  trualp=trualp[,c(1:300)]

  train_data<-list('dat'=datTr,'ytr'=ytr,'clinical'=clinical_tr)
  test_data<-list('dat'=datTe,'yte'=yte,'clinical'=clinical_test)

  list("train_data"=train_data,'test_data'=test_data,'B_true'=B_true,'A'=A,'s'=s,'phi'=phi)

}
