mqnorm<-function(q,mu=c(0,0),Sigma=diag(2),C=diag(2)){
  ##C is the points generating cone with default R+.
  P<-mvrnorm(10000,mu=mu,Sigma=Sigma)
  P<-as.data.frame(P)
  lam<-seq(0,1,by=0.01)
  lam<-cbind(1-lam,lam)
  if (is.matrix(C)){
    if (any(C[,2]>=0)&&any(C[,2]<0)){
      C1<-C[which(C[,2]>=0),]
      C1<-matrix(C1,ncol=2,byrow=F)
      C2<-C[which(C[,2]<0),]
      C2<-matrix(C2,ncol=2,byrow=F)
      Cos1<-C1[,1]/sqrt(C1[,1]^2+C1[,2]^2)
      Cos2<-C2[,1]/sqrt(C2[,1]^2+C2[,2]^2)
      Cosm1<-which.min(Cos1)
      Cosm2<-which.min(Cos2)
      sinm1<-C1[Cosm1,2]/sqrt(C1[Cosm1,1]^2+C1[Cosm1,2]^2)
      sinm2<-C2[Cosm2,2]/sqrt(C2[Cosm2,1]^2+C2[Cosm2,2]^2)
      sinab<-Cos1[Cosm1]*sinm2+Cos2[Cosm2]*sinm1
      if(sinab>=0){
        w1<-c(C1[Cosm1,2],-C1[Cosm1,1])
        w2<-c(-C2[Cosm2,2],C2[Cosm2,1])
      }
      else{
        w1<-c(-C1[Cosm1,2],C1[Cosm1,1])
        w2<-c(C2[Cosm2,2],-C2[Cosm2,1])
      }
      w<-cbind(w1[1]*lam[,1]+w2[1]*lam[,2],w1[2]*lam[,1]+w2[2]*lam[,2])
      wmax<-w[which(w[,2]<0),]
      wmin<-w[which(w[,2]>0),]
    }
    else if (C[1,2]>=0){
      Cos<-C[,1]/sqrt(C[,1]^2+C[,2]^2)
      Cosmin<-which.min(Cos)
      Cosmax<-which.max(Cos)
      w1<-c(C[Cosmin,2],-C[Cosmin,1])
      w2<-c(-C[Cosmax,2],C[Cosmax,1])
      w<-cbind(w1[1]*lam[,1]+w2[1]*lam[,2],w1[2]*lam[,1]+w2[2]*lam[,2])
      wmax<-w[which(w[,2]<0),]
      wmin<-w[which(w[,2]>0),]
    }
    else {
      Cos<-C[,1]/sqrt(C[,1]^2+C[,2]^2)
      Cosmin<-which.max(Cos)
      Cosmax<-which.min(Cos)
      w1<-c(C[Cosmin,2],-C[Cosmin,1])
      w2<-c(-C[Cosmax,2],C[Cosmax,1])
      w<-cbind(w1[1]*lam[,1]+w2[1]*lam[,2],w1[2]*lam[,1]+w2[2]*lam[,2])
      wmax<-w[which(w[,2]<0),]
      wmin<-w[which(w[,2]>0),]
    }
  }
  
  else {
    w1<-c(-C[2],C[1])
    w2<-c(C[2],-C[1])
    w<-cbind(w1[1]*lam[,1]+w2[1]*lam[,2],w1[2]*lam[,1]+w2[2]*lam[,2])
    wmax<-w[which(w[,2]<0),]
    wmin<-w[which(w[,2]>0),]
  }
  
  if (w1[2]*w2[2]<=0) xmin<-qnorm(q,mean=mu[1],sd=sqrt(Sigma[1,1]))
  else xmin<-(-5)
  # for w=(1,0), we have x>=q_{X_1}(q)
  x<-seq(xmin,5,by=0.1)
  
  ymax<-matrix(5,nrow=length(x),ncol = dim(wmax)[1])
  i<-dim(wmax)[1]
  while (i>0) {
    #qw<-qnorm(q,mean=0,sd=sqrt(wmax[i,1]^2+wmax[i,2]^2))
    qw<-qnorm(q,mean=sum(wmax[i,]*mu),sd=sqrt(wmax[i,1]^2*Sigma[1,1]+wmax[i,2]^2*Sigma[2,2]+2*wmax[i,1]*wmax[i,2]*Sigma[1,2]))
    ymax[,i]<-(qw-wmax[i,1]*x)/wmax[i,2]
    i<-i-1
  }
  ymax<-rowMins(ymax)
  
  ymin<-matrix(-5,nrow=length(x),ncol = dim(wmin)[1])
  i<-dim(wmin)[1]
  while (i>0) {
    #qw<-qnorm(q,mean=0,sd=sqrt(wmin[i,1]^2+wmin[i,2]^2))
    qw<-qnorm(q,mean=sum(wmin[i,]*mu),sd=sqrt(wmin[i,1]^2*Sigma[1,1]+wmin[i,2]^2*Sigma[2,2]+2*wmin[i,1]*wmin[i,2]*Sigma[1,2]))
    ymin[,i]<-(qw-wmin[i,1]*x)/wmin[i,2]
    i<-i-1
  }
  ymin<-rowMaxs(ymin)
  
  
  mdt<-data.frame(x,ymax,ymin)
  ggplot(mdt, aes(x = x)) + 
    xlab("z1")+ylab("z2")+
    geom_ribbon(aes(ymin=ymin,ymax =ymax), fill = 'gray60')+
    geom_point(data=P,aes(x=P[,1],y=P[,2]),shape=".")
}


