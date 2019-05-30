mqnorm_2<-function(q,mu=c(0,0),Sigma=diag(2),C=diag(2)){
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
  g<-ggplot(mdt, aes(x = x)) + 
    xlab("z1")+ylab("z2")+
    xlim(c(-2,5))+ylim(c(-5,5))+
    geom_ribbon(aes(ymin=ymin,ymax =ymax), fill = 'grey60')+
    geom_point(data=P,aes(x=P[,1],y=P[,2]),shape=".")
  return(mdt)
}

mqnorm_3<-function(q,mu_1=c(0,0),Sigma_1=diag(2),mu_2=c(0,0),Sigma_2=diag(2),Lam,C=diag(2)){
  ##C is the points generating cone with default R+.
  P_1<-mvrnorm(10000,mu=mu_1,Sigma=Sigma_1)
  P_1<-as.data.frame(P_1)
  P_2<-mvrnorm(10000,mu=mu_2,Sigma=Sigma_2)
  P_2<-as.data.frame(P_2)
  
  mdt_1<-mqnorm_2(q,mu_1,Sigma_1,C)
  mdt_2<-mqnorm_2(q,mu_2,Sigma_2,C)
  
  mu_s<-Lam*mu_1+(1-Lam)*mu_2
  Sigma_s<-Lam^2*Sigma_1+(1-Lam)^2*Sigma_2
  mdt_s<-mqnorm_2(q,mu_s,Sigma_s,C)
  
  ggplot(mdt_1, aes(x = x)) + 
    xlab("z1")+ylab("z2")+
    xlim(c(-2,5))+ylim(c(-5,5))+
    geom_ribbon(aes(ymin=ymin,ymax =ymax), fill = 'red',alpha=0.5)+
    geom_point(data=P_1,aes(x=P_1[,1],y=P_1[,2]),shape=".",col='red')+
    geom_ribbon(data=mdt_s,aes(ymin=ymin,ymax =ymax), fill = 'black',alpha=0.9)+
    geom_ribbon(data=mdt_2,aes(x=x,ymin=ymin,ymax =ymax), fill = 'blue',alpha=0.5)+
    geom_point(data=P_2,aes(x=P_2[,1],y=P_2[,2]),shape=".",col='blue')
}
mqnorm_3(0.5,mu_1,Sig_1,mu_2,Sig_2,0.5)
mqnorm_3(0.05,mu_1,Sig_1,mu_2,Sig_2,0.5)
mqnorm_3(0.95,mu_1,Sig_1,mu_2,Sig_2,0.5)
#library('Hmisc')
mqdis_2<-function(q,D,we=NA,C=diag(2)){
  ## 2-dimensional q-quantile for discrete distribution on points D.
  ## w is weight on each point. sum(w)=1.
  ##C is the points generating cone with default R+.
  
  if(is.matrix(D)) {
    D<-as.data.frame(D)
    if(anyNA(we)) we=rep(1/dim(D)[1],dim(D)[1])  
  }
  else {
    D<-t(as.data.frame(D))
    we<-1
  }
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
  
  if (w1[2]*w2[2]<=0) xmin<-wtd.quantile(D[,1],probs=q,weights = we,type="i/n")
  else xmin<-(-5)
  # for w=(1,0), we have x>=q_{X_1}(q)
  x<-seq(xmin,5,by=0.1)
  
  ymax<-matrix(5,nrow=length(x),ncol = dim(wmax)[1])
  i<-dim(wmax)[1]
  while (i>0) {
    qw<-wtd.quantile(wmax[i,1]*D[,1]+wmax[i,2]*D[,2],probs=q,weights = we,type="i/n")
    ymax[,i]<-(qw-wmax[i,1]*x)/wmax[i,2]
    i<-i-1
  }
  ymax<-rowMins(ymax)
  
  ymin<-matrix(-5,nrow=length(x),ncol = dim(wmin)[1])
  i<-dim(wmin)[1]
  while (i>0) {
    qw<-wtd.quantile(wmin[i,1]*D[,1]+wmin[i,2]*D[,2],probs=q,weights = we,type="i/n")
    ymin[,i]<-(qw-wmin[i,1]*x)/wmin[i,2]
    i<-i-1
  }
  ymin<-rowMaxs(ymin)
  
  
  mdt<-data.frame(x,ymax,ymin)
  D<-cbind(D,we)
  colnames(D)<-c("","","Weights")
  g<-ggplot(mdt, aes(x = x)) + 
    xlab("z1")+ylab("z2")+
    geom_ribbon(aes(ymin=ymin,ymax =ymax), fill = 'gray60')+
    geom_point(data=D,aes(x=D[,1],y=D[,2],size=Weights))
  return(mdt)
}

mqemp_2<-function(q,F,C=diag(2)){
  ### 2-dimensional q-quantile for empirical distribution.
  ##C is the points generating cone with default R+.
  ##F is the points in empirical distribution.
  F<-as.data.frame(F)
  lam<-seq(0,1,by=0.01)
  lam<-cbind(1-lam,lam)
  if (is.matrix(C)){
    if (any(C[,2]>=0)&&any(C[,2]<0)){
      C1<-C[which(C[,2]>=0),]
      C1<-matrix(C1,ncol=2,byrow=FALSE)
      C2<-C[which(C[,2]<0),]
      C2<-matrix(C2,ncol=2,byrow=FALSE)
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
  
  if (w1[2]*w2[2]<=0) xmin<-quantile(F[,1],q,type=1)
  else xmin<-(-5)
  # for w=(1,0), we have x>=q_{X_1}(q)
  x<-seq(xmin,5,by=0.1)
  
  ymax<-matrix(5,nrow=length(x),ncol = dim(wmax)[1])
  i<-dim(wmax)[1]
  while (i>0) {
    #qw<-qnorm(q,mean=0,sd=sqrt(wmax[i,1]^2+wmax[i,2]^2))
    qw<-quantile(wmax[i,1]*F[,1]+wmax[i,2]*F[,2],q,type=1)
    ymax[,i]<-(qw-wmax[i,1]*x)/wmax[i,2]
    i<-i-1
  }
  ymax<-rowMins(ymax)
  
  ymin<-matrix(-5,nrow=length(x),ncol = dim(wmin)[1])
  i<-dim(wmin)[1]
  while (i>0) {
    #qw<-qnorm(q,mean=0,sd=sqrt(wmin[i,1]^2+wmin[i,2]^2))
    qw<-quantile(wmin[i,1]*F[,1]+wmin[i,2]*F[,2],q,type=1)
    ymin[,i]<-(qw-wmin[i,1]*x)/wmin[i,2]
    i<-i-1
  }
  ymin<-rowMaxs(ymin)
  
  
  mdt<-data.frame(x,ymax,ymin)
  g<-ggplot(mdt, aes(x = x)) + 
    xlab("z1")+ylab("z2")+
    geom_ribbon(aes(ymin=ymin,ymax =ymax), fill = 'gray60')+
    geom_point(data=F,aes(x=F[,1],y=F[,2]),shape=".")
  return(mdt)
}

CxLSnorm<-function(q,mu_1=c(0,0),Sigma_1=diag(2),mu_2=c(0,0),Sigma_2=diag(2),Lam,C=diag(2)){
  P_1<-mvrnorm(10000,mu=mu_1,Sigma=Sigma_1)
  P_2<-mvrnorm(10000,mu=mu_2,Sigma=Sigma_2)
  n_1<-dim(P_1)[1]
  n_2<-dim(P_2)[1]
  P_s<-rbind(P_1,P_2)
  P_1<-as.data.frame(P_1)
  P_2<-as.data.frame(P_2)
  mdt_1<-mqnorm_2(q,mu_1,Sigma_1,C)
  mdt_2<-mqnorm_2(q,mu_2,Sigma_2,C)
  we<-c(rep(Lam,n_1),rep(1-Lam,n_2))
  mdt_s<-mqdis_2(q,P_s,we,C)
  g<-ggplot(mdt_1, aes(x = x)) + 
    xlab("z1")+ylab("z2")+
    xlim(c(-2,5))+ylim(c(-5,5))+
    geom_ribbon(aes(ymin=ymin,ymax =ymax), fill = 'red',alpha=0.5)+
    geom_point(data=P_1,aes(x=P_1[,1],y=P_1[,2]),shape=".",col='red')+
    geom_ribbon(data=mdt_s,aes(ymin=ymin,ymax =ymax), fill = 'black',alpha=0.9)+
    geom_ribbon(data=mdt_2,aes(x=x,ymin=ymin,ymax =ymax), fill = 'blue',alpha=0.5)+
    geom_point(data=P_2,aes(x=P_2[,1],y=P_2[,2]),shape=".",col='blue')
  return(g)
}

CxLSemp<-function(q,F_1,F_2,Lam,C=diag(2)){
  F_s<-rbind(F_1,F_2)
  n_1<-dim(F_1)[1]
  n_2<-dim(F_2)[1]
  we<-c(rep(Lam/n_1,n_1),rep((1-Lam)/n_2,n_2))
  mdt_1<-mqemp_2(q,F_1,C)
  mdt_2<-mqemp_2(q,F_2,C)
  mdt_s<-mqdis_2(q,F_s,we,C)
  F_1<-as.data.frame(F_1)
  F_2<-as.data.frame(F_2)
  g<-ggplot(mdt_1, aes(x = x)) + 
    xlab("z1")+ylab("z2")+
    xlim(c(-2,5))+ylim(c(-5,5))+
    geom_ribbon(aes(ymin=ymin,ymax =ymax), fill = 'red',alpha=0.5)+
    geom_point(data=F_1,aes(x=F_1[,1],y=F_1[,2],size=Lam,alpha=0.5),col='red')+
    geom_ribbon(data=mdt_s,aes(ymin=ymin,ymax =ymax), fill = 'black',alpha=0.9)+
    geom_ribbon(data=mdt_2,aes(x=x,ymin=ymin,ymax =ymax), fill = 'blue',alpha=0.5)+
    geom_point(data=F_2,aes(x=F_2[,1],y=F_2[,2],size=1-Lam,alpha=0.5),col='blue')+
    scale_size_identity()
  return(g)
}

#CxLSemp(0.95,P_1,P_2,0.9)
#CxLSnorm(0.95,mu_1,Sig_1,mu_2,Sig_2,0.9)
#C<-matrix(c(0,1,1,-1),byrow=T,nrow=2)
#CxLSnorm(0.95,mu_1,Sig_1,mu_2,Sig_2,0.9,C)
#CxLSemp(0.95,P_1,P_2,0.9,C)
