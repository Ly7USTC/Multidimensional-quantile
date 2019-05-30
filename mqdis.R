#library('Hmisc')
mqdis<-function(q,D,we=NA,C=diag(2)){
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
  ggplot(mdt, aes(x = x)) + 
    xlab("z1")+ylab("z2")+
    geom_ribbon(aes(ymin=ymin,ymax =ymax), fill = 'gray60')+
    geom_point(data=D,aes(x=D[,1],y=D[,2],size=Weights))
}


