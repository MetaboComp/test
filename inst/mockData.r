# Mock data
rm(list=ls())

nvar=17
nsamp=40
X=matrix(runif(nvar*nsamp),nrow=nsamp)
colnames(X)=paste('var',1:nvar,sep='')
rownames(X)=paste('samp',1:nsamp,sep='')
Yr=runif(nsamp)
Y2=sample(c(-1,1),nsamp,replace=TRUE)
Y3=sample(c('A','B','C'),nsamp,replace=TRUE)
ID=sample(1:20,nsamp,replace=TRUE)
Y=Y2

DA=T
method='PLS'
nRep=5
nOuter=6
nInner=5
varRatio=0.75
fitness='AUROC'
ML=FALSE
pred=FALSE
XP=NULL
methParam=list(compMax=3,mode='regression')
