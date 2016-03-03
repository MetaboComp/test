library(doParallel)
library(mixOmics)
library(randomForest)
library(pROC)
library(MUVR)

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

DA=T
method='PLS'
nRep=4
nOuter=6
nInner=5
varRatio=0.75
fitness='AUROC'
ML=FALSE
methParam=list(compMax=3,mode='regression')
pred=FALSE
modReturn=FALSE


cl=makeCluster(detectCores())
registerDoParallel(cl)
# PLSModDA=MVWrap(X,Y2,ID,DA=T,nRep=nRep,method=method,fitness=fitness)
PLSModreg=MVWrap(X,Yr,ID,DA=F,nRep=nRep,method=method,fitness="RMSEP")
stopCluster(cl)

## <--------- Test zone for RF
xTrain=X[1:20,]
yTrain=Yr[1:20]
xVal=X[21:40,]
yVal=Yr[21:40]

rfmod=randomForest(xTrain,yTrain,xVal,yVal)

yClass=factor(sample(c('A','B'),40,replace=T))
