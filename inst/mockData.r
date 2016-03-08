library(doParallel)
library(MUVR)

# Mock data
rm(list=ls())

nvar=17
nsamp=40
X=matrix(runif(nvar*nsamp,0.9,1.1),nrow=nsamp,dimnames=list(paste('samp',1:nsamp,sep=''),paste('var',1:nvar,sep='')))
Yr=runif(nsamp)
Y2=factor(sample(c('A','B'),nsamp,replace=TRUE))
Y3=factor(sample(c('A','B','C'),nsamp,replace=TRUE))
# X2 Manipulation
X2=matrix(runif(nvar*nsamp,0.9,1.1),nrow=nsamp,dimnames=list(paste('samp',1:nsamp,sep=''),paste('var',1:nvar,sep='')))
X2[Y2=='A',1:2]=1.1*X2[Y2=='A',1:2]
X2[Y2=='A',3:4]=0.9*X2[Y2=='A',3:4]
MT=sample(c(TRUE,FALSE),nsamp,replace=TRUE)
X2[MT,5:6]=0.5*X2[MT,5:6]
X2[MT & Y2=='A',7:8]=1.3*X2[MT & Y2=='A',7:8]
# X3 manipulation
X3=matrix(runif(nvar*nsamp,0.9,1.1),nrow=nsamp,dimnames=list(paste('samp',1:nsamp,sep=''),paste('var',1:nvar,sep='')))
X3[Y3=='A',1:2]=1.1*X3[Y3=='A',1:2]
X3[Y3=='C',1:2]=0.9*X3[Y3=='C',1:2]
X3[Y3=='A',3:4]=0.9*X3[Y3=='A',3:4]
X3[Y3=='B',3:4]=1.1*X3[Y3=='B',3:4]
MT=sample(c(TRUE,FALSE),nsamp,replace=TRUE)
X3[MT,5:6]=0.5*X3[MT,5:6]
X3[MT & Y3=='B',7:8]=1.3*X3[MT & Y3=='B',7:8]
X3[MT & Y3=='C',7:8]=0.7*X3[MT & Y3=='C',7:8]
# ID=sample(1:20,nsamp,replace=TRUE)

cl=makeCluster(4)
registerDoParallel(cl)
test1=testWrap(X2,Y2,nRep=5,method='PLS',modReturn=F)
test2=testWrap(X2,Y2,nRep=5,method='RF',modReturn=F)
test3=testWrap(X3,Y3,nRep=5,method='PLS',modReturn=F)
test4=testWrap(X3,Y3,nRep=5,method='RF',modReturn=F)
stopCluster(cl)

plotMV(test1)
plotMV(test2)
plotMV(test3)
plotMV(test4)

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
