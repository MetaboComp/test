library(doParallel)
library(mixOmics)
library(randomForest)
library(pROC)
library(MUVR)

load(file='freelivedata.rdata')

## PLS regression
cl=makeCluster(3)
registerDoParallel(cl)
R.pls=MVWrap(X=XRVIP,Y=YR,ID=IDR,nRep=6,method='PLS')
stopCluster(cl)

plot(YR,R.pls$yPred[,2])
cor(YR,R.pls$yPred[,2])
R.pls$nVar

## RF regression
cl=makeCluster(3)
registerDoParallel(cl)
R.rf=MVWrap(X=XRVIP,Y=YR,ID=IDR,nRep=6,method='RF')
stopCluster(cl)

plot(YR,R.rf$yPred[,2])
cor(YR,R.rf$yPred[,2])
R.rf$nVar

plot(R.pls$yPred[,2],R.rf$yPred[,2])
cor(R.pls$yPred[,2],R.rf$yPred[,2])

