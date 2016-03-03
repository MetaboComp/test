library(doParallel)
library(mixOmics)
library(randomForest)
library(pROC)
library(MUVR)

load(file='freelivedata.rdata')

## PLS regression
cl=makeCluster(4)
registerDoParallel(cl)
R.pls=MVWrap(X=XRVIP,Y=YR,ID=IDR,nRep=40,method='PLS',varRatio=0.9)
stopCluster(cl)

plot(YR,R.pls$yPred[,2])
cor(YR,R.pls$yPred[,2])
R.pls$nVar

## RF regression
cl=makeCluster(4)
registerDoParallel(cl)
R.rf=MVWrap(X=XRVIP,Y=YR,ID=IDR,nRep=40,method='RF',varRatio=0.9)
stopCluster(cl)

plot(YR,R.rf$yPred[,2])
cor(YR,R.rf$yPred[,2])
R.rf$nVar

plot(R.pls$yPred[,2],R.rf$yPred[,2])
cor(R.pls$yPred[,2],R.rf$yPred[,2])

