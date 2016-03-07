library(doParallel)
library(MUVR)

load(file='freelivedata.rdata')

## Test VAL new Index function
cl=makeCluster(4)
registerDoParallel(cl)
RVAL.pls=testWrap(X=XRVIP,Y=YR,ID=IDR,nRep=40,method='PLS',varRatio=0.7)
stopCluster(cl)
plotVAL(RVAL.pls)
plot(YR,RVAL.pls$yPred[,2])
cor(YR,RVAL.pls$yPred[,2])

## PLS regression
cl=makeCluster(3)
registerDoParallel(cl)
R.pls=MVWrap(X=XRVIP,Y=YR,ID=IDR,nRep=30,method='PLS',varRatio=0.9)
stopCluster(cl)

plotVAL(R.pls)

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



rm(list=ls())
load(file='otudata.rdata')

## RF classification
cl=makeCluster(3)
registerDoParallel(cl)
M.rf=testWrap(X=Xotu,Y=Yotu,nRep=12,method='RF',varRatio=0.7)
stopCluster(cl)

## PLS classification
cl=makeCluster(3)
registerDoParallel(cl)
M.pls=testWrap(X=Xotu2,Y=Yotu,nRep=12,method='PLS',varRatio=0.7)
stopCluster(cl)
