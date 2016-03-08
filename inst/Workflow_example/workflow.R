library(doParallel)
library(MUVR)

load(file='freelivedata.rdata')

## Test VAL new Index function
cl=makeCluster(3)
registerDoParallel(cl)
Rtest.pls=testWrap(X=XRVIP,Y=YR,ID=IDR,nRep=30,method='PLS',varRatio=0.9)
stopCluster(cl)
plotVAL(Rtest.pls)
plot(YR,Rtest.pls$yPred[,2])
cor(YR,Rtest.pls$yPred[,2])

## PLS regression
cl=makeCluster(3)
registerDoParallel(cl)
R.pls=MVWrap(X=XRVIP,Y=YR,ID=IDR,nRep=30,method='PLS',varRatio=0.9)
stopCluster(cl)
plotMV(R.pls)
plotVAL(R.pls)

R.pls$nComp
R.pls$nVar
VIPs=R.pls$VIP[rank(R.pls$VIP[,1])<85,1]
head(sort(VIPs),10)

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
M.rf=testWrap(X=Xotu,Y=Yotu,nRep=30,method='RF',varRatio=0.9)
stopCluster(cl)
plotMV(M.rf)
plotVAL(M.rf)

## PLS classification
cl=makeCluster(3)
registerDoParallel(cl)
M.pls=testWrap(X=Xotu2,Y=Yotu,nRep=12,method='PLS',varRatio=0.7)
stopCluster(cl)

## Multilevel ClinDiff
load(file='clinDiff.rdata')
cl=makeCluster(4)
registerDoParallel(cl)
ML.pls=testWrap(X=clinDiff,ML=T,nRep=16,method='PLS',varRatio=0.9)
stopCluster(cl)
plotMV(ML.pls)
plotVAL(ML.pls)
class(ML.pls)
