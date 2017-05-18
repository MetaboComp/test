# Regression examples using "freelive" data
rm(list=ls())
library(doParallel)
library(MUVR)
data("freelive")
#  XRVIP=LCMS metabolomics of urine samples (selected metabolite features)
#  YR=Dietary consumption of whole grain rye in a free living population
#  IDR=Individual identifier (due to resampling after 2-3 months -> Dependent samples)

nCore=detectCores()-1
cl=makeCluster(nCore)
registerDoParallel(cl)
Regr_PLS_Quick=MUVR(X=XRVIP,Y=YR,ID=IDR,nRep=nCore,nOuter=5,varRatio=0.75,method='PLS') # Quick'N'Dirty
Regr_PLS_Full=MUVR(X=XRVIP,Y=YR,ID=IDR,nRep=5*nCore,nOuter=8,varRatio=0.9,method='PLS') # More proper model - Also more time consuming
Regr_RF_Quick=MUVR(X=XRVIP,Y=YR,ID=IDR,nRep=nCore,nOuter=5,varRatio=0.75,method='RF') # Quick'N'Dirty
Regr_RF_Full=MUVR(X=XRVIP,Y=YR,ID=IDR,nRep=5*nCore,nOuter=8,varRatio=0.9,method='RF') # More proper model - Also more time consuming
stopCluster(cl)


# Classification examples using "mosquito" data
rm(list=ls())
library(doParallel)
library(MUVR)
data("mosquito")

## RF classification
nCore=detectCores()-1
cl=makeCluster(nCore)
registerDoParallel(cl)
Class_PLS_Quick=MUVR(X=Xotu2,Y=Yotu,nRep=nCore,nOuter=5,varRatio=0.75,method='PLS',parallel=F) # Quick'N'Dirty
Class_PLS_Full=MUVR(X=Xotu2,Y=Yotu,nRep=5*nCore,nOuter=8,varRatio=0.9,method='PLS') # More proper model - Also more time consuming
Class_RF_Quick=MUVR(X=Xotu,Y=Yotu,nRep=nCore,nOuter=5,varRatio=0.75,method='RF',parallel=F) # Quick'N'Dirty
Class_RF_Full=MUVR(X=Xotu,Y=Yotu,nRep=5*nCore,nOuter=8,varRatio=0.9,method='RF') # More proper model - Also more time consuming
stopCluster(cl)

plotVAL(Class_RF)
plotMV(Class_RF)

## PLS classification
cl=makeCluster(3)
registerDoParallel(cl)
MA.pls=MVWrap(X=Xotu2[Yotu!='VK7',],Y=factor(Yotu[Yotu!='VK7']),nRep=12,method='PLS',varRatio=0.7,fitness='AUROC')
MA.test=testWrap(X=Xotu2[Yotu!='VK7',],Y=factor(Yotu[Yotu!='VK7']),nRep=12,method='PLS',varRatio=0.7,fitness='AUROC')
MM.pls=MVWrap(X=Xotu2[Yotu!='VK7',],Y=factor(Yotu[Yotu!='VK7']),nRep=12,method='PLS',varRatio=0.7,fitness='MISS')
MM.test=testWrap(X=Xotu2[Yotu!='VK7',],Y=factor(Yotu[Yotu!='VK7']),nRep=12,method='PLS',varRatio=0.7,fitness='MISS')
stopCluster(cl)
png(filename='Exclude/class.png')
par(mfrow=c(2,2))
par(mar=c(4,4,0,0)+.5)
plotVAL(MM.pls)
plotVAL(MA.pls)
plotVAL(MM.test)
plotVAL(MA.test)
dev.off()

## Multilevel ClinDiff
data("clinDiff")
cl=makeCluster(3)
registerDoParallel(cl)
ML.pls=MVWrap(X=clinDiff,ML=T,nRep=15,method='PLS',varRatio=0.9)
stopCluster(cl)
cl=makeCluster(3)
registerDoParallel(cl)
ML.test=testWrap(X=clinDiff,ML=T,nRep=15,method='PLS',varRatio=0.9)
stopCluster(cl)

plotVAL(ML.pls)
plotVAL(ML.test)


plotMV(ML.pls)
plotVAL(ML.pls)
class(ML.pls)


