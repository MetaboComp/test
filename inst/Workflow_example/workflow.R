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
# Regr_PLS_Full=MUVR(X=XRVIP,Y=YR,ID=IDR,nRep=5*nCore,nOuter=8,varRatio=0.9,method='PLS') # More proper model - Also more time consuming
# Regr_RF_Quick=MUVR(X=XRVIP,Y=YR,ID=IDR,nRep=nCore,nOuter=5,varRatio=0.75,method='RF') # Quick'N'Dirty
# Regr_RF_Full=MUVR(X=XRVIP,Y=YR,ID=IDR,nRep=5*nCore,nOuter=8,varRatio=0.9,method='RF') # More proper model - Also more time consuming
stopCluster(cl)
plotVAL(Regr_PLS_Quick)
plotMV(Regr_PLS_Quick,model='mid')
plotVIP(Regr_PLS_Quick,model='mid')

# Classification examples using "mosquito" data
rm(list=ls())
library(doParallel)
library(MUVR)
data("mosquito")
#  Xotu=Microbiota OTU (16S rDNA) from mosquitos captured in 3 different villages in Burkina Faso
#  Yotu=Villages of capture

## RF classification
nCore=detectCores()-1
cl=makeCluster(nCore)
registerDoParallel(cl)
Class_RF_Quick=MUVR(X=Xotu,Y=Yotu,nRep=nCore,nOuter=5,varRatio=0.75,method='RF') # Quick'N'Dirty
# Class_RF_Full=MUVR(X=Xotu,Y=Yotu,nRep=5*nCore,nOuter=8,varRatio=0.9,method='RF') # More proper model - Also more time consuming
# Class_PLS_Quick=MUVR(X=Xotu,Y=Yotu,nRep=nCore,nOuter=5,varRatio=0.75,method='PLS') # Quick'N'Dirty
# Class_PLS_Full=MUVR(X=Xotu,Y=Yotu,nRep=5*nCore,nOuter=8,varRatio=0.9,method='PLS') # More proper model - Also more time consuming
stopCluster(cl)
plotVAL(Class_PLS_Quick)
plotMV(Class_PLS_Quick)


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


