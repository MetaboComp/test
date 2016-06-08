# Call in relevant libraries

library(doParallel)
library(MUVR)

# Freelive - Regression

rm(list=ls())
data("freelive")

## PLS regression

### MUVR: Full data

cl=makeCluster(3)
registerDoParallel(cl)
R.pls=MVWrap(X=XR,Y=YR,ID=IDR,nRep=15,method='PLS',varRatio=0.9)
stopCluster(cl)

### MUVR: sPLS-filter

cl=makeCluster(3)
registerDoParallel(cl)
R2.pls=MVWrap(X=XRVIP,Y=YR,ID=IDR,nRep=15,method='PLS',varRatio=0.9)
stopCluster(cl)

### rdCV: Full data

cl=makeCluster(3)
registerDoParallel(cl)
Rrdcv.pls=rdCV(X=XR,Y=YR,ID=IDR,nRep=15,method='PLS')
stopCluster(cl)

### rdCV: sPLS filter

cl=makeCluster(3)
registerDoParallel(cl)
Rrdcv2.pls=rdCV(X=XRVIP,Y=YR,ID=IDR,nRep=15,method='PLS')
stopCluster(cl)

### Take out R2 and Q2

R.pls$fitMetric
R2.pls$fitMetric
Rrdcv.pls$fitMetric
Rrdcv2.pls$fitMetric

## RF regression

### MUVR: Full data

cl=makeCluster(3)
registerDoParallel(cl)
R.rf=MVWrap(X=XR,Y=YR,ID=IDR,nRep=15,method='RF',varRatio=0.9)
stopCluster(cl)

### MUVR: sPLS filter

cl=makeCluster(3)
registerDoParallel(cl)
R2.rf=MVWrap(X=XRVIP,Y=YR,ID=IDR,nRep=15,method='RF',varRatio=0.9)
stopCluster(cl)

### rdCV: Full data

cl=makeCluster(3)
registerDoParallel(cl)
Rrdcv.rf=rdCV(X=XR,Y=YR,ID=IDR,nRep=15,method='RF')
stopCluster(cl)

### rdCV: sPLS filter

cl=makeCluster(3)
registerDoParallel(cl)
Rrdcv2.rf=rdCV(X=XRVIP,Y=YR,ID=IDR,nRep=15,method='RF')
stopCluster(cl)

### Take out R2 and Q2

R.rf$fitMetric
R2.rf$fitMetric
Rrdcv.rf$fitMetric
Rrdcv2.rf$fitMetric

# save(R.pls,R2.pls,Rrdcv.pls,Rrdcv2.pls,R.rf,R2.rf,Rrdcv.rf,Rrdcv2.rf,file='Exclude/freeliveModels.rda')

# Mosquito - Classification

rm(list=ls())
data("mosquito")

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
