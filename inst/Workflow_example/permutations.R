# Call in relevant libraries

library(doParallel)
library(MUVR)

# Freelive - Regression - Permutation

rm(list=ls())
data("freelive")
nPerm=100

### Univ
UP_corYR=UP_R2=UP_Q2=numeric(nPerm)
for (p in 48:nPerm) {
  YPerm=sample(YR)
  UP_corYR[p]=cor(YR,YPerm)
  pUni=pUniFDR=numeric(ncol(XR))
  for(v in 1:ncol(XR)) {
    Xv=XR[,v]
    rye.df=data.frame(X=Xv,Y=YPerm)
    pUni[v]=anova(lm(Y~X,data=rye.df))[1,5]
  }
  pUniFDR=p.adjust(pUni,method='fdr')
  UV=colnames(XR)[pUniFDR<0.05]
  if (length(UV)<100) UV=colnames(XR)[pUniFDR<sort(pUniFDR)[101]]
  cl=makeCluster(3)
  registerDoParallel(cl)
  permMod=MUVR(X=subset(XR,select=UV),Y=YPerm,ID=IDR,nRep=15,method='PLS',varRatio=0.75)
  stopCluster(cl)
  UP_R2[p]=permMod$fitMetric$R2[1]
  UP_Q2[p]=permMod$fitMetric$Q2[1]
}

save(UP_Q2,UP_R2,UP_corYR,file='UP.rda')

### sPLS-1
sP1P_corYR=sP1P_R2=sP1P_Q2=numeric(nPerm)
for (p in 1:nPerm) {
  YPerm=sample(YR)
  sP1P_corYR[p]=cor(YR,YPerm)
  splsMod=mixOmics::spls(XR,YR,ncomp=3,mode='regression',keepX=c(500,500,500))
  splsVIP=mixOmics::vip(splsMod)[,3]
  sP1=names(splsVIP)[splsVIP>1]
  cl=makeCluster(3)
  registerDoParallel(cl)
  permMod=MUVR(X=subset(XR,select=sP1),Y=YPerm,ID=IDR,nRep=15,method='PLS',varRatio=0.75)
  stopCluster(cl)
  sP1P_R2[p]=permMod$fitMetric$R2[1]
  sP1P_Q2[p]=permMod$fitMetric$Q2[1]
}

### sPLS-2
sP2P_corYR=sP2P_R2=sP2P_Q2=numeric(nPerm)
for (p in 1:nPerm) {
  YPerm=sample(YR)
  sP2P_corYR[p]=cor(YR,YPerm)
  sp=cv.spls(XR,YR,eta=seq(.3,.9,.1),K=1:5)
  sp2=spls::spls(XR,YPerm,eta=sp$eta.opt,K=sp$K.opt)
  sP2=colnames(sp2$x)[sp2$A]
  cl=makeCluster(3)
  registerDoParallel(cl)
  permMod=MUVR(X=subset(XR,select=sP2),Y=YPerm,ID=IDR,nRep=15,method='PLS',varRatio=0.75)
  stopCluster(cl)
  sP2P_R2[p]=permMod$fitMetric$R2[1]
  sP2P_Q2[p]=permMod$fitMetric$Q2[1]
}

### Boruta
BP_corYR=BP_R2=BP_Q2=numeric(nPerm)
for (p in 1:nPerm) {
  YPerm=sample(YR)
  BP_corYR[p]=cor(YR,YPerm)
  bor=Boruta::Boruta(XR,YPerm,holdHistory=F)
  Bor=colnames(XR)[bor$finalDecision=='Confirmed']
  cl=makeCluster(3)
  registerDoParallel(cl)
  permMod=MUVR(X=subset(XR,select=Bor),Y=YPerm,ID=IDR,nRep=15,method='PLS',varRatio=0.75)
  stopCluster(cl)
  BP_R2[p]=permMod$fitMetric$R2[1]
  BP_Q2[p]=permMod$fitMetric$Q2[1]
}

### Full data
FP_corYR=FP_R2=FP_Q2=numeric(nPerm)
for (p in 1:nPerm) {
  YPerm=sample(YR)
  FP_corYR[p]=cor(YR,YPerm)
  cl=makeCluster(3)
  registerDoParallel(cl)
  permMod=MUVR(X=XR,Y=YPerm,ID=IDR,nRep=15,method='PLS',varRatio=0.75)
  stopCluster(cl)
  FP_R2[p]=permMod$fitMetric$R2[1]
  FP_Q2[p]=permMod$fitMetric$Q2[1]
}

### rdCV
rP_corYR=rP_R2=rP_Q2=numeric(nPerm)
for (p in 1:nPerm) {
  YPerm=sample(YR)
  rP_corYR[p]=cor(YR,YPerm)
  cl=makeCluster(3)
  registerDoParallel(cl)
  permMod=rdCV(X=XR,Y=YPerm,ID=IDR,nRep=15,method='PLS')
  stopCluster(cl)
  rP_R2[p]=permMod$fitMetric$R2[1]
  rP_Q2[p]=permMod$fitMetric$Q2[1]
}

### Full RF
FR_corYR=FR_R2=FR_Q2=numeric(nPerm)
for (p in 1:nPerm) {
  YPerm=sample(YR)
  FR_corYR[p]=cor(YR,YPerm)
  cl=makeCluster(3)
  registerDoParallel(cl)
  permMod=MUVR(X=XR,Y=YPerm,ID=IDR,nRep=15,method='RF',varRatio=0.75)
  stopCluster(cl)
  FR_R2[p]=permMod$fitMetric$R2[1]
  FR_Q2[p]=permMod$fitMetric$Q2[1]
}

