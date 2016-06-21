# Call in relevant libraries

library(doParallel)
library(MUVR)

# Freelive - Regression

rm(list=ls())
data("freelive")

## PLS regression

### Variable selections

### Univ
pUni=pUniFDR=numeric(ncol(XR))
for(v in 1:ncol(XR)) {
  Xv=XR[,v]
  rye.df=data.frame(X=Xv,Y=YR)
  pUni[v]=anova(lm(Y~X,data=rye.df))[1,5]
}
pUniFDR=p.adjust(pUni,method='fdr')
UV=colnames(XR)[pUniFDR<0.05]
### sPLS-1

library(mixOmics)
splsMod=mixOmics::spls(XR,YR,ncomp=3,mode='regression',keepX=c(1000,1000,1000))
splsVIP=vip(splsMod)[,3]
VIP=names(splsVIP)[splsVIP>1]

### sPLS-2

library(spls)
a=proc.time()[3]
sp=cv.spls(XR,YR,eta=seq(.3,.9,.1),K=1:5)
sp2=spls::spls(XR,YR,eta=sp$eta.opt,K=sp$K.opt)
varsp=colnames(sp2$x)[sp2$A]
a=(proc.time()[3]-a)/60

### MUVR: Full data

cl=makeCluster(3)
registerDoParallel(cl)
RMP_full=MVWrap(X=XR,Y=YR,ID=IDR,nRep=15,method='PLS',varRatio=0.9)
stopCluster(cl)

### MUVR: Uni-variate

cl=makeCluster(3)
registerDoParallel(cl)
RMP_uni=MVWrap(X=subset(XR,select = UV),Y=YR,ID=IDR,nRep=15,method='PLS',varRatio=0.9)
stopCluster(cl)

### MUVR: sPLS-filter 1

cl=makeCluster(3)
registerDoParallel(cl)
RMP_sp1=MVWrap(X=subset(XR,select = VIP),Y=YR,ID=IDR,nRep=15,method='PLS',varRatio=0.9)
stopCluster(cl)

### MUVR: sPLS filter 2

cl=makeCluster(3)
registerDoParallel(cl)
RMP_sp2=MVWrap(X=subset(XR,select = varsp),Y=YR,ID=IDR,nRep=15,method='PLS',varRatio=0.9)
stopCluster(cl)

### rdCV: Full data

cl=makeCluster(3)
registerDoParallel(cl)
RRP_full=rdCV(X=XR,Y=YR,ID=IDR,nRep=15,method='PLS')
stopCluster(cl)

### rdCV: Uni-variate

cl=makeCluster(3)
registerDoParallel(cl)
RRP_Uni=rdCV(X=subset(XR,select = UV),Y=YR,ID=IDR,nRep=15,method='PLS')
stopCluster(cl)

### rdCV: sPLS filter 1

cl=makeCluster(3)
registerDoParallel(cl)
RRP_sp1=rdCV(X=subset(XR,select = VIP),Y=YR,ID=IDR,nRep=15,method='PLS')
stopCluster(cl)

### rdCV: sPLS filter 2

cl=makeCluster(3)
registerDoParallel(cl)
RRP_sp2=rdCV(X=subset(XR,select = varsp),Y=YR,ID=IDR,nRep=15,method='PLS')
stopCluster(cl)

mods=ls(pattern='RRP')
nMod=length(mods)
R2=Q2=tVal=nC=nV=numeric(nMod)
for (m in 1:nMod) {
  eval(parse(text=paste('R2[m]=',mods[m],'$fitMetric$R2',sep='')))
  eval(parse(text=paste('Q2[m]=',mods[m],'$fitMetric$Q2',sep='')))
  eval(parse(text=paste('tVal[m]=',mods[m],'$calcMins',sep='')))
  eval(parse(text=paste('nC[m]=',mods[m],'$nComp',sep='')))
  eval(parse(text=paste('nV[m]=length(',mods[m],'$VIP)',sep='')))
}
RRPMat=data.frame(nComp=nC,nVar=nV,R2=R2,Q2=Q2,tVal=tVal,row.names=mods)

mods=ls(pattern='RMP')
nMod=length(mods)
R2=Q2=tVal=nC=nV=numeric(nMod)
for (m in 1:nMod) {
  eval(parse(text=paste('R2[m]=',mods[m],'$fitMetric$R2[2]',sep='')))
  eval(parse(text=paste('Q2[m]=',mods[m],'$fitMetric$Q2[2]',sep='')))
  eval(parse(text=paste('tVal[m]=',mods[m],'$calcMins',sep='')))
  eval(parse(text=paste('nC[m]=',mods[m],'$nComp[2]',sep='')))
  eval(parse(text=paste('nV[m]=',mods[m],'$nVar[2]',sep='')))
}
RMPMat=data.frame(nComp=nC,nVar=nV,R2=R2,Q2=Q2,tVal=tVal,row.names=mods)

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

### MUVR: sPLS filter 2

cl=makeCluster(3)
registerDoParallel(cl)
R3.rf=MVWrap(X=subset(XR,select = varsp),Y=YR,ID=IDR,nRep=15,method='RF',varRatio=0.9)
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

### rdCV: sPLS filter

cl=makeCluster(3)
registerDoParallel(cl)
Rrdcv3.rf=rdCV(X=subset(XR,select = varsp),Y=YR,ID=IDR,nRep=15,method='RF')
stopCluster(cl)

### Take out R2 and Q2

R.rf$fitMetric
R2.rf$fitMetric
R3.rf$fitMetric
R3.rf$nVar
Rrdcv.rf$fitMetric
Rrdcv2.rf$fitMetric
Rrdcv3.rf$fitMetric
Rrdcv3.rf$nVar

# save(R.pls,R2.pls,Rrdcv.pls,Rrdcv2.pls,R.rf,R2.rf,Rrdcv.rf,Rrdcv2.rf,file='Exclude/freeliveModels.rda')

### Permutation test sPLS2

nPerm=100
h0_sp2_pls=eta_sp2_pls=K_sp2_pls=nVar_sp2_pls=numeric(nPerm)
t_sp2_pls=proc.time()[3]
for (p in 1:nPerm) {
  cat ('permutation',p,'of',nPerm,'\n')
  yPerm=sample(YR)
  spp=cv.spls(XR,yPerm,eta=seq(.3,.9,.1),K=1:5,plot.it=F)
  spp2=spls::spls(XR,yPerm,eta=spp$eta.opt,K=spp$K.opt)
  eta_sp2_pls[p]=spp$eta.opt
  K_sp2_pls[p]=spp$K.opt
  nVar_sp2_pls[p]=length(spp2$A)
  varspp=colnames(spp2$x)[spp2$A]
  cl=makeCluster(3)
  registerDoParallel(cl)
  permMod=MVWrap(X=subset(XR,select=varspp),Y=yPerm,ID=IDR,nRep=12,method='PLS',varRatio=0.75)
  h0_sp2_pls[p]=permMod$fitMetric$Q2[2]
  stopCluster(cl)
}
t_sp2_pls=(proc.time()[3]-t_sp2_pls)/60


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
