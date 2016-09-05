# Call in relevant libraries

library(doParallel)
library(MUVR)

# Freelive - Regression

rm(list=ls())
data("freelive")

## PLS regression
#  Different types of variable selections are performed prior to main analysis to 
#  examine impact of performance
#  Data for regression is
#  X=LCMS metabolomics of urine samples
#  Y=Dietary consumption of whole grain rye in a free living population
#  Individuals are resampled after 2-3 months
#  ID variable is added to perform subsampling only on iondependent samples


### Variable selections

### Univ

tUni=proc.time()[3]
pUni=pUniFDR=numeric(ncol(XR))
for(v in 1:ncol(XR)) {
  Xv=XR[,v]
  rye.df=data.frame(X=Xv,Y=YR)
  pUni[v]=anova(lm(Y~X,data=rye.df))[1,5]
}
pUniFDR=p.adjust(pUni,method='fdr')
UV=colnames(XR)[pUniFDR<0.05]
tUni=(proc.time()[3]-tUni)/60

### sPLS-1

library(mixOmics)
tsp1=proc.time()[3]
splsMod=mixOmics::spls(XR,YR,ncomp=3,mode='regression',keepX=c(500,500,500))
splsVIP=vip(splsMod)[,3]
VIP=names(splsVIP)[splsVIP>1]
tsp1=(proc.time()[3]-tsp1)/60

### sPLS-2

library(spls)
tsp2=proc.time()[3]
sp=cv.spls(XR,YR,eta=seq(.3,.9,.1),K=1:5)
sp2=spls::spls(XR,YR,eta=sp$eta.opt,K=sp$K.opt)
varsp=colnames(sp2$x)[sp2$A]
tsp2=(proc.time()[3]-tsp2)/60

### VSURF

library(VSURF)
vs=VSURF(XR,YR,parallel=T)
varVSI=colnames(XR)[vs$varselect.interp]
varVSP=colnames(XR)[vs$varselect.pred]
tVSURF=as.numeric(vs$overall.time*60)

### Boruta

library(Boruta)
bor=Boruta(XR,YR,holdHistory=F)
varBor=colnames(XR)[bor$finalDecision=='Confirmed']
tBor=bor$timeTaken

tVS=c(tBor,0,tsp1,tsp2,tUni,tVSURF,tVSURF)
names(tVS)=c('Boruta','Full','sPLS-1','sPLS-2','Univ','VSURF-I','VSURF-P')

### MUVR: Full data

cl=makeCluster(3)
registerDoParallel(cl)
RMP_full=MUVR(X=XR,Y=YR,ID=IDR,nRep=15,method='PLS',varRatio=0.9)
stopCluster(cl)

### MUVR: Uni-variate

cl=makeCluster(3)
registerDoParallel(cl)
RMP_uni=MUVR(X=subset(XR,select = UV),Y=YR,ID=IDR,nRep=15,method='PLS',varRatio=0.9)
stopCluster(cl)

### MUVR: sPLS-filter 1

cl=makeCluster(3)
registerDoParallel(cl)
RMP_sp1=MUVR(X=subset(XR,select = VIP),Y=YR,ID=IDR,nRep=15,method='PLS',varRatio=0.9)
stopCluster(cl)

### MUVR: sPLS filter 2

cl=makeCluster(3)
registerDoParallel(cl)
RMP_sp2=MUVR(X=subset(XR,select = varsp),Y=YR,ID=IDR,nRep=15,method='PLS',varRatio=0.9)
stopCluster(cl)

### MUVR: VSURF-I

cl=makeCluster(3)
registerDoParallel(cl)
RMP_VSI=MUVR(X=subset(XR,select = varVSI),Y=YR,ID=IDR,nRep=15,method='PLS',varRatio=0.9)
stopCluster(cl)

### MUVR: VSURF-P

cl=makeCluster(3)
registerDoParallel(cl)
RMP_VSP=MUVR(X=subset(XR,select = varVSP),Y=YR,ID=IDR,nRep=15,method='PLS',varRatio=0.9)
stopCluster(cl)

### MUVR: Boruta

cl=makeCluster(3)
registerDoParallel(cl)
RMP_Bor=MUVR(X=subset(XR,select = varBor),Y=YR,ID=IDR,nRep=15,method='PLS',varRatio=0.9)
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

### rdCV: VSURF-I

cl=makeCluster(3)
registerDoParallel(cl)
RRP_VSI=rdCV(X=subset(XR,select = varVSI),Y=YR,ID=IDR,nRep=15,method='PLS')
stopCluster(cl)

### rdCV: VSURF-P

cl=makeCluster(3)
registerDoParallel(cl)
RRP_VSP=rdCV(X=subset(XR,select = varVSP),Y=YR,ID=IDR,nRep=15,method='PLS')
stopCluster(cl)

### rdCV: Boruta

cl=makeCluster(3)
registerDoParallel(cl)
RRP_Bor=rdCV(X=subset(XR,select = varBor),Y=YR,ID=IDR,nRep=15,method='PLS')
stopCluster(cl)

mods=ls(pattern='RRP_')
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

mods=ls(pattern='RMP_')
nMod=length(mods)
R2=Q2=tVal=nC=nV=numeric(nMod)
for (m in 1:nMod) {
  eval(parse(text=paste('R2[m]=',mods[m],'$fitMetric$R2[2]',sep='')))
  eval(parse(text=paste('Q2[m]=',mods[m],'$fitMetric$Q2[2]',sep='')))
  eval(parse(text=paste('tVal[m]=',mods[m],'$calcMins',sep='')))
  eval(parse(text=paste('nC[m]=',mods[m],'$nComp[2]',sep='')))
  eval(parse(text=paste('nV[m]=',mods[m],'$nVar[2]',sep='')))
}
RMPMat=data.frame(nComp=round(nC),nVar=round(nV),R2=round(R2,2),Q2=round(Q2,2),tVal=round(tVal,2),row.names=mods)
write.csv2(RMPMat,file='RMPMat.csv')

## RF regression

### MUVR: Full data

cl=makeCluster(3)
registerDoParallel(cl)
RMR_full=MUVR(X=XR,Y=YR,ID=IDR,nRep=15,method='RF',varRatio=0.9)
stopCluster(cl)

### MUVR: Uni-variate

cl=makeCluster(3)
registerDoParallel(cl)
RMR_uni=MUVR(X=subset(XR,select = UV),Y=YR,ID=IDR,nRep=15,method='RF',varRatio=0.9)
stopCluster(cl)

### MUVR: sPLS filter

cl=makeCluster(3)
registerDoParallel(cl)
RMR_sp1=MUVR(X=XRVIP,Y=YR,ID=IDR,nRep=15,method='RF',varRatio=0.9)
stopCluster(cl)

### MUVR: sPLS filter 2

cl=makeCluster(3)
registerDoParallel(cl)
RMR_sp2=MUVR(X=subset(XR,select = varsp),Y=YR,ID=IDR,nRep=15,method='RF',varRatio=0.9)
stopCluster(cl)

### MUVR: VSURF-I

cl=makeCluster(3)
registerDoParallel(cl)
RMR_VSI=MUVR(X=subset(XR,select = varVSI),Y=YR,ID=IDR,nRep=15,method='RF',varRatio=0.9)
stopCluster(cl)

### MUVR: VSURF-P

cl=makeCluster(3)
registerDoParallel(cl)
RMR_VSP=MUVR(X=subset(XR,select = varVSP),Y=YR,ID=IDR,nRep=15,method='RF',varRatio=0.9)
stopCluster(cl)

### MUVR: Boruta

cl=makeCluster(3)
registerDoParallel(cl)
RMR_Bor=MUVR(X=subset(XR,select = varBor),Y=YR,ID=IDR,nRep=15,method='RF',varRatio=0.9)
stopCluster(cl)

### rdCV: Full data

cl=makeCluster(3)
registerDoParallel(cl)
RRR_full=rdCV(X=XR,Y=YR,ID=IDR,nRep=15,method='RF',varRatio=0.9)
stopCluster(cl)

### rdCV: Uni-variate

cl=makeCluster(3)
registerDoParallel(cl)
RRR_uni=rdCV(X=subset(XR,select = UV),Y=YR,ID=IDR,nRep=15,method='RF',varRatio=0.9)
stopCluster(cl)

### rdCV: sPLS filter

cl=makeCluster(3)
registerDoParallel(cl)
RRR_sp1=rdCV(X=XRVIP,Y=YR,ID=IDR,nRep=15,method='RF',varRatio=0.9)
stopCluster(cl)

### rdCV: sPLS filter 2

cl=makeCluster(3)
registerDoParallel(cl)
RRR_sp2=rdCV(X=subset(XR,select = varsp),Y=YR,ID=IDR,nRep=15,method='RF',varRatio=0.9)
stopCluster(cl)

### rdCV: VSURF-I

cl=makeCluster(3)
registerDoParallel(cl)
RRR_VSI=rdCV(X=subset(XR,select = varVSI),Y=YR,ID=IDR,nRep=15,method='RF',varRatio=0.9)
stopCluster(cl)

### rdCV: VSURF-P

cl=makeCluster(3)
registerDoParallel(cl)
RRR_VSP=rdCV(X=subset(XR,select = varVSP),Y=YR,ID=IDR,nRep=15,method='RF',varRatio=0.9)
stopCluster(cl)

### rdCV: Boruta

cl=makeCluster(3)
registerDoParallel(cl)
RRR_Bor=rdCV(X=subset(XR,select = varBor),Y=YR,ID=IDR,nRep=15,method='RF',varRatio=0.9)
stopCluster(cl)



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


