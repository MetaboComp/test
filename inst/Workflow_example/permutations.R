library(doParallel)
library(MUVR)

###################################################
# Regression example using "freelive" data
# This example is for illustration only, since this X data was obtained after variable selection

rm(list=ls())
data("freelive")
#  XRVIP=LCMS metabolomics of urine samples (selected metabolite features)
#  YR=Dietary consumption of whole grain rye in a free living population
#  IDR=Individual identifier (due to resampling after 2-3 months -> Dependent samples)
nCore=detectCores()-1
cl=makeCluster(nCore)
registerDoParallel(cl)
# Declare modelling parameters
nRep=2*nCore   # Number of repetitions per actual model and permutations
nOuter=5       # Number of validation segments
varRatio=0.75  # Proportion of variables to keep per iteration during variable selection
method='PLS'   # Core modelling technique
size=1         # 1 for min, 2 for mid and 3 for max
nPerm=10       # Number of permutations (here set to 10 for illustration; normally set to 100 and inspected afterwards (see below))
permFit=numeric(nPerm)   # Allocate vector for permutation fitness
# Compute actual model and extract fitness metric
actual=MUVR(X=XRVIP,Y=YR,ID=IDR,nRep=nRep,nOuter=nOuter,varRatio=varRatio,method=method) # Quick'N'Dirty
actualFit=actual$fitMetric$Q2[size]
for (p in 1:nPerm) {
  cat('\nPermutation',p,'of',nPerm)
  YPerm=sample(YR)
  perm=MUVR(X=XRVIP,Y=YPerm,ID=IDR,nRep=nRep,nOuter=nOuter,varRatio=varRatio,method=method) # Quick'N'Dirty
  permFit[p]=perm$fitMetric$Q2[size]
}
stopCluster(cl)
pPerm(actual = actualFit, h0 = permFit)
plotPerm(actual = actualFit, h0 = permFit) # Look at histogram to assess whether Gaussian in shape 
pPerm(actual = actualFit, h0 = permFit, type = 'non') # If not Gaussian, make a non-parametric test instead
plotPerm(actual = actualFit, h0 = permFit, type = 'non') # If not Gaussian, make a non-parametric test instead

###################################################
# Classification examples using "mosquito" data
rm(list=ls())
data("mosquito")
#  Xotu: Microbiota OTU (16S rDNA) from mosquitos captured in 3 different villages in Burkina Faso
#  Yotu: One of three villages of capture
nCore=detectCores()-1
cl=makeCluster(nCore)
registerDoParallel(cl)
Class_RF_Quick=MUVR(X=Xotu,Y=Yotu,nRep=nCore,nOuter=5,varRatio=0.75,method='RF') # Quick'N'Dirty
# Class_RF_Full=MUVR(X=Xotu,Y=Yotu,nRep=5*nCore,nOuter=8,varRatio=0.9,method='RF') # More proper model - Also more time consuming
# Class_PLS_Quick=MUVR(X=Xotu,Y=Yotu,nRep=nCore,nOuter=5,varRatio=0.75,method='PLS') # Quick'N'Dirty
# Class_PLS_Full=MUVR(X=Xotu,Y=Yotu,nRep=5*nCore,nOuter=8,varRatio=0.9,method='PLS') # More proper model - Also more time consuming
stopCluster(cl)
plotVAL(Class_RF_Quick)
plotMV(Class_RF_Quick)
plotVIP(Class_RF_Quick,model='min')
plotStability(Class_RF_Quick,model='min')



###################################################
# Multilevel example using "crisp" data
rm(list=ls())
data("crisp")
#  crispEM: Effect matrix of 
nCore=detectCores()-1
cl=makeCluster(nCore)
registerDoParallel(cl)
ML_RF_Quick=MUVR(X=crispEM,ML=TRUE,nRep=nCore,nOuter=5,varRatio=0.75,method='RF') # Quick'N'Dirty
# ML_RF_Full=MUVR(X=crispEM,ML=TRUE,nRep=5*nCore,nOuter=8,varRatio=0.9,method='RF') # More proper model - Also more time consuming
# ML_PLS_Quick=MUVR(X=crispEM,ML=TRUE,nRep=nCore,nOuter=5,varRatio=0.75,method='PLS') # Quick'N'Dirty
# ML_PLS_Full=MUVR(X=crispEM,ML=TRUE,nRep=5*nCore,nOuter=8,varRatio=0.9,method='PLS') # More proper model - Also more time consuming
stopCluster(cl)
plotVAL(ML_RF_Quick)
plotMV(ML_RF_Quick)
plotVIP(ML_RF_Quick,model='min')
plotStability(ML_RF_Quick,model='min')


