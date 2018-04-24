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
# Declare modelling parameters
nRep=2*nCore   # Number of repetitions per actual model and permutations
nOuter=5       # Number of validation segments
varRatio=0.75  # Proportion of variables to keep per iteration during variable selection
method='PLS'   # Core modelling technique
size=1         # 1 for min, 2 for mid and 3 for max
nPerm=10       # Number of permutations (here set to 10 for illustration; normally set to 100 and inspected afterwards (see below))
permFit=numeric(nPerm)   # Allocate vector for permutation fitness
# Compute actual model and extract fitness metric
nCore=detectCores()-1
cl=makeCluster(nCore)
registerDoParallel(cl)
actual=MUVR(X=XRVIP,Y=YR,ID=IDR,nRep=nRep,nOuter=nOuter,varRatio=varRatio,method=method) 
actualFit=actual$fitMetric$Q2[size]
for (p in 1:nPerm) {
  cat('\nPermutation',p,'of',nPerm)
  YPerm=sample(YR)
  perm=MUVR(X=XRVIP,Y=YPerm,ID=IDR,nRep=nRep,nOuter=nOuter,varRatio=varRatio,method=method)
  permFit[p]=perm$fitMetric$Q2[size]
}
stopCluster(cl)
pPerm(actual = actualFit, h0 = permFit)
plotPerm(actual = actualFit, h0 = permFit) # Look at histogram to assess whether Gaussian in shape 
pPerm(actual = actualFit, h0 = permFit, type = 'non') # If not Gaussian, make a non-parametric test instead
plotPerm(actual = actualFit, h0 = permFit, type = 'non') # And plot with non-parametric p-value instead

###################################################
# Classification examples using "mosquito" data
rm(list=ls())
data("mosquito")
#  Xotu: Microbiota OTU (16S rDNA) from mosquitos captured in 3 different villages in Burkina Faso
#  Yotu: One of three villages of capture
# Declare modelling parameters
nRep=2*nCore   # Number of repetitions per actual model and permutations
nOuter=5       # Number of validation segments
varRatio=0.75  # Proportion of variables to keep per iteration during variable selection
method='RF'    # Core modelling technique
model=1        # 1 for min, 2 for mid and 3 for max
nPerm=25       # Number of permutations (here set to 10 for illustration; normally set to 100 and inspected afterwards (see below))
permFit=numeric(nPerm)   # Allocate vector for permutation fitness
# Compute actual model and extract fitness metric
nCore=detectCores()-1
cl=makeCluster(nCore)
registerDoParallel(cl)
actual=MUVR(X=Xotu,Y=Yotu,nRep=nRep,nOuter=nOuter,varRatio=varRatio,method=method) 
actualFit=actual$miss[model]
for (p in 1:nPerm) {
  cat('\nPermutation',p,'of',nPerm)
  YPerm=sample(Yotu)
  perm=MUVR(X=Xotu,Y=YPerm,nRep=nRep,nOuter=nOuter,varRatio=varRatio,method=method)
  permFit[p]=perm$miss[model]
}
stopCluster(cl)
par(mfrow=c(1,2),mar=c(4,4,2,0)+.5)
pPerm(actual = actualFit, h0 = permFit)
plotPerm(actual = actualFit, h0 = permFit, xlab = 'Misclassifications') # Look at histogram to assess whether Gaussian in shape 
title(main = 'Parametric')
pPerm(actual = actualFit, h0 = permFit, type = 'non') # If not Gaussian, make a non-parametric test instead
plotPerm(actual = actualFit, h0 = permFit, type = 'non', xlab = 'Misclassifications') # And plot with non-parametric p-value instead
title(main = 'Non-parametric')


###################################################
# Multilevel example using "crisp" data
rm(list=ls())
data("crisp")
#  crispEM: Effect matrix of difference (Treatment1 - Treatment2) in crossover intervention
# Declare modelling parameters
nRep=2*nCore   # Number of repetitions per actual model and permutations
nOuter=5       # Number of validation segments
varRatio=0.75  # Proportion of variables to keep per iteration during variable selection
method='RF'    # Core modelling technique
size=1         # 1 for min, 2 for mid and 3 for max
nPerm=10       # Number of permutations (here set to 10 for illustration; normally set to 100 and inspected afterwards (see below))
permFit=numeric(nPerm)   # Allocate vector for permutation fitness
# Compute actual model and extract fitness metric
nCore=detectCores()-1
cl=makeCluster(nCore)
registerDoParallel(cl)
actual=MUVR(X=crispEM,ML=TRUE,nRep=nRep,nOuter=nOuter,varRatio=varRatio,method=method) 
actualFit=actual$miss[size]
for (p in 1:nPerm) {
  cat('\nPermutation',p,'of',nPerm)
  YPerm=sample(c(-1,1),size = nrow(crispEM),replace=TRUE) # Make permuted classes per individual
  perm=MUVR(X=crispEM,Y=YPerm,ML=TRUE,nRep=nRep,nOuter=nOuter,varRatio=varRatio,method=method)
  permFit[p]=perm$miss[size]
}
stopCluster(cl)
pPerm(actual = actualFit, h0 = permFit)
plotPerm(actual = actualFit, h0 = permFit) # Look at histogram to assess whether Gaussian in shape 
pPerm(actual = actualFit, h0 = permFit, type = 'non') # If not Gaussian, make a non-parametric test instead
plotPerm(actual = actualFit, h0 = permFit, type = 'non') # And plot with non-parametric p-value instead

