library(doParallel)
library(MUVR)


###################################################
# Regression example using "freelive" data

rm(list=ls())
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
plotMV(Regr_PLS_Quick,model='min')
plotVIP(Regr_PLS_Quick,model='min')
plotStability(Regr_PLS_Quick,model='min')

head(Regr_PLS_Quick$yPred) # min, mid and max predictions side-by-side
head(cbind(YR,Regr_PLS_Quick$yPred))

Regr_PLS_Quick$fitMetric   # Look at fitness metrics for min, mid and max models
Regr_PLS_Quick$nComp       # Number of components for min, mid and max models
Regr_PLS_Quick$nVar        # Number of variables for min, mid and max models

getVIRank(Regr_PLS_Quick,model='min')   # Extract most informative variables: Lower rank is better



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

head(Class_RF_Quick$yClass) # min, mid and max predictions side-by-side
head(cbind(Yotu,Class_RF_Quick$yClass))

Class_RF_Quick$nVar        # Number of variables for min, mid and max models

getVIRank(Class_RF_Quick,model='min')   # Extract most informative variables: Lower rank is better


###################################################
# Multilevel example using "crisp" data
rm(list=ls())
data("crisp")
#  crispEM: Effect matrix of difference (Treatment1 - Treatment2) in crossover intervention
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

head(ML_RF_Quick$yClass) # min, mid and max predictions side-by-side
head(cbind(ML_RF_Quick$inData$Y,ML_RF_Quick$yClass))

ML_RF_Quick$nVar        # Number of variables for min, mid and max models

getVIRank(ML_RF_Quick,model='min')   # Extract most informative variables: Lower rank is better

