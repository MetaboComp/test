rm(list=ls())
library(doParallel)
library(MUVR)
# Classification examples using "mosquito" data and PLS modelling
data("mosquito")
#  Xotu: Microbiota OTU (16S rDNA) from mosquitos captured in 3 different villages in Burkina Faso (Buck et al 2016)
#  Yotu: One of three villages of capture

#####################################
#
#  Stability over nRep

ncore=detectCores()-1 # 7 cores on my computer
cl=makeCluster(ncore)
registerDoParallel(cl)
Mosq=MUVR(X=Xotu,Y=Yotu,nRep=15*ncore,nOuter=6,varRatio=0.9,method='PLS') #
Mosq.75=MUVR(X=Xotu,Y=Yotu,nRep=15*ncore,nOuter=6,varRatio=0.75,method='PLS') #
Mosq.5=MUVR(X=Xotu,Y=Yotu,nRep=15*ncore,nOuter=6,varRatio=0.5,method='PLS') #
stopCluster(cl)

# Validation plots of Classification with different varRatios
par(mfrow=c(3,1))
plotVAL(Mosq)
plotVAL(Mosq.75)
plotVAL(Mosq.5)

# Stability with nRep
plotStability(Mosq)
plotStability(Mosq.75)
plotStability(Mosq.5)

# Function to investigate how variables are distributed in orignal and new models
tabVars=function(Org,New) {
  c1=varClass(Org)
  v1=data.frame(varsOrg=as.character(c(c1$Optimal,c1$Redundant,c1$Noisy)),classOrg=rep(c('Optimal','Redundant','Noisy'),c(length(c1$Optimal),length(c1$Redundant),length(c1$Noisy))))
  c2=varClass(New)
  v2=data.frame(varsNew=as.character(c(c2$Optimal,c2$Redundant,c2$Noisy)),classNew=rep(c('Optimal','Redundant','Noisy'),c(length(c2$Optimal),length(c2$Redundant),length(c2$Noisy))))
  v=merge(x=v1,y=v2,by.x='varsOrg',by.y='varsNew',all=T)
  if(any(is.na(v$classNew))){
    v$classNew=`levels<-`(addNA(v$classNew), c(levels(v$classNew), 'Not included'))
  }
  v$classOrg=droplevels(factor(v$classOrg,levels=c('Optimal','Redundant','Noisy')))
  v$classNew=droplevels(factor(v$classNew,levels=c('Optimal','Redundant','Noisy')))
  table(Original=v$classOrg,New=v$classNew,useNA='ifany')
}

# Look at classification of variables into Optimal, Redundant and Noisy as a function of varRatio
tabVars(Mosq,Mosq.75)
tabVars(Mosq,Mosq.5)
tabVars(Mosq.75,Mosq.5)



#####################################
#
#  Stability over varRatio

cl=makeCluster(ncore)
registerDoParallel(cl)
vRs=seq(0.5,0.95,by=0.05) # Different variable ratios
vR_pls1=list()
for (v in 1:length(vRs)) {
  vR=vRs[v]
  vR_pls1[[v]]=MUVR(X=Xotu,Y=Yotu,nRep=7*ncore,nOuter=6,varRatio=vR,method='PLS')
}
stopCluster(cl)

# plot calculation times
calcPLS=numeric(length(vRs))
for (v in 1:length(vRs)) {
  calcPLS[v]=vR_pls1[[v]]$calcMins
}
plot(vRs,calcPLS,log='y',ylim=c(1,35),xlab='varRatio',ylab='Calculation time (min)',type='l',main="Classification of 'Mosquito' data")

# Plot performance
vR_pls=vR_pls1
identical(vR_pls,vR_pls1)
nVar=round(vR_pls[[10]]$nVar[1])
VAll=names(sort(vR_pls[[10]]$VIP[,1])[1:nVar])
nV=VA=miss=numeric(10)
for (i in 1:length(vRs)) {
  nV[i]=round(mean(vR_pls[[i]]$nVar[1]))
  VA[i]=sum(names(sort(vR_pls[[i]]$VIP[,1])[1:nV[i]])%in%VAll)
  preds=vR_pls[[i]]$yClass[,1]
  miss[i]=sum(preds!=Yotu)
}
VA=VA/length(VAll)
par(mar=c(4,8,0,4)+.5)
plot(vRs,VA,type='l',ylim=c(0,1),lty=2,col='red',xlab='',ylab='')
title(ylab='Proportion of selected variables',line = 2,col.lab="red")
par(new=T)
plot(vRs,nV,ylim=c(0,85),type='l',xlab='varRatio (nRep=15)',ylab='',axes=F)
axis(2,line=4)
title(ylab='Number of selected variables',line = 6)
par(new=T)
plot(vRs,miss,ylim=c(0,length(Yotu)),type='l',col='blue',lty=3,axes=F,xlab='',ylab='')
axis(4)
mtext('Number of misclassifications',side=4,line = 2,col="blue",cex=par()$cex)



#####################################
#
#  Look at permutations of variables 


# Function to investigate how variables are distributed in model with original and permuted variables (starting with 'name') are cross-tabularized
tabVarsPerm=function(Org,Perm,name='P') {
  c1=varClass(Org)
  v1=data.frame(varsOrg=as.character(c(c1$Optimal,c1$Redundant,c1$Noisy)),classOrg=rep(c('Optimal','Redundant','Noisy'),c(length(c1$Optimal),length(c1$Redundant),length(c1$Noisy))))
  c2=varClass(Perm)
  v2=data.frame(varsPerm=as.character(c(c2$Optimal,c2$Redundant,c2$Noisy)),classPerm=rep(c('Optimal','Redundant','Noisy'),c(length(c2$Optimal),length(c2$Redundant),length(c2$Noisy))))
  vp=v2$varsPerm[substring(v2$varsPerm,1,1)==name]
  vp=data.frame(varsOrg=vp,classOrg="Permuted")
  v1=rbind(v1,vp)
  v=merge(x=v1,y=v2,by.x='varsOrg',by.y='varsPerm',all=T)
  if(any(is.na(v$classPerm))){
    v$classPerm=`levels<-`(addNA(v$classPerm), c(levels(v$classPerm), 'Not included'))
  }
  v$classOrg=droplevels(factor(v$classOrg,levels=c('Optimal','Redundant','Noisy','Permuted')))
  v$classPerm=droplevels(factor(v$classPerm,levels=c('Optimal','Redundant','Noisy','Not included')))
  table(Original=v$classOrg,Permuted=v$classPerm,useNA='ifany')
}

plPerms=function(Org,Perm,name='P') {
  whP=grep(name,varClass(Perm)$Optimal,value=T)
  for(i in 1:length(whP)) {
    XO=Org$inData$X
    XP=Perm$inData$X
    permVar=subset(XP,select=whP[i])[,1]
    orgname=substring(whP[i],nchar(name)+1)
    orgVar=subset(XO,select=orgname)[,1]
    par(mfrow=c(1,1))
    plot(orgVar,permVar,xlab='Original',ylab='Permuted',main=orgname)
    legend('topright',paste('r=',signif(cor(orgVar,permVar),3)),bty='n')
    par(mfrow=c(2,1))
    plot(Yotu,orgVar,ylab='original',main=orgname)
    plot(Yotu,permVar,ylab='permuted',main=paste('Permuted',orgname))
  }
}

# Variable permutations PLS classification:
# O for Optimal, R for Redundant, N for Noise, P for permutation
Morg=Mosq
mclass=varClass(Morg)
O=mclass$Optimal
R=mclass$Redundant
N=mclass$Noisy

# 1st permutation: Addition of randomised "O" -> N?
P1=subset(Xotu,select = O)
colnames(P1)=paste('P1',colnames(P1),sep='_')
for(i in 1:ncol(P1)) P1[,i]=sample(P1[,i])
V1=cbind(Xotu,P1)
cl=makeCluster(ncore)
registerDoParallel(cl)
MP1P=MUVR(X=V1,Y=Yotu,nRep=7*ncore,nOuter=6,varRatio=0.9,method='PLS') # More proper model - Also more time consuming
stopCluster(cl)
tabVarsPerm(Morg,MP1P)
# plPerms(Morg,MP1P,'P1_') # If any permuted variables are (by random) classified as Optimal 

# 2nd permutation: Randomise N -> N?
P2=subset(Xotu,select = N)
colnames(P2)=paste('P2',colnames(P2),sep='_')
for(i in 1:ncol(P2)) P2[,i]=sample(P2[,i])
V2=cbind(subset(Xotu,select = c(O,R)),P2)
cl=makeCluster(ncore)
registerDoParallel(cl)
MP2P=MUVR(X=V2,Y=Yotu,nRep=7*ncore,nOuter=6,varRatio=0.9,method='PLS') # More proper model - Also more time consuming
stopCluster(cl)
tabVarsPerm(Morg,MP2P)
plPerms(Morg,MP2P,'P2_') # If any permuted variables are (by random) classified as Optimal 

# 3rd permutation: Randomise R -> N?
P3=subset(Xotu,select = R)
colnames(P3)=paste('P3',colnames(P3),sep='_')
for(i in 1:ncol(P3)) P3[,i]=sample(P3[,i])
V3=cbind(subset(Xotu,select = c(O,N)),P3)
cl=makeCluster(ncore)
registerDoParallel(cl)
MP3P=MUVR(X=V3,Y=Yotu,nRep=7*ncore,nOuter=6,varRatio=0.9,method='PLS') # More proper model - Also more time consuming
stopCluster(cl)
tabVarsPerm(Morg,MP3P)
plPerms(Morg,MP3P,'P3_') # If any permuted variables are (by random) classified as Optimal 

# 4th permutation: Randomise O -> N?; R -> O?
P4=subset(Xotu,select = O)
colnames(P4)=paste('P4',colnames(P4),sep='_')
for(i in 1:ncol(P4)) P4[,i]=sample(P4[,i])
V4=cbind(subset(Xotu,select = c(R,N)),P4)
cl=makeCluster(ncore)
registerDoParallel(cl)
MP4P=MUVR(X=V4,Y=Yotu,nRep=7*ncore,nOuter=6,varRatio=0.9,method='PLS') # More proper model - Also more time consuming
stopCluster(cl)
tabVarsPerm(Morg,MP4P)
plPerms(Morg,MP4P,'P4_') # If any permuted variables are (by random) classified as Optimal 
