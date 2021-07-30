rm(list=ls())
library(doParallel)
library(MUVR)
# Classification examples using "mosquito" data and PLS modelling
data("mosquito")
#  Xotu: Microbiota OTU (16S rDNA) from mosquitos captured in 3 different villages in Burkina Faso (Buck et al 2016)
#  Yotu: One of three villages of capture

#####################################
# For classification
#  Stability over nRep

ncore=detectCores()-1 # 7 cores on my computer
cl=makeCluster(ncore)
registerDoParallel(cl)

Mosq=MUVR(X=Xotu,Y=Yotu,nRep=15*ncore,nOuter=6,varRatio=0.9,method='PLS') #varRatio=0.9
Mosq.75=MUVR(X=Xotu,Y=Yotu,nRep=15*ncore,nOuter=6,varRatio=0.75,method='PLS') #varRatio=0.75
Mosq.5=MUVR(X=Xotu,Y=Yotu,nRep=15*ncore,nOuter=6,varRatio=0.5,method='PLS') #varRatio=0.5

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
###label each variable as optimal redundant or nosisy
tabVars=function(Org,New) {
  c1=varClass(Org)
  v1=data.frame(varsOrg=as.character(c(c1$Optimal,
                                       c1$Redundant,
                                       c1$Noisy)),
                classOrg=rep(c('Optimal','Redundant','Noisy'),
                             c(length(c1$Optimal),
                               length(c1$Redundant),
                               length(c1$Noisy))))
  c2=varClass(New)
  v2=data.frame(varsNew=as.character(c(c2$Optimal,
                                       c2$Redundant,
                                       c2$Noisy)),
                classNew=rep(c('Optimal','Redundant','Noisy'),
                             c(length(c2$Optimal),
                               length(c2$Redundant),
                               length(c2$Noisy))))

##########Real merge start
  v=merge(x=v1,    ###x and y have the equal position
          y=v2,
          by.x='varsOrg', ##to match this column with columns in y
          by.y='varsNew', ##to match this column with columns in x
##by.x,by.y specifications of the columns used for merging
##By default the data frames are merged on the columns with names they both have,
##but separate specifications of the columns can be given by by.x and by.y. The rows in the two data frames that match on the specified columns are extracted, and joined together. If there is more than one match, all possible matches contribute one row each.

          all=T)  ##only common nams is included in the merge data frame
##logical; all = L is shorthand for all.x = L and all.y = L, where L is either TRUE or FALSE.
##all.x	logical; if TRUE, then extra rows will be added to the output, one for each row in x that has no matching
##row in y. These rows will have NAs in those columns that are usually filled with values from y.
##The default is FALSE, so that only rows with data from both x and y are included in the output.


  if(any(is.na(v$classNew))){                             ####add this in the level structure of NA
    v$classNew=`levels<-`(addNA(v$classNew),   ###include NA in the level output
                          c(levels(v$classNew), 'Not included')) ####new level that is called "not included" is added
  }
  v$classOrg=droplevels(factor(v$classOrg,
                               levels=c('Optimal','Redundant','Noisy')))
  v$classNew=droplevels(factor(v$classNew,
                               levels=c('Optimal','Redundant','Noisy')))

  #######droplevel is to drop unused level from a factor

  table(Original=v$classOrg,
        New=v$classNew,
        useNA='ifany')
###useNA
  ##if no,ignores NA value in the table inputs
  ##if ifany, includes NA value in the table
  ##if always, adds an extra NA value added to each factor level if the NA does not really exist
###controls if the table includes counts of NA values: the allowed values correspond to never ("no"),
###only if the count is positive ("ifany") and even for zero counts ("always")
}

# Look at classification of variables into Optimal, Redundant and Noisy as a function of varRatio
tabVars(Mosq,Mosq.75)
tabVars(Mosq,Mosq.5)
tabVars(Mosq.75,Mosq.5)


##############################################################################################################
#############################################################################################################
#
#  Stability over varRatio

cl=makeCluster(ncore)
registerDoParallel(cl)

vRs=seq(0.5,0.95,by=0.05) # Different variable ratios

vR_pls1=list()

for (v in 1:length(vRs)) {
  vR=vRs[v]
  vR_pls1[[v]]=MUVR(X=Xotu,Y=Yotu,nRep=7*ncore,nOuter=6,varRatio=vR,method='PLS')   ####a list for MUVR objects
}
stopCluster(cl)

#######################################################################################
# plot calculation times
calcPLS=numeric(length(vRs))

for (v in 1:length(vRs)) {
  calcPLS[v]=vR_pls1[[v]]$calcMins
}
plot(vRs,
     calcPLS,
     log='y',
     ylim=c(1,35),
     xlab='varRatio',
     ylab='Calculation time (min)',
     type='l',
     main="Classification of 'Mosquito' data")

#########################################################################################
# Plot performance
vR_pls=vR_pls1

identical(vR_pls,vR_pls1)

nVar=round(vR_pls[[10]]$nVar[1])   ###The 10th MUVR subject under min model

VAll=names(sort(vR_pls[[10]]$VIP[,1])[1:nVar])    ####names of the variables kept in the 10th MUVR model

nV=VA=miss=numeric(10)     ### a vector of length of 10

for (i in 1:length(vRs)) {          ###length of MUVR
  nV[i]=round(mean(vR_pls[[i]]$nVar[1]))                            ##### a vector of number of variables kept in each MUVR model
  VA[i]=sum(names(sort(vR_pls[[i]]$VIP[,1])[1:nV[i]])%in%VAll)      ##### variables kept that is included in VALL
  preds=vR_pls[[i]]$yClass[,1]                              #####this is for classification problems regardless of method
  miss[i]=sum(preds!=Yotu)                                      ###miss classification number of
}

VA=VA/length(VAll)
par(mar=c(4,8,0,4)+.5)
plot(vRs,                  ####0.05,0.1,0.15......
     VA,                   ####percentage of variables included in the 10th MUVR model for each model
     type='l',
     ylim=c(0,1),
     lty=2,               ##line type
     col='red',
     xlab='',
     ylab='')
title(ylab='Proportion of selected variables',
      line = 2,      ####specifying a value for line overrides the default placement of labels,
                     ###and places them this many lines outwards from the plot edge.
      col.lab="red")

par(new=T)   ###This is to add the following plot to the plot below


plot(vRs,         ####0.05,0.1,0.15......
     nV,          ####number of selected variables
     ylim=c(0,85),
     type='l',
     xlab='varRatio (nRep=15)',
     ylab='',
     axes=F)     ###remove axis

axis(2,         ###y axis
     line=4)    ###The position of the line

title(ylab='Number of selected variables',
      line = 6)                          ####position of the line

par(new=T)

plot(vRs,
     miss,
     ylim=c(0,length(Yotu)),
     type='l',
     col='blue',
     lty=3,
     axes=F,
     xlab='',
     ylab='')
axis(4)
###	an integer specifying which side of the plot the axis is to be drawn on. The axis is placed as follows: 1=below, 2=left, 3=above and 4=right.
mtext('Number of misclassifications',
      side=4,      ####on which side of the plot (1=bottom, 2=left, 3=top, 4=right).
      line = 2,    ########on which MARgin line, starting at 0 counting outwards.
      col="blue",
      cex=par()$cex)   #### ###custom text size, if this is not added, the text size is bigger



#####################################################################################
####################################################################################
#  Look at permutations of variables in X


##Function to investigate how variables are distributed in model with original and permuted variables
##(starting with 'name') are cross-tabularized
tabVarsPerm=function(Org,
                     Perm,
                     name='P') {
  ####This could only be understood after see after the code line 295 where permutation variables are added with P
  ####What is the original variables has the name P inside
  c1=varClass(Org)
  v1=data.frame(varsOrg=as.character(c(c1$Optimal,c1$Redundant,c1$Noisy)),
                classOrg=rep(c('Optimal','Redundant','Noisy'),
                             c(length(c1$Optimal),length(c1$Redundant),length(c1$Noisy))))
  ###a data frame with 2 columns, one is variable names , the other is type(optimal redundant, noisy)
  c2=varClass(Perm)
  v2=data.frame(varsPerm=as.character(c(c2$Optimal,c2$Redundant,c2$Noisy)),
                classPerm=rep(c('Optimal','Redundant','Noisy'),
                              c(length(c2$Optimal),length(c2$Redundant),length(c2$Noisy))))

  vp=v2$varsPerm[substring(v2$varsPerm,1,1)==name]    ####Extract or replace sub strings in a character vector.
  ## Find variables that are permuted

  vp=data.frame(varsOrg=vp,                    ###vp is a new data frame of 2 column, vp of 0 value, class Org is permuted
                classOrg="Permuted")

  v1=rbind(v1,vp)         ###original variables and permuted variables

  v=merge(x=v1,                               ###a new data frame of 3 columns, name, classOrg, classPerm
          y=v2,
          by.x='varsOrg',
          by.y='varsPerm',
          all=T)

  if(any(is.na(v$classPerm))){
    v$classPerm=`levels<-`(addNA(v$classPerm),
                           c(levels(v$classPerm), 'Not included'))
  }

  v$classOrg=droplevels(factor(v$classOrg,                                                 ###It has more values since its rbind the values of
                               levels=c('Optimal','Redundant','Noisy','Permuted')))

  v$classPerm=droplevels(factor(v$classPerm,
                                levels=c('Optimal','Redundant','Noisy','Not included')))
  table(Original=v$classOrg,
        Permuted=v$classPerm,
        useNA='ifany')           ###NA is included in the table
}




########One plot for one variable
########plot Perm
plPerms=function(Org,
                 Perm,
                 name='P') {
  whP=grep(name,
           varClass(Perm)$Optimal,  ##What if the variable itself has the value of P
           value=T)   ###(value = TRUE) returns a character vector containing the selected elements of x (after coercion, preserving names but no other attributes).
 ###grep("[a,C,x,z]", letters,value=T)      [1] "a" "x" "z"
 ### grep("[a,C,x,z]", letters,value=F)     [1]  1 24 26

###For each variable that is in Perm$Optimal
   for(i in 1:length(whP)) {
    XO=Org$inData$X
    XP=Perm$inData$X
    permVar=subset(XP,select=whP[i])[,1]   ### the first column  only one column
    orgname=substring(whP[i],nchar(name)+1)  ###output names from the first letter after "P"
    ##x <- c("asfef", "qwerty", "yuiop[", "b", "stuff.blah.yech")
    ###nchar(x)
    # 5  6  6  1 15
    ###substring("abcdefg",3)
    ###[1] "cdefg"
    orgVar=subset(XO,select=orgname)[,1]   ###variables  that is in Perm$optimal . select those variables in Original data the first column

    par(mfrow=c(1,1))
####plot X1 and X2 for each vairables
    plot(orgVar,
         permVar,
         xlab='Original',
         ylab='Permuted',
         main=orgname)
    legend('topright',
           paste('r=',
                 signif(cor(orgVar,permVar),3)),    ## round to 3 digit
          bty='n')       ##border type
    par(mfrow=c(2,1))
####plot X1 Y and X2 Y for each variable
     plot(Yotu,
         orgVar,
         ylab='original',
         main=orgname)
    plot(Yotu,
         permVar,
         ylab='permuted',
         main=paste('Permuted',orgname))
  }
}



###############################################################################################
# Variable permutations PLS classification:
# O for Optimal, R for Redundant, N for Noise, P for permutation
Morg=Mosq
mclass=varClass(Morg)   ###2 column one column variable name, the other column is name of class(optimal, redundant,nosiy)
O=mclass$Optimal        ###O is the name of optimal variables
R=mclass$Redundant      ###R is the name of redundant variables
N=mclass$Noisy          ###N is the name of noisy variables

# 1st permutation: Addition of randomised "O" -> N?

P1=subset(Xotu,select = O)   ###variable columns that is included in O
colnames(P1)=paste('P1',colnames(P1),sep='_')

for(i in 1:ncol(P1)) P1[,i]=sample(P1[,i])   ###all the variables are resampled
#########################################################################################
V1=cbind(Xotu,P1)  ###########################################P1 variables is sampled 2 times?
####################################################################################
cl=makeCluster(ncore)
registerDoParallel(cl)

MP1P=MUVR(X=V1,     ##cbind Xotu and O vairables
          Y=Yotu,
          nRep=7*ncore,
          nOuter=6,
          varRatio=0.9,
          method='PLS') # More proper model - Also more time consuming
stopCluster(cl)

tabVarsPerm(Morg,MP1P)
# plPerms(Morg,MP1P,'P1_') # If any permuted variables are (by random) classified as Optimal

########### 2nd permutation: Randomise N -> N?
P2=subset(Xotu,select = N)
colnames(P2)=paste('P2',colnames(P2),sep='_')

for(i in 1:ncol(P2)) P2[,i]=sample(P2[,i])
V2=cbind(subset(Xotu,select = c(O,R)),P2)

cl=makeCluster(ncore)
registerDoParallel(cl)
MP2P=MUVR(X=V2,
          Y=Yotu,
          nRep=7*ncore,
          nOuter=6,
          varRatio=0.9,
          method='PLS') # More proper model - Also more time consuming
stopCluster(cl)

tabVarsPerm(Morg,MP2P)

plPerms(Morg,
        MP2P,
        'P2_') # If any permuted variables are (by random) classified as Optimal

########### 3rd permutation: Randomise R -> N?
P3=subset(Xotu,select = R)
colnames(P3)=paste('P3',colnames(P3),sep='_')
for(i in 1:ncol(P3)) P3[,i]=sample(P3[,i])

V3=cbind(subset(Xotu,select = c(O,N)),P3)

cl=makeCluster(ncore)
registerDoParallel(cl)
MP3P=MUVR(X=V3,
          Y=Yotu,
          nRep=7*ncore,
          nOuter=6,
          varRatio=0.9,
          method='PLS') # More proper model - Also more time consuming
stopCluster(cl)

tabVarsPerm(Morg,MP3P)

plPerms(Morg,
        MP3P,
        'P3_') # If any permuted variables are (by random) classified as Optimal


############ 4th permutation: Randomise O -> N?; R -> O?
P4=subset(Xotu,select = O)
colnames(P4)=paste('P4',colnames(P4),sep='_')
for(i in 1:ncol(P4)) P4[,i]=sample(P4[,i])
V4=cbind(subset(Xotu,select = c(R,N)),P4)

cl=makeCluster(ncore)
registerDoParallel(cl)
MP4P=MUVR(X=V4,
          Y=Yotu,
          nRep=7*ncore,
          nOuter=6,
          varRatio=0.9,
          method='PLS') # More proper model - Also more time consuming
stopCluster(cl)

tabVarsPerm(Morg,MP4P)

plPerms(Morg,
        MP4P,
        'P4_') # If any permuted variables are (by random) classified as Optimal
