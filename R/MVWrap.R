#' MVWrap: Wrapper for Multivariate modelling with Variable selection
#' 
#' Repeated double cross validation with tuning of variables in the inner loop.
#' @param X Independent variables. NB: Variables (columns) must have names/unique identifiers. NAs not allowed in data.
#' @param Y Response vector (Dependent variable). For PLS-DA, values should be -1 vs 1 for the two classes.
#' @param ID Subject identifier (for sampling by subject; Assumption of independence if not specified)
#' @param nRep Number of repetitions of double CV..
#' @param nOuter Number of outer CV loop segments.
#' @param nInner Number of inner CV loop segments.
#' @param varRatio Ratio of variables to include in subsequent inner loop iteration.
#' @param DA Logical for Classification (discriminant analysis) (Defaults do FALSE, i.e. regression). PLS is limited to two-class problems (see `Y` above).
#' @param fitness Fitness function for model tuning (choose either 'AUROC' or 'misClass' for classification; or 'RMSEP' (default) for regression.)
#' @param method Multivariate method. Supports 'PLS' and 'RF' (default)
#' @param methParam List with parameter settings for specified MV method (defaults to ???)
#' @param ML Logical for multilevel analysis (defaults to FALSE)
#' @param modReturn Logical for returning outer segment models (defaults to FALSE)
#' @param newdata New data matrix ONLY for prediction NOT modelling
#' @return An object containing stuff...
#' @export
MVWrap=function(X,Y,ID,nRep=5,nOuter=6,nInner,varRatio=0.75,DA=FALSE,fitness=c('AUROC','misClass','RMSEP'),method=c('PLS','RF'),methParam,ML=FALSE,pred=FALSE,newdata=NULL){
  # Start timer
  start.time=proc.time()[3]
  # Check indata
  if (is.null(dim(X))){
    cat('\nError: Wrong format of X matrix.\n')
    return(NULL)
  }
  nSamp=nrow(X)
  nVar=nVar0=ncol(X)
  if (missing(ID)) ID=1:nSamp
  if (missing(nInner)) nInner=nOuter-1
  if (missing(method)) method='RF'
  if (missing(methParam)) {
    if (method=='PLS') {
      methParam=list(compMax=ifelse(nVar<3,nVar,3),mode='regression')
    } else {
      methParam=NULL
    }
  }
  if (ML) {
    X=rbind(X,-X)
    Y=rep(c(-1,1),each=nSamp)
    nSamp=2*nSamp
    ID=c(ID,ID)
  }
  if (missing(fitness)) {
    if (DA) {
      fitness='misClass'
    } else {
      fitness='RMSEP'
    }
  }
  pred=ifelse(missing(newdata),FALSE,TRUE)
  if (pred) {  # If model is also used to predict newdata matrix, allocate a YP response array
    YP=matrix(nrow=nrow(newdata),ncol=nRep*nOuter)  # Predictions for all samples (rows) across outer segments (cols) and repetitions (dim3)
  }
  if (nrow(X)!=length(Y)) {
    cat('\nError: Must have same nSamp in X and Y.\n')
    return(NULL)
  }
  ## Sort sampling based in subjects and not index
  unik=!duplicated(ID)  # boolean of unique IDs
  unikID=ID[unik]  
  if (DA) {
    unikY=Y[unik]  # Counterintuitive, but needed for groupings by Ynames
    Ynames=sort(unique(Y))  # Find groups
    groups=length(Ynames) # Number of groups
    groupID=list()  # Allocate list for indices of groups
    for (g in 1:groups) { 
      groupID[[g]]=unikID[unikY==Ynames[g]]  # Find indices per group
    }
  }
  ### Allocate variables
  VIPMin=VIPMid=VIPMax=matrix(nrow=nVar0,ncol=3) # Allocate matrix for variable importance averaged over repetition
  rownames(VIPMin)=rownames(VIPMid)=rownames(VIPMax)=colnames(X)
  colnames(VIPMin)=colnames(VIPMid)=colnames(VIPMax)=c('Mean','SD','CV')
  # Allocate matrix for prediction of outer segments Y per repetition
  yPredMin=yPredMid=yPredMax=matrix(nrow=length(Y),ncol=nRep)
  colnames(yPredMin)=colnames(yPredMid)=colnames(yPredMax)=paste('Rep',1:nRep,sep='')
  yPredMinR=yPredMidR=yPredMaxR=numeric(length(Y))
  # Allocate response vectors and matrices for var's, nComp and VIP ranks over repetitions
  varRepMin=varRepMid=varRepMax=nCompRepMin=nCompRepMid=nCompRepMax=missRep=numeric(nRep)
  names(varRepMin)=names(varRepMid)=names(varRepMax)=names(nCompRepMin)=names(nCompRepMid)=names(nCompRepMax)=names(missRep)=paste(rep('rep',nRep),1:nRep,sep='')
  VIPRepMin=VIPRepMid=VIPRepMax=matrix(data=nVar0,nrow=nVar0,ncol=nRep)
  rownames(VIPRepMin)=rownames(VIPRepMid)=rownames(VIPRepMax)=colnames(X)
  colnames(VIPRepMin)=colnames(VIPRepMid)=colnames(VIPRepMax)=paste(rep('rep',nRep),1:nRep,sep='')
  var=numeric()
  cnt=0
  while (nVar>2) {  
    cnt=cnt+1
    var=c(var,nVar)
    nVar=floor(varRatio*nVar)
  }
  ## Choose package/core algorithm according to chosen method
  packs=c(ifelse(method=='PLS','mixOmics','randomForest'),'pROC')
  ## Start repetitions
  reps=foreach(r=1:nRep, .packages=packs, .export='vectSamp') %dopar% {
    # r=1
    # r=r+1
    sink('log.txt',append=TRUE)
    if (pred) YPR=matrix(nrow=nrow(newdata),ncol=nOuter)
    cat('\n','   Repetition ',r,' of ',nRep,':',sep='')
    if (DA) {
      groupTest=list()  ## Allocate list for samples within group
      for (gT in 1:groups) { 
        groupTest[[gT]]=vectSamp(groupID[[gT]],n=nOuter)  # Draw random samples within group
      }
      allTest=groupTest[[1]] # Add 1st groups to 'Master' sample of all groups
      for (gT in 2:groups) {  # Add subsequent groups
        allTest=allTest[order(sapply(allTest,length))]
        for (aT in 1:nOuter) {
          allTest[[aT]]=sort(c(allTest[[aT]],groupTest[[gT]][[aT]]))
        }
      }
    } else {
      allTest=vectSamp(unikID,n=nOuter)
    }
    varOutMin=varOutMid=varOutMax=nCompOutMin=nCompOutMid=nCompOutMax=numeric(nOuter)
    names(varOutMin)=names(varOutMid)=names(varOutMax)=names(nCompOutMin)=names(nCompOutMax)=paste(rep('outSeg',nOuter),1:nOuter,sep='')
    VIPOutMin=VIPOutMid=VIPOutMax=matrix(data=nVar0,nrow=nVar0,ncol=nOuter)
    rownames(VIPOutMin)=rownames(VIPOutMid)=rownames(VIPOutMax)=colnames(X)
    colnames(VIPOutMin)=colnames(VIPOutMid)=colnames(VIPOutMax)=paste(rep('outSeg',nOuter),1:nOuter,sep='')
    ## Perform outer loop segments -> one "majority vote" MV model per segment
    for (i in 1:nOuter) {   
      # i=1
      # i=i+1
      cat('\n Segment ',i,' (variables):',sep='') # Counter
      ## Draw out test set
      testID=allTest[[i]] # Draw out segment = holdout set BASED ON UNIQUE ID
      testIndex=ID%in%testID # Logical for samples included
      xTest=X[testIndex,]
      yTest=Y[testIndex]
      inID=unikID[!unikID%in%testID]  # IDs not in test set
      inY=unikY[!unikID%in%testID]  # Counterintuitive, but needed for grouping by Ynames
      ## Allocate variables for later use
      missIn=aucIn=rmsepIn=nCompIn=matrix(nrow=nInner,ncol=cnt)
      rownames(rmsepIn)=rownames(missIn)=rownames(aucIn)=rownames(nCompIn)=paste(rep('inSeg',nInner),1:nInner,sep='')
      colnames(rmsepIn)=colnames(missIn)=colnames(aucIn)=colnames(nCompIn)=var
      VIPInner=array(data=nVar0,dim=c(nVar0,cnt,nInner))
      rownames(VIPInner)=colnames(X)
      colnames(VIPInner)=var
      dimnames(VIPInner)[[3]]=paste(rep('inSeg',nInner),1:nInner,sep='')
      # Restart variables
      incVar=colnames(X)
      ## Perform steps with successively fewer features
      for (count in 1:cnt) {  # Build models with successively fewer variables. Quality metric = number of missclassifications for Validation set
        # count=1 # for testing
        # count=count+1
        nVar=var[count]
        cat(nVar)
        # if (resampInner==TRUE) inSamp=vectSamp(inInd,n=nInner)  # Resample inner segments for each variable elimination step
        if (method=='PLS') comp=ifelse(nVar<methParam$compMax,nVar,methParam$compMax)
        if (DA==T) {
          groupIDVal=list()
          for (g in 1:groups) { 
            groupIDVal[[g]]=inID[inY==Ynames[g]]  # Find indices per group
          }
          groupVal=list()  ## Allocate list for samples within group
          for (gV in 1:groups) { 
            groupVal[[gV]]=vectSamp(groupIDVal[[gV]],n=nInner)  # Draw random samples within group
          }
          allVal=groupVal[[1]] # Add 1st groups to 'Master' sample of all groups
          for (gV in 2:groups) {  # Add subsequent groups
            allVal=allVal[order(sapply(allVal,length))]
            for (aV in 1:nInner) {
              allVal[[aV]]=sort(c(allVal[[aV]],groupVal[[gV]][[aV]]))
            }
          }
        } else {
          allVal=vectSamp(inID,n=nInner)
        }
        ## Inner CV loop
        for (j in 1:nInner) {
          # j=1 
          # j=j+1
          cat('.') # Counter
          valID=allVal[[j]] # Draw out segment = validation set
          valIndex=ID%in%valID
          xVal=X[valIndex,]
          xVal=subset(xVal,select=incVar)
          yVal=Y[valIndex]
          trainID=inID[!inID%in%valID]
          trainIndex=ID%in%trainID # Define Training segment
          xTrain=X[trainIndex,]
          xTrain=subset(xTrain,select=incVar)
          yTrain=Y[trainIndex]
          # sum(trainIndex,valIndex,testIndex)
          # trainIndex|valIndex|testIndex
          ## Make inner model
          if (method=='PLS') {
            inMod=plsInner(xTrain,yTrain,xVal,yVal,fitness,comp,methParam$mode)
            nCompIn[j,count]=inMod$nComp
          } else {
            cat('Inner RFMod yet to be coded')
            return(NULL)
            inReturn=rfInner(xTrain,yTrain,xVal,yVal,fitness,methParam)
          }
          # Store fitness metric
          if (fitness=='misClass') {
            missIn[j,count]=inMod$miss
          } else if (fitness=='AUROC') {
            aucIn[j,count]=inMod$auc
          } else {
            rmsepIn[j,count]=inMod$rmsep
          }
          # Store VIPs
          VIPInner[match(names(inMod$vip),rownames(VIPInner)),count,j]=inMod$vip
        }
        ## Average inner VIP ranks before variable elimination
        VIPInAve=apply(VIPInner[,count,],1,mean)
        if (count<cnt) {
          incVar=names(VIPInAve[order(VIPInAve)])[1:var[count+1]]
        }
      }
      if (fitness=='AUROC') {
        minIndex=max(which(apply(t(apply(aucIn,1,rank)),2,mean)==max(apply(t(apply(aucIn,1,rank)),2,mean))))
        maxIndex=min(which(apply(t(apply(aucIn,1,rank)),2,mean)==max(apply(t(apply(aucIn,1,rank)),2,mean))))
        # midIndex=which.min(abs(var-mean(c(var[minIndex],var[maxIndex]))))  # Arithmetic mean
        midIndex=which.min(abs(var-exp(mean(log(c(var[minIndex],var[maxIndex]))))))  # Geometric mean
      } else if (fitness=='misClass') {
        minIndex=max(which(apply(t(apply(missIn,1,rank)),2,mean)==min(apply(t(apply(missIn,1,rank)),2,mean))))
        maxIndex=min(which(apply(t(apply(missIn,1,rank)),2,mean)==min(apply(t(apply(missIn,1,rank)),2,mean))))
        # midIndex=which.min(abs(var-mean(c(var[minIndex],var[maxIndex]))))  # Arithmetic mean
        midIndex=which.min(abs(var-exp(mean(log(c(var[minIndex],var[maxIndex]))))))  # Geometric mean
      }else {
        minIndex=max(which(apply(t(apply(rmsepIn,1,rank)),2,mean)==min(apply(t(apply(rmsepIn,1,rank)),2,mean))))
        maxIndex=min(which(apply(t(apply(rmsepIn,1,rank)),2,mean)==min(apply(t(apply(rmsepIn,1,rank)),2,mean))))
        # midIndex=which.min(abs(var-mean(c(var[minIndex],var[maxIndex]))))  # Arithmetic mean
        midIndex=which.min(abs(var-exp(mean(log(c(var[minIndex],var[maxIndex]))))))  # Geometric mean
      }
      # Per outer segment: Average inner loop variables, nComp and VIP ranks 
      varOutMin[i]=var[minIndex]
      varOutMid[i]=var[midIndex]
      varOutMax[i]=var[maxIndex]
      if (method=='PLS') {
        nCompOutMin[i]=round(mean(nCompIn[,minIndex]))
        nCompOutMid[i]=round(mean(nCompIn[,midIndex]))
        nCompOutMax[i]=round(mean(nCompIn[,maxIndex]))
      }
      VIPOutMin[,i]=apply(VIPInner[,minIndex,],1,mean)
      VIPOutMid[,i]=apply(VIPInner[,midIndex,],1,mean)
      VIPOutMax[,i]=apply(VIPInner[,maxIndex,],1,mean)
      # Build outer model for min and max nComp and predict YTEST
      xIn=X[!testIndex,] # Perform Validation on all samples except holdout set
      yIn=Y[!testIndex]
      incVarMin=rownames(VIPOutMin)[rank(VIPOutMin[,i])<=varOutMin[i]]
      incVarMid=rownames(VIPOutMid)[rank(VIPOutMid[,i])<=varOutMid[i]]
      incVarMax=rownames(VIPOutMax)[rank(VIPOutMax[,i])<=varOutMax[i]]
      if (method=='PLS'){
        plsOutMin=pls(subset(xIn,select=incVarMin),yIn,ncomp=nCompOutMin[i],mode="classic")
        yPredMinR[testIndex]=predict(plsOutMin,newdata=subset(xTest,select=incVarMin))$predict[,,nCompOutMin[i]]  # 	
        plsOutMid=pls(subset(xIn,select=incVarMid),yIn,ncomp=nCompOutMid[i],mode="classic")
        yPredMidR[testIndex]=predict(plsOutMid,newdata=subset(xTest,select=incVarMid))$predict[,,nCompOutMid[i]]  # 	
        plsOutMax=pls(subset(xIn,select=incVarMax),yIn,ncomp=nCompOutMax[i],mode="classic")
        yPredMaxR[testIndex]=predict(plsOutMax,newdata=subset(xTest,select=incVarMax))$predict[,,nCompOutMax[i]]  # 	
        # If models complain about near zero variance, this code can be adapted and added before prediction!
        # if (length(plsOutMax$nzv$Position)>0) {
          # removeVar=rownames(plsOutMax$nzv$Metrics)
          # xTestMax=xTestMax[,!colnames(xTestMax)%in%removeVar]
        # }
        if (pred) {  # Predict extra prediction samples (newdata)
          YPR[,i]=predict(plsOutMid,newdata=subset(newdata,select=incVarMid))$predict[,,nCompOutMid[i]]
        }
      } else {
        cat('Outer RFMod yet to be coded')
        return(NULL)
      }
    }
    # Per repetition: Average outer loop variables, nComp and VIP ranks 
    parReturn=list(yPredMin=yPredMinR,yPredMid=yPredMidR,yPredMax=yPredMaxR)
    parReturn$varRepMin=round(mean(varOutMin))
    parReturn$varRepMid=round(mean(varOutMid))
    parReturn$varRepMax=round(mean(varOutMax))
    parReturn$VIPRepMin=apply(VIPOutMin,1,mean)
    parReturn$VIPRepMid=apply(VIPOutMid,1,mean)
    parReturn$VIPRepMax=apply(VIPOutMax,1,mean)
    if (method=='PLS'){
      parReturn$nCompRepMin=round(mean(nCompOutMin))
      parReturn$nCompRepMid=round(mean(nCompOutMid))
      parReturn$nCompRepMax=round(mean(nCompOutMax))
    }
    if (pred) parReturn$YP=YPR
    return(parReturn)
  }
  for (r in 1:nRep) {
    yPredMin[,r]=reps[[r]]$yPredMin
    varRepMin[r]=reps[[r]]$varRepMin
    nCompRepMin[r]=reps[[r]]$nCompRepMin
    VIPRepMin[,r]=reps[[r]]$VIPRepMin
    yPredMid[,r]=reps[[r]]$yPredMid
    varRepMid[r]=reps[[r]]$varRepMid
    nCompRepMid[r]=reps[[r]]$nCompRepMid
    VIPRepMid[,r]=reps[[r]]$VIPRepMid
    yPredMax[,r]=reps[[r]]$yPredMax
    varRepMax[r]=reps[[r]]$varRepMax
    nCompRepMax[r]=reps[[r]]$nCompRepMax
    VIPRepMax[,r]=reps[[r]]$VIPRepMax
    if (pred) YP[,(nOuter*(r-1)+1):(nOuter*r)]=reps[[r]]$YP
  }
  yMin=yMid=yMax=matrix(nrow=nrow(X),ncol=3)
  colnames(yMin)=colnames(yMid)=colnames(yMax)=c('Mean','SD','CV')
  yMin[,1]=apply(yPredMin,1,mean)[1:nrow(X)]
  yMin[,2]=apply(yPredMin,1,sd)[1:nrow(X)]
  yMin[,3]=abs(100*yMin[,2]/yMin[,1])
  yMid[,1]=apply(yPredMid,1,mean)[1:nrow(X)]
  yMid[,2]=apply(yPredMid,1,sd)[1:nrow(X)]
  yMid[,3]=abs(100*yMid[,2]/yMid[,1])
  yMax[,1]=apply(yPredMax,1,mean)[1:nrow(X)]
  yMax[,2]=apply(yPredMax,1,sd)[1:nrow(X)]
  yMax[,3]=abs(100*yMax[,2]/yMax[,1])
  # Make ROCs
  if (DA==T) {
    ROCMin=roc(Y,yMin[,1])
    ROCMid=roc(Y,yMid[,1])
    ROCMax=roc(Y,yMax[,1])
    # Classify predictions
    yClassMin=ifelse(yMin[,1]>coords(ROCMin,'b',ret='t')[1],1,-1)
    yClassMid=ifelse(yMid[,1]>coords(ROCMid,'b',ret='t')[1],1,-1)
    yClassMax=ifelse(yMax[,1]>coords(ROCMax,'b',ret='t')[1],1,-1)
    # Bind together all Y data
    yMat=cbind(yMin[,1],yMid[,1],yMax[,1],yClassMin,yClassMax,Y)
    colnames(yMat)[1:3]=c('yMin','yMid','yMax')
    # Calculate misclassifications
    missMin=sum(abs((Y-yClassMin))/2)
    missMid=sum(abs((Y-yClassMid))/2)
    missMax=sum(abs((Y-yClassMax))/2)
  }
  cat(X,Y,ID,nRep,nOuter,nInner,varRatio,fitness,method,methParam,ML)
  return(list(X,Y,ID))
}

