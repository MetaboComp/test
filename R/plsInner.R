#' PLS model in inner CV loop 
#'
#' Called from Wrapper
#' @param xTrain Training data (samples as rows; variables as columns)
#' @param yTrain Training response
#' @param xVal Validation data
#' @param yVal Validation response (for tuning)
#' @param fitness Fitness function ('misClass', 'AUROC' or 'RMSEP')
#' @param comp Max number of components to try
#' @param mode PLS method (defaults to regression, see 'mixOmics' for details)
#'
#' @return An object containing:
#' @return (`miss`, `auc` or `rmsep`) A fitness metric 
#' @return `nComp` Optimised number of components within range (1:comp)
#' @return `vip` VIP rankings
#' @export
#'
plsInner=function(xTrain,yTrain,xVal,yVal,fitness,comp,mode='regression') {
  plsModIn=pls(xTrain,yTrain,ncomp=comp,mode=mode)
  if (length(plsModIn$nzv$Position)>0) {
    removeVar=rownames(plsModIn$nzv$Metrics)
    xVal=xVal[,!colnames(xVal)%in%removeVar]
  }
  yValInner=predict(plsModIn,newdata=xVal)$predict[,,]  # Store  prediction estimates per validation segment 
  returnIn=list()
  if (fitness=='misClass') {
    # cat(' miss',count)
    yClassInner=ifelse(yValInner>0,1,-1)
    misClass=apply(yClassInner,2,function(x) sum(x/yVal!=1))
    returnIn$miss=min(misClass)
    nComp=which.min(misClass)
  } 
  if (fitness=='AUROC') {
    # cat(' auc',count)
    auc=apply(yValInner,2,function(x) roc(yVal,x)$auc)
    returnIn$auc=max(auc)
    nComp=which.max(auc)
  }
  if (fitness=='RMSEP') {
    # cat(' rmsep',count)
    rmsep=apply(yValInner,2,function(x) sqrt(sum((yVal-x)^2,na.rm=T)/(length(yVal)-sum(is.na(x)))))
    returnIn$rmsep=min(rmsep)
    nComp=which.min(rmsep)
  }
  returnIn$nComp=nComp
  returnIn$vip=rank(-vip(plsModIn)[,nComp])
  return(returnIn)
}