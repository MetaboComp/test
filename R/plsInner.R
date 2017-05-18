#' PLS model in inner CV loop 
#'
#' Called from Wrapper
#' @param xTrain Training data (samples as rows; variables as columns)
#' @param yTrain Training response
#' @param xVal Validation data
#' @param yVal Validation response (for tuning)
#' @param DA Logical for discriminant analysis (classification)
#' @param fitness Fitness function ('MISS', 'AUROC' or 'RMSEP')
#' @param comp Max number of components to try
#' @param mode PLS method (defaults to regression, see 'mixOmics' for details)
#'
#' @return An object containing:
#' @return (`miss`, `auc` or `rmsep`) A fitness metric 
#' @return `nComp` Optimised number of components within range (1:comp)
#' @return `vip` VIP rankings
#' @export
#'
plsInner=function(xTrain,yTrain,xVal,yVal,DA,fitness,comp,mode='regression') {
  cond=TRUE
  while(cond) {
    if (DA) plsModIn=plsda(xTrain,yTrain,ncomp=comp,near.zero.var=TRUE) else 
      plsModIn=pls(xTrain,yTrain,ncomp=comp,mode=mode,near.zero.var=TRUE)
    yValInner=tryCatch(predict(plsModIn,newdata=xVal)$predict[,,,drop=F], error=function(e) return('error'))
    if (any(yValInner=='error') | any(is.na(yValInner))) comp=comp-1 else cond=FALSE
  }
  returnIn=list()
  if (DA) {
    if (fitness=='MISS') {
      classes=apply(yValInner,c(1,3),which.max)
      misClass=apply(classes,2,function(x) sum(x!=as.numeric(yVal)))
      returnIn$miss=min(misClass,na.rm=T)
      nComp=which.min(misClass)
    } else {
      auc=apply(yValInner[,1,],2,function(x) roc(yVal,x)$auc)
      returnIn$auc=max(auc,na.rm=T)
      nComp=which.max(auc)
    }
  } else {
    if (fitness=='MISS') {
      # cat(' miss',count)
      yClassInner=ifelse(yValInner>0,1,-1)
      misClass=apply(yClassInner,2,function(x) sum(x!=yVal))
      returnIn$miss=min(misClass,na.rm=T)
      nComp=which.min(misClass)
    } 
    if (fitness=='AUROC') {
      # cat(' auc',count)
      auc=apply(yValInner,2,function(x) roc(yVal,x)$auc)
      returnIn$auc=max(auc,na.rm=T)
      nComp=which.max(auc)
    }
    if (fitness=='RMSEP') {
      # cat(' rmsep',count)
      rmsep=apply(yValInner,2,function(x) sqrt(sum((yVal-x)^2,na.rm=T)/(length(yVal)-sum(is.na(x)))))
      returnIn$rmsep=min(rmsep,na.rm=T)
      nComp=which.min(rmsep)
    }
  }
  returnIn$nComp=nComp
  returnIn$vip=rank(-vip(plsModIn)[,nComp])
  return(returnIn)
}
