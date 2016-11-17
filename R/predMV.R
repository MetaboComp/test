#' Predict MV object
#'
#' At present, this function only supports predictions for PLS regression type problems
#'
#' @param MVObj An `MVobject` obtained from the MVWrap function
#' @param newdata New data for which to predict outcomes
#' @param model What type of model to plot ('min', 'mid' or 'max'). Defaults to 'mid'.
#'
#' @return A pdf with plots of results from multivariate predictions
#' @export
predMV=function(MVObj,newdata,model='mid') {
  if (!any(class(MVObj)=='MVObject')) {
    cat('\nWrong object class: Return NULL')
    return(NULL)
  }
  modNum=ifelse(model=='min',1,ifelse(model=='mid',2,3))
  method=MVObj$inData$method
  nRep=MVObj$inData$nRep
  nOuter=MVObj$inData$nOuter
  yPred=matrix(nrow=length(MVObj$outModels),ncol=nrow(newdata),dimnames=list(paste('model',1:length(MVObj$outModels),sep=''),paste('observation',1:nrow(newdata),sep='')))
  # par(mar=c(4,4,0,0)+.5)
  if (class(MVObj)[3]=='Regression') {
    if (method=='PLS') {
      library(mixOmics)
      nComps=MVObj$nCompPerSeg[[modNum]]
      n=0
      for(r in 1:nRep) {
        for(i in 1:nOuter) {
          n=n+1
          mod=MVObj$outModels[[n]][[modNum]]
          X=subset(newdata,select=mod$names$X)
          nComp=nComps[r,i]
          yPred[n,]=predict(mod,newdata=X)$predict[,,nComp]  # 
        }
      }
    } else {
      cat('\nNot yet implemented')
    }
  } else if (class(MVObj)[3]=='Classification') {
    cat('\nNot yet implemented')
  } else {
    cat('\nNot yet implemented')
  }
  return(yPred)
}
