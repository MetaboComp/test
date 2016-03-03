#' RF model in inner CV loop 
#'
#' Called from Wrapper
#' @param xTrain Training data (samples as rows; variables as columns)
#' @param yTrain Training response
#' @param xVal Validation data
#' @param yVal Validation response (for tuning)
#' @param fitness Fitness function ('misClass', 'AUROC' or 'RMSEP')
#' @param ntree See original function (`randomForest`). Passed from wrapper.
#' @param mtry See original function (`randomForest`). Passed from wrapper.
#' @return An object containing:
#' @return (`miss`, `auc` or `rmsep`) A fitness metric 
#' @return `vip` VIP rankings
#' @export
#'
rfInner=function(xTrain,yTrain,xVal,yVal,fitness,ntree,mtry) {
  rfModIn=randomForest(xTrain,yTrain,xVal,yVal,ntree=ntree,mtry=mtry)
  yValInner=rfModIn$test$predicted 
  returnIn=list()
  returnIn$vip=rank(-rfModIn$importance)
  names(returnIn$vip)=rownames(rfModIn$importance)
  if (fitness=='misClass') {
    # cat(' miss',count)
    returnIn$miss=sum(rfInner$test$predicted!=yVal)
  } 
  if (fitness=='AUROC') {
    # cat(' auc',count)
    returnIn$auc=auc=roc(yVal,yValInner)$auc
  }
  if (fitness=='RMSEP') {
    # cat(' rmsep',count)
    returnIn$rmsep=sqrt(sum((yVal-yValInner)^2,na.rm=T)/(length(yVal)-sum(is.na(yValInner))))
  }
  return(returnIn)
}