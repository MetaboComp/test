#' Make a confusion matrix from a MUVR object
#'
#' @param MUVRclassObject A MUVR class object (classification analysis)
#' @param model min, mid or max model
#'
#' @return A confusion matrix of actual vs predicted class
#' @export
#'
confusionMatrix=function(MUVRclassObject,model='mid'){
  if(!any(class(MUVRclassObject)=='Classification')) stop ('The MUVR object needs to be from a classification analysis')
  nMod=ifelse(model=='min',1,ifelse(model=='mid',2,3))
  actual=MUVRclassObject$inData$Y
  predicted=MUVRclassObject$yClass[,nMod]
  table(actual=actual,predicted=predicted)
}
