#' MVWrap: Wrapper for Multivariate modelling with Variable selection
#' 
#' Repeated double cross validation with tuning of variables in the inner loop.
#' @param X Independent variables (NAs not allowed)
#' @param Y Response vector (Dependent variable)
#' @param ID Sample identifier (to ensure same sample is always in same segment)
#' @param nRep Number of repetitions of double CV..
#' @param nOuter Number of outer CV loop segments.
#' @param nInner Number of inner CV loop segments.
#' @param featRatio Ratio of features to include in subsequent inner loop iteration.
#' @param DA Logical for Classification (discriminant analysis) (Defaults do FALSE, i.e. regression)
#' @param fitness Fitness function for model tuning (choose either 'AUROC' or 'misClass' for classification; or 'RMSEP' (default) for regression.)
#' @param method Multivariate method. Supports 'PLS' and 'RF' (default)
#' @param methParam List with parameter settings for specified MV method (defaults to ???)
#' @param ML Logical for multilevel analysis (defaults to FALSE)
#' @return An object containing stuff...
#' @export
MVWrap=function(X,Y,ID,nRep=5,nOuter=6,nInner,featRatio=0.75,DA=FALSE,fitness=c('AUROC','misClass','RMSEP'),method=c('PLS','RF'),methParam,ML=FALSE){
  start.time=proc.time()[3]
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
      methParam=list(nComp=3)
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
      Yname=unique(Y)
    } else {
      fitness='RMSEP'
    }
  }
  if (nrow(X)!=length(Y)) {
    cat('\nError: Must have same nSamp in X and Y.\n')
    return(NULL)
  }
  # Allocate object to return from function
  modelReturn=list()
  
  cat(X,Y,ID,nRep,nOuter,nInner,featRatio,fitness,method,methParam,ML)
  return(list(X,Y,ID))
}

X=matrix(1:6,ncol=2)
Y=1:3
MVWrap(X,Y,ML=T)
MVWrap(X,Y)
