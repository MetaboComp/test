#' Perform matrix pre-processing
#'
#' @param X Data matrix with samples in rows and variables in columns
#' @param offset Add offset to all data points (defaults to 0)
#' @param zeroOffset Add offset to zero data (defaults to 0)
#' @param trans If specified, accomodates 'log' or 'sqrt' transforms (default is no transformation)
#' @param center Either a logical value or a numeric vector of length equal to the number of columns of X (defaults to TRUE for mean centering).
#' @param scale Either 'UV', 'Pareto', 'none' or a numeric vector of length equal to the number of columns of X (defaults to scaling to unit variance 'UV').
#'
#' @return A pre-processed data matrix
#' @export
preProcess=function(X,offset=0,zeroOffset=0,trans,center,scale) {
  nVar=nrow(X)
  # Add offset
  X=X+offset
  cat('Offset by',offset)
  X[X==0]=zeroOffset
  cat('\nZero offset by',zeroOffset)
  # Perform transformation
  if(!missing(trans)) {
    trans=match.arg(trans,c('log','sqrt'))
    cat('\nTransformation:',trans)
    if(trans=='log') {
      if(any(X<=0)) stop('no zero or negative values allowed when doing log transformation')
      X=apply(X,2,log) 
    }
    if(trans=='sqrt') {
      if(any(X<0)) stop('no negative values allowed when doing sqrt transformation')
      X=apply(X,2,sqrt) 
    }
  }
  if (missing(center)) center=TRUE 
  if(!(is.logical(center) | length(center)==nVar)) stop('center should be either a logical value or a numeric vector of length equal to the number of columns of x')
  if(is.logical(center)) cat('\nCenter:',center) else cat('\nCenter: By vector')
  if (missing(scale)) scale='UV' 
  if (length(scale)!=nVar) {
    scale=match.arg(scale,c('UV','Pareto','none'))
    cat('\nScale:',scale)
    if (scale=='UV') scale=TRUE else if (scale=='none') scale=FALSE else scale=apply(X,2,function(x) sqrt(sd(x)))
  } else cat('\nScale: By vector')
  if(!(is.logical(scale) | length(scale)==nVar)) stop('Error with scaling')
  X=scale(X,center,scale)
}
