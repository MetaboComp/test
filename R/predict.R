#' Predict pls
#'
#' Adapted and stripped down from mixOmics v 5.2.0 (https://cran.r-project.org/web/packages/mixOmics/)
#'
#' @param object a plsMUVR object
#' @param newdata new data
#' @param onlyPred Boolean for whether to report back predictions only (defaults to FALSE)
#' @param ... 
#'
#' @return pls prediction
#' @export
predict.plsMUVR <-function(object, newdata, onlyPred=FALSE, ...){
  #-- validation des arguments --#
  if (missing(newdata)) stop("No new data available.")
  
  x = object$X
  y = object$Y
  q = ncol(y)
  p = ncol(x)
  
  if (length(dim(newdata)) == 2) {
    if (ncol(newdata) != p)
      stop("'newdata' must be a numeric matrix with ncol = ", p, 
           " or a vector of length = ", p, ".")
  }
  if (length(dim(newdata)) == 0) {
    if (length(newdata) != p)
      stop("'newdata' must be a numeric matrix with ncol = ", p, 
           " or a vector of length = ", p, ".")
    dim(newdata) = c(1, p) 
  }
  
  #-- initialisation des matrices --#	
  ncomp = object$ncomp
  a = object$loadings$X
  b = object$loadings$Y
  c = object$mat.c
  
  means.x = attr(x, "scaled:center")
  means.y = attr(y, "scaled:center")
  sigma.x = attr(x, "scaled:scale")
  sigma.y = attr(y, "scaled:scale")
  
  newdata = as.matrix(newdata)
  ##- coeff de regression 
  B.hat = array(0, dim = c(p, q, ncomp))
  ##- prediction
  y.hat = y.hat2 = array(0, dim = c(nrow(newdata), q, ncomp))
  ##- variates
  t.pred = array(0, dim = c(nrow(newdata), ncomp))
  
  variates.x = object$variates$X
  betay = list()
  
  #-- prediction --#
  for(h in 1:ncomp){
    
    dd= coefficients(lm(y~variates.x[,1:h,drop=FALSE])) #regression of y on variates.global.x => =loadings.global.y at a scale factor
    if (q==1) betay[[h]]=(dd[-1]) else betay[[h]]=(dd[-1,])
    
    W = a[, 1:h,drop=FALSE] %*% solve(t(c[, 1:h,drop=FALSE]) %*% a[, 1:h,drop=FALSE])
    B = W %*% drop(betay[[h]])
    
    y.temp=scale(newdata,center=means.x,scale=sigma.x) %*% as.matrix(B) #so far: gives a prediction of y centered and scaled
    y.temp=scale(y.temp,center=FALSE,scale=1/sigma.y) #so far: gives a prediction of y centered, with the right scaling
    y.temp=scale(y.temp,center=-means.y,scale=FALSE) #so far: gives a prediction of y with the right centering and scaling
    
    y.hat[, , h] = y.temp # we add the variance and the mean of y used in object to predict
    t.pred[, h] = scale(newdata, center = means.x, scale = sigma.x) %*% W[, h]
    B.hat[, , h] = B
  }  #end h
  
  #-- valeurs sortantes --#
  rownames(t.pred) = rownames(newdata)
  colnames(t.pred) = paste("dim", c(1:ncomp), sep = " ")
  rownames(y.hat) = rownames(newdata)
  colnames(y.hat) = colnames(y)
  
  if (onlyPred) return(invisible(list(predict = y.hat))) else return(invisible(list(predict = y.hat, variates = t.pred, B.hat = B.hat,betay=betay)))
}


#' Predict plsda
#'
#' Adapted and stripped down from mixOmics v 5.2.0 (https://cran.r-project.org/web/packages/mixOmics/)
#'
#' @param object a plsdaMUVR object
#' @param newdata new data
#' @param onlyPred Boolean for whether to report back predictions only (defaults to FALSE)
#' @param ... 
#'
#' @return plsda predictions
#' @export
predict.plsdaMUVR=function(object, newdata, onlyPred=FALSE, ...)  {
  #-- validation des arguments --#
  if (missing(newdata)) stop("No new data available.")
  
  x = object$X
  y = object$Y 
  yprim = object$ind.mat   
  q = ncol(yprim)          
  p = ncol(x)
  
  if (length(dim(newdata)) == 2) {
    if (ncol(newdata) != p)
      stop("'newdata' must be a numeric matrix with ncol = ", p, 
           " or a vector of length = ", p, ".")
  }
  if (length(dim(newdata)) == 0) {
    if (length(newdata) != p)
      stop("'newdata' must be a numeric matrix with ncol = ", p, 
           " or a vector of length = ", p, ".")
    dim(newdata) = c(1, p) 
  }
  
  #-- initialisation des matrices --#	
  ncomp = object$ncomp
  a = object$loadings$X
  b = object$loadings$Y
  c = object$mat.c
  
  means.x = attr(x, "scaled:center")
  means.y = attr(y, "scaled:center")
  sigma.x = attr(x, "scaled:scale")
  sigma.y = attr(y, "scaled:scale")
  
  newdata = as.matrix(newdata)
  ##- coeff de regression 
  B.hat = array(0, dim = c(p, q, ncomp))
  ##- prediction
  y.hat = array(0, dim = c(nrow(newdata), q, ncomp))
  ##- variates
  t.pred = array(0, dim = c(nrow(newdata), ncomp))
  variates.x = object$variates$X
  betay = list()
  
  #-- prediction --#
  for(h in 1:ncomp){
    dd= coefficients(lm(y~variates.x[,1:h,drop=FALSE])) #regression of y on variates.global.x => =loadings.global.y at a scale factor
    if(q==1) betay[[h]]=(dd[-1]) else betay[[h]]=(dd[-1,])
    
    W = a[, 1:h,drop=FALSE] %*% solve(t(c[, 1:h,drop=FALSE]) %*% a[, 1:h,drop=FALSE])
    B = W %*% drop(betay[[h]])
    
    y.temp=scale(newdata,center=means.x,scale=sigma.x) %*% as.matrix(B) #so far: gives a prediction of y centered and scaled
    y.temp=scale(y.temp,center=FALSE,scale=1/sigma.y) #so far: gives a prediction of y centered, with the right scaling
    y.temp=scale(y.temp,center=-means.y,scale=FALSE) #so far: gives a prediction of y with the right centering and scaling
    
    y.hat[, , h] = y.temp # we add the variance and the mean of y used in object to predict
    t.pred[, h] = scale(newdata, center = means.x, scale = sigma.x) %*% W[, h]
    B.hat[, , h] = B
  }  #end h
  
  #-- valeurs sortantes --#
  rownames(t.pred) = rownames(newdata)
  colnames(t.pred) = paste("dim", c(1:ncomp), sep = " ")
  rownames(y.hat) = rownames(newdata)
  colnames(y.hat) = colnames(y)
  
  if (onlyPred) return(invisible(list(predict = y.hat))) else return(invisible(list(predict = y.hat, variates = t.pred, B.hat = B.hat,betay=betay)))
}
