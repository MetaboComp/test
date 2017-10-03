#' PLS regression 
#' 
#' Adapted and stripped down from mixOmics v 5.2.0 (https://cran.r-project.org/web/packages/mixOmics/)
#' @param X 
#' @param Y 
#' @param ncomp 
#' @param max.iter 
#' @param tol 
#' @param near.zero.var 
#'
#' @return
#' @export
pls <- function(X, Y, ncomp = 2, max.iter = 500, tol = 1e-06, near.zero.var = TRUE) {
  X = as.matrix(X)
  Y = as.matrix(Y)
  n = nrow(X)
  p = ncol(X)
  q = ncol(Y)
  
  # Remove variables with near zero variance 
  if(near.zero.var) {
    nzv = MUVR::nearZeroVar(X)
    if (length(nzv$Position > 0)) {
      warning("Zero- or near-zero variance predictors.\nReset predictors matrix to not near-zero variance predictors.\nSee $nzv for problematic predictors.")
      X = X[, -nzv$Position,drop=FALSE]
      if(ncol(X)==0) stop("No more predictors after Near Zero Var has been applied!")
    }
  }
  
  # Names
  X.names = colnames(X)
  if (is.null(X.names)) {
    X.names = paste("X", 1:p, sep = "")
    colnames(X)=X.names
  }
  if (dim(Y)[2] == 1) {
    Y.names = "Y"
  } else {
    Y.names = colnames(Y)
    if (is.null(Y.names)) {
      Y.names = paste("Y", 1:q, sep = "")
      colnames(Y)=Y.names
    }
  }
  ind.names = rownames(X)
  if (is.null(ind.names)) {
    ind.names = rownames(Y)
    rownames(X) = ind.names
  }
  if (is.null(ind.names)) {
    ind.names = 1:n
    rownames(X) = rownames(Y) = ind.names
  }	
  
  # Center and scale indata
  X.temp = X = scale(X, center = TRUE, scale = TRUE)
  Y.temp = Y = scale(Y, center = TRUE, scale = TRUE) 
  
  # Allocate matrices
  mat.t = mat.u = matrix(nrow = n, ncol = ncomp)
  mat.a = mat.c = matrix(nrow = p, ncol = ncomp)
  mat.b = mat.d = matrix(nrow = q, ncol = ncomp)
  
  # Iterate pls components h
  iter=NULL
  for (h in 1:ncomp) {
    #-- initialisation --#
    M = crossprod(X.temp, Y.temp)
    svd.M = svd(M, nu = 1, nv = 1)
    a.old = svd.M$u
    b.old = svd.M$v
    #-- latent variables --#
    t = X.temp %*% a.old / drop(crossprod(a.old))
    u = Y.temp %*% b.old / drop(crossprod(b.old))
    iterh = 1
    #-- convergence of a  --#
    repeat{
      a = t(X.temp) %*% u 
      a = a / drop(sqrt(crossprod(a)))
      t = X.temp %*% a / drop(crossprod(a))
      b = t(Y.temp) %*% t 
      b = b / drop(sqrt(crossprod(b)))
      u = Y.temp %*% b / drop(crossprod(b))
      if (crossprod(a - a.old) < tol) break
      if (iterh == max.iter) break
      a.old = a
      b.old = b
      iterh = iterh + 1
    }
    #-- deflation --#
    c = crossprod(X.temp, t) / drop(crossprod(t))
    X.temp = X.temp - t %*% t(c)   
    #-- mode regression --#
    d = crossprod(Y.temp, t) / drop(crossprod(t))
    Y.temp = Y.temp - t %*% t(d)
    
    mat.t[, h] = t
    mat.u[, h] = u
    mat.a[, h] = a
    mat.b[, h] = b
    mat.c[, h] = c
    mat.d[, h] = d
    iter=c(iter,iterh) #save the number of iteration per component
  } 
  #-- valeurs sortantes --#
  rownames(mat.a) = rownames(mat.c) = X.names
  rownames(mat.b) = Y.names
  rownames(mat.t) = rownames(mat.u) = ind.names
  comp = paste("comp", 1:ncomp)
  colnames(mat.t) = colnames(mat.u) = comp
  colnames(mat.a) = colnames(mat.b) = colnames(mat.c) = comp 
  result = list(X = X, Y = Y, ncomp = ncomp, mat.c = mat.c, variates = list(X = mat.t, Y = mat.u), 
                loadings = list(X = mat.a, Y = mat.b), tol = tol, max.iter = max.iter, iter=iter)
  if (near.zero.var == TRUE) result$nzv = nzv
  class(result) = "pls"
  return(invisible(result))
}

