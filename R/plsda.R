#' PLS-DA
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
plsda <- function(X, Y, ncomp = 2, max.iter = 500, tol = 1e-06, near.zero.var = TRUE) {
  Y = as.factor(Y)	

  n <- length(Y)
  groups <- sort(unique(Y))
  levels =  levels(Y)### Add levels
  cgroups <- as.character(groups)
  groups <- as.numeric(factor(cgroups, levels = unique(cgroups)))
  classification <- as.numeric(factor(as.character(Y), levels = unique(cgroups)))
  k <- length(groups)
  nam <- levels(groups)
  ind.mat <- matrix(0, n, k, dimnames = c(names(classification), nam))
  for (j in 1:k) ind.mat[classification == groups[j], j] <- 1
  attr(ind.mat, "levels") = levels
  
  result = MUVR::pls(X, ind.mat, ncomp = ncomp, max.iter = max.iter, tol = tol,near.zero.var=near.zero.var)
  
  result$ind.mat = ind.mat
  result$names$Y = levels(Y)
  return(invisible(result))	
}
