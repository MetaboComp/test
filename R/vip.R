#' Extract PLS VIP values
#'
#' Adapted and stripped down from mixOmics v 5.2.0 (https://cran.r-project.org/web/packages/mixOmics/)
#' @param object pls(da) object
#'
#' @return vip object
#' @export
vip <- function(object) {
  #-- initialisation des matrices --#
  W = object$loadings$X
  H = object$ncomp
  q = ncol(object$Y)
  p = ncol(object$X)
  
  cor2 = cor(object$Y, object$variates$X, use = "pairwise")^2
  cor2 = as.matrix(cor2, nrow = q)
  
  VIP = matrix(0, nrow = p, ncol = H)
  VIP[, 1] = W[, 1]^2
  if (H > 1) {
    for (h in 2:H) {
      if (q == 1) {
        Rd = cor2[, 1:h] 
      } else {
        Rd = colSums(cor2[, 1:h])
      }
      VIP[, h] = Rd %*% t(W[, 1:h]^2) / sum(Rd)
    }
  }
  VIP = sqrt(p * VIP)
  rownames(VIP) = rownames(W)
  colnames(VIP)= paste("comp", 1:H)
  return(invisible(VIP))
}
