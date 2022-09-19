#' Make custom parameters for rdcvNet internal modelling
#'
#' Custom parameters can be set in the function call or by manually setting "slots" in the resulting methParam object
#'
#' @param robust Robustness (slack) criterion for determining min and max knees (defaults to 0.05)
#' @param family
#'
#' @return a `methParam` object
#' @export
#'
#' @examples
#' # Standard parameters for rdcvNet
#' methParam <- rdcvNetParams()
rdcvNetParams <- function(robust=0.05, family='gaussian', nRepInner=1) {
  methParam <- list(robust=robust, family=family, nRepInner=nRepInner)
  return(methParam)
}
