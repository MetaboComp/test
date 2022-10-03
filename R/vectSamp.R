#' Sampling of a vector
#'
#' Random sampling of a vector into n groups
#' Bug: Sampling is faulty when length(vect) < n
#' @param vect A vector to be sampled into groups
#' @param n Number of groups
## @param sampLen A vector with custom number of samples per group. Is calculated if missing (best choice).
#' @return a list with n groups containing sub sampled `vect`
#' @export
vectSamp=function(vect,
                  n=4  #,
                  #sampLen
                  ) { # Pick 'n' random samples within vector 'vect'
  # sampLen is a vector of number of observations within sample
  # If sampLen is not specified it is automatically calculated

  library(caret)
  fold<-caret::createFolds(y=vect, k = n)  ###stratified spliting

  data<-list()
  for(i in 1:n){
    data[[i]]<-vect[fold[[i]]]
  }
  sorteddata<-data[order(sapply(data,length))]
  return(sorteddata)
}




