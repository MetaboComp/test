#' Support vector machine model in inner CV loop
#'
#'
#'
#' @param xTrain Training data (samples as rows; variables as columns)
#' @param yTrain Training response
#' @param xVal Validation data
#' @param yVal Validation response (for tuning)
#' @param DA Logical for discriminant analysis (classification)
#' @param fitness Fitness function ('MISS', 'AUROC' or 'RMSEP')
#' @param ntree See original function (`randomForest`). Passed from wrapper.
#' @param method Choice of Random Forest implementation (randomForest, ranger or Rborist). Passed from wrapper.
#' @param cost cost
#'
#' @return An object containing:
#' @return (`miss`,`ber`,`auc` or `rmsep`) A fitness metric
#' @return `virank` variable importance rankings
#' @export
#'
#'
#'
#'
#'
#'
#'
#'
#'

svmInner<-function(x,
                   y,
                   epsilon,
                   kernal,
                   cost,
                   gamma,
                   degree){



}

####tune for gamma and cost























