#' RF model in inner CV loop 
#'
#' Called from Wrapper
#'
#' @param xTrain Training data (samples as rows; variables as columns)
#' @param yTrain Training response
#' @param xVal Validation data
#' @param yVal Validation response (for tuning)
#' @param DA Logical for discriminant analysis (classification)
#' @param fitness Fitness function ('MISS', 'AUROC' or 'RMSEP')
#' @param ntree See original function (`randomForest`). Passed from wrapper.
#' @param method Choice of Random Forest implementation (randomForest, ranger or Rborist). Passed from wrapper.
#' @param mtry See original function (`randomForest`). Passed from wrapper.
#'
#' @return An object containing:
#' @return (`miss`, `auc` or `rmsep`) A fitness metric 
#' @return `vi` variable importance rankings
#' @export
#'
rfInner <- function(xTrain,
                    yTrain,
                    xVal,
                    yVal,
                    DA,
                    fitness,
                    ntree,
                    mtry, 
                    method) {
  # Allocate return object
  returnIn <- list()
  
  # Different functions for different RF implementations
  if(method == 'randomForest') {
    rfModIn <- randomForest(xTrain,
                            yTrain,
                            xVal,
                            yVal,
                            ntree = ntree,
                            mtry = mtry)
    # Extract predictions
    yValInner <- rfModIn$test$predicted 
    # Variable importance rank
    returnIn$vi <- rank(-rfModIn$importance)
    names(returnIn$vi) <- rownames(rfModIn$importance)
  } else if(method == 'ranger') {
    rfModIn <- ranger(x = xTrain,
                      y = yTrain,
                      num.trees = ntree,
                      importance = 'impurity',
                      # respect.unordered.factors = 'order', # Sort this out as "ifany"
                      mtry = mtry)
    # Extract predictions
    yValInner <- predict(rfModIn, data = xVal)$predictions
    # Variable importance rank
    returnIn$vi <- rank(-rfModIn$variable.importance)
  } else {
    stop('other RF methods not yet implemented')
  }
  if (fitness == 'MISS') {
    # cat(' miss',count)
    if (DA) returnIn$miss <- sum(yValInner != yVal) else {
      yClassInner <- ifelse(yValInner > 0, 1, -1)
      returnIn$miss <- sum(yClassInner != yVal)
    }
  } 
  if (fitness == 'BER') {
    returnIn$ber <- getBER(actual = yVal, predicted = yValInner)
  }
  if (fitness == 'AUROC') {
    returnIn$auc <- roc(yVal,rfModIn$test$votes[,1])$auc
  }
  if (fitness == 'RMSEP') {
    returnIn$rmsep <- sqrt(sum((yVal - yValInner) ^ 2, na.rm = T) / (length(yVal) - sum(is.na(yValInner))))
  }
  return(returnIn)
}
