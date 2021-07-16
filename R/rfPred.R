#' RF predictions for outer segments and consensus models
#'
#' @param xTrain 
#' @param yTrain 
#' @param xTest 
#' @param yTest 
#' @param ntree 
#' @param keep.forest 
#' @param method 
#' @param DA 
#'
#' @return
#' @export
rfPred <- function(xTrain,
                   yTrain,
                   xTest,
                   yTest,
                   ntree = 500,
                   DA,
                   keep.forest = FALSE,
                   method) {
  
  # Allocate return object
  return <- list()
  
  # Use "Train" for "Testing" if lacking (for fit-predict)
  if(missing(xTest)) {
    xTest <- xTrain
    yTest <- yTrain
  }
  
  if(method == 'randomForest') {
    return$model <- randomForest(x = xTrain,
                                 y = yTrain,
                                 xtest = xTest,
                                 ytest = yTest,
                                 ntree = ntree,
                                 keep.forest = keep.forest)
    if(DA) {
      return$fit <- return$model$votes
      return$predicted <- return$model$test$votes 
    } else {
      return$fit <- return$model$predicted
      return$predicted <- return$model$test$predicted
    }
  } else if(method == 'ranger') {
    probability <- ifelse(DA, TRUE, FALSE)
    return$model <- ranger(x = xTrain,
                           y = yTrain,
                           num.trees = ntree,
                           importance = 'impurity',
                           probability = probability,
                           # respect.unordered.factors = 'order', # Sort out criterion for "if there are any"
                           write.forest = keep.forest)
    # Extract predictions
    return$fit <- return$model$predictions
    if (keep.forest) return$predicted <- predict(return$model, data = xTest)$predictions
  } else {
    stop('other RF methods not yet implemented')
  }
  return(return)
}