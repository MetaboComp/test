#' Predict MV object
#'
#' At present, this function only supports predictions for PLS regression type problems
#'
#' @param MVObj An `MVobject` obtained from the MVWrap function
#' @param newdata New data for which to predict outcomes
#' @param model What type of model to plot ('min', 'mid' or 'max'). Defaults to 'mid'.
#'
#' @return A pdf with plots of results from multivariate predictions
#' @export
predMV <- function(MVObj,newdata,model='mid') {
  
  # Check for correct class of object
  if (!any(class(MVObj) == 'MVObject')) {
    cat('\nWrong object class: Return NULL')
    return(NULL)
  }
  
  # Extract model parameters
  modNum <- ifelse(model == 'min',
                   1,
                   ifelse(model == 'mid',
                          2,
                          3))
  method <- MVObj$inData$method
  nRep <- MVObj$inData$nRep
  nOuter <- MVObj$inData$nOuter
  
  if (method == 'PLS') {
    nComps <- MVObj$nCompPerSeg[[modNum]]
  } else {
    library(randomForest)
  }
  
  ###############
  # Regression
  ###############
  
  if (class(MVObj)[2] == 'Regression') {
    
    # Allocate object containing predictions per segment model 
    yPredPerMod <- matrix(ncol = length(MVObj$outModels),
                          nrow = nrow(newdata),
                          dimnames = list(paste('observation', 1:nrow(newdata), sep = ''),
                                          paste('model', 1:length(MVObj$outModels), sep = '')))
    
    # Make predictions for each segment
    n <- 0 # Set segment model counter to zero
    
    for(r in 1:nRep) {
      
      for(i in 1:nOuter) {
        
        n <- n + 1 # Increase segment model counter
        
        mod <- MVObj$outModels[[n]][[modNum]] # Extract model
        
        if (method == 'PLS') {
          
          # Sanity check predictor variables
          if (any(!colnames(mod$X) %in% colnames(newdata))) {
            cat('\nMismatch variable names in model', n, ': Return NULL')
            return(NULL)
            
          } else { # If ok then run
            
            # Extract model-specific params
            X <- subset(newdata, select = colnames(mod$X)) # Subset variables to those in the actual model
            nComp <- nComps[r, i] 
            
            # Make actual model predictions
            yPredPerMod[, n] <- predict(mod, newdata = X)$predict[, , nComp]  # 
            
          }
          
        } else { # If Random Forest
          
          # Sanity check predictor variables
          if (any(!rownames(mod$importance) %in% colnames(newdata))) {
            cat('\nMismatch variable names in model', n, ': Return NULL')
            return(NULL)
            
          } else { # If ok then run
            
            # Extract model-specific params
            X <- subset(newdata, select = mod$names$X)
            X <- subset(newdata, select = rownames(mod$importance))
            
            # Make actual model predictions
            yPredPerMod[, n] <- predict(mod, newdata=X)  # 
            
          }
          
        }
        
      }
      
    }
    
    # Average predictions over all models
    yPred <- apply(yPredPerMod,1,mean)
    
    # Return
    return(list(yPred = yPred,
                yPredPerMod = yPredPerMod))
    
    
  } else if (class(MVObj)[2]=='Classification') {

    ###################
    # Classification
    ###################
    
    # Allocate object containing predictions per segment model 
    yPredPerMod <- array(dim = c(nrow(newdata), 
                                 length(levels(MVObj$inData$Y)), 
                                 length(MVObj$outModels)),
                         dimnames = list(paste('observation', 1:nrow(newdata), sep=''), 
                                         levels(MVObj$inData$Y), 
                                         paste('model', 1:length(MVObj$outModels), sep='')))
    
    # Make predictions for each segment
    n <- 0 # Set segment model counter to zero
    
    for(r in 1:nRep) {
      
      for(i in 1:nOuter) {
        
        n <- n + 1 # Increase segment model counter
        
        mod <- MVObj$outModels[[n]][[modNum]] # Extract model
        
        if (method == 'PLS') {
          cat('\nNot yet implemented')
          return(NULL)
          
        } else {
          # RF
          # Sanity check predictor variables
          if (any(!rownames(mod$importance) %in% colnames(newdata))) {
            cat('\nMismatch variable names in model',n,': Return NULL')
            return(NULL)
            
          } else { # If ok then run
            
            # Extract model-specific params
            X <- subset(newdata, select = rownames(mod$importance))
            
            # Make actual model predictions
            yPredPerMod[,,n] <- predict(mod, newdata = X, type='vote')  # 
            
          }
          
        }
        
      }
      
    }
    
    # Average predictions over all models
    yPred <- apply(yPredPerMod, c(1, 2), mean)
    
    # Convert numeric predictions to class
    yClass <- levels(MVObj$inData$Y)[apply(yPred, 1, which.max)]
    names(yClass) <- paste('observation', 1:nrow(newdata), sep = '')
    
    # Return
    return(list(yClass = yClass,
                yPred = yPred,
                yPredPerMod = yPredPerMod))
    
  } else {
    
    ###################
    # Multilevel
    ###################
    
    cat('\nNot yet implemented')
    
  }
  
}
