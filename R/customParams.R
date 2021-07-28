#' Make custom parameters for MUVR internal modelling,not rdCV
#'
#' Custom parameters can be set in the function call or by manually setting "slots" in the resulting methParam object
#' Pls note that, at present, there is no mtryMax for the outer (consensus) loop in effect
#'
#' @param method PLS or RF (default)
#' @param robust Robustness (slack) criterion for determining min and max knees (defaults to 0.05)
#' @param ntreeIn RF parameter: Number of trees in inner cross-validation loop models (defaults to 150)
#' @param ntreeOut RF parameter: Number of trees in outer (consensus) cross-validation loop models (defaults to 300)
#' @param mtryMaxIn RF parameter: Max number of variables to sample from at each node in the inner CV loop (defaults to 150). Will be further limited by standard RF rules (see randomForest documentation)
#' @param compMax PLS parameter: Maximum number of PLS components (defaults to 5)
#' @param rfMethod
#' @param oneHot
#' @param NZV
#'
#' @return a `methParam` object
#' @export
#'
#' @examples
#' # Standard parameters for random forest
#' methParam <- customParams() # or
#' methParam <- customParams('RF')
#' # Custom `ntreeOut` parameters for random forest
#' methParam <- customParams('RF',ntreeOut=50) # or
#' methParam <- customParams('RF')
#' methParam$ntreeOut <- 50
#' methParam
customParams <- function(method = c('RF','PLS','randomForest','ranger'),
                         robust = 0.05,
                         ntreeIn = 150,
                         ntreeOut = 300,
                         mtryMaxIn = 150,
                         compMax = 5,
                         oneHot=T,
                         NZV=T,
                         rfMethod=c('randomForest','ranger')
                           ) {


  # Allocate methParam object
  methParam <- list(robust=robust)

  # Random Forest as default method
  ##when method is missing
  if(missing(method)) method <- 'RF'
  ##when method is not missing but is not the method included in the function
  if(!missing(method))
     {if(method!="RF"&method!="PLS"&method !="randomForest"&method!="ranger")
       {stop('other methods not implemented')}}
  ##when rfMethod is not missing but is not the method included in the function
  if(!missing(rfMethod))
     {if(rfMethod !="randomForest"&rfMethod!="ranger")
        {stop('other rfMethods not implemented')}}
  ##when oneHot i snot missing but not T& F
  if(!missing(oneHot))
      {if(oneHot!=T&oneHot!=F)
         {stop('oneHot can only be defined as TRUE or FAlSE')}}
  ##when NZV i snot missing but not T& F
  if(!missing(NZV))
      {if(NZV!=T&NZV!=F)
         {stop('NZV can only be defined as TRUE or FAlSE')}}


# Default oneHot values per method
#    if(missing(oneHot)) {
#      if (method == 'PLS') {
#        oneHot <- TRUE
#      } else if(method == 'RF') {
#        oneHot <- FALSE
#      } }
    methParam$oneHot <- oneHot
# Default oneHot values per method
#    if(missing(NZV)) {
#      if (method == 'PLS') {
#        NZV <- TRUE
#      } else if(method == 'RF') {
#        NZV <- FALSE
#      } }
     methParam$NZV <- NZV


  # Fix shorthand for different RF implementations
  if(method == 'randomForest') {
    if(rfMethod=="ranger"){stop('Contradictory order')
      }else{
    method <- 'RF'
    methParam$method<-'RF'
    methParam$rfMethod <- 'randomForest'
  }}
  if(method == 'ranger') {
    if(rfMethod=="randomForest"){stop('Contradictory order')
      }else{
    method <- 'RF'
    methParam$method<-'RF'
    methParam$rfMethod <- 'ranger'
  }}

if (method=="PLS"&oneHot==F){stop("PLS method must use oneHot encoding. ")}
if (method=="PLS"&NZV==F){stop("PLS method must use near zero variance. ")}

  #########################
  # Default RF parameters
  # For PLS oneHot and NZV can only be true, for RF, the default value is set as T but can be changed to F
  #########################

  if (method=='RF'){
    methParam$ntreeIn <- ntreeIn
    methParam$ntreeOut <- ntreeOut
    methParam$mtryMaxIn <- mtryMaxIn
    if(missing(rfMethod)){methParam$rfMethod <- 'randomForest'}
    else{methParam$rfMethod<-rfMethod}
   }
    # default RF implementation

    else if (method=='PLS') {
    if(!is.null(rfMethod)){stop('Method is PLS. There should not be rfMethod')}
    #########################
    # Default PLS parameters
    #########################

    methParam$compMax <- compMax
  }
  return(methParam)
}
