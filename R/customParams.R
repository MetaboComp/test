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
#' @param oneHot  T or F using onehot endcoding or not
#' @param scoring_matrix  A matrix that its score can be adjusted
#' @param NZV T or F using NZV or not
#' @param rfMethod randomforest method, which includes randomForest and ranger
#' @param svmMethod support vector machine method which includes 3 different support vector machine methods
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
customParams <- function(method = c('RF','PLS','randomForest','ranger',"SVM"),
                         robust = 0.05,
                         ntreeIn = 150,
                         ntreeOut = 300,
                         mtryMaxIn = 150,
                         compMax = 5,
                         scoring_matrix=F,
                         oneHot,
                         NZV,
                         rfMethod=c('randomForest','ranger'),
                         svmMethod=c("svm","ksvm","svmlight"),
                         annMethod=c("nnet","neuralnet")
                           ) {


  # Allocate methParam object
  methParam <- list(robust=robust)

  # Random Forest as default method
  ##when method is missing
  if(missing(method)) {
    if(!missing(rfMethod)&!missing(svmMethod)){stop("There could only be one method")}
    if(!missing(rfMethod))
      {if(rfMethod !="randomForest"&rfMethod!="ranger")
                                    {stop('other rfMethods not implemented')}
                           methParam$method<-"RF"
                           methParam$ntreeIn <- ntreeIn
                           methParam$ntreeOut <- ntreeOut
                           methParam$mtryMaxIn <- mtryMaxIn
                           methParam$rfMethod<-rfMethod

    }else if(!missing(svmMethod)){
      if(svmMethod !="svm"&svmMethod!="ksvm"&svmMethod!="svmlight")
      {stop('other svmMethods not implemented')}
      methParam$method<-"SVM"
      methParam$svmMethod<-svmMethod
 #     methParam$
#      methParam$
#        methParam$
    ####Here new method could be added
    }else{methParam$method <- 'RF'                            ####when method is missing rfmethod and svmmethod are all missing
        methParam$ntreeIn <- ntreeIn
       methParam$ntreeOut <- ntreeOut
        methParam$mtryMaxIn <- mtryMaxIn
         methParam$rfMethod<-"randomForest"}
  }



  ##when method is not missing but is not the method included in the function
  if(!missing(method))
     {if(method!="RF"&method!="PLS"&method !="randomForest"&method!="ranger"&method!="SVM")
       {stop('other methods not implemented')}
      methParam$method<-method
      if(method=="RF")
        {if(!missing(svmMethod)){stop('Method is RF. There should not be svmMethod')}
        if(missing(rfMethod))
          {methParam$rfMethod<-"randomForest"
        methParam$ntreeIn <- ntreeIn
        methParam$ntreeOut <- ntreeOut
        methParam$mtryMaxIn <- mtryMaxIn
        }else{
          if(rfMethod !="randomForest"&rfMethod!="ranger")
          {stop('other rfMethods not implemented')}
        methParam$rfMethod<-rfMethod
        methParam$ntreeIn <- ntreeIn
        methParam$ntreeOut <- ntreeOut
        methParam$mtryMaxIn <- mtryMaxIn}
      }

      if(method=="randomForest")
      {if(!missing(svmMethod)){stop('Method is randomForest. There should not be svmMethod')}
       if(missing(rfMethod))
         {methParam$rfMethod<-"randomForest"
       } else if(rfMethod !="randomForest"&rfMethod!="ranger")
       {stop('other rfMethods not implemented')
       } else if (rfMethod=="randomForest")
       {methParam$rfMethod<-rfMethod
       } else if(rfMethod=="ranger"){
         stop("could only choose one method")
         }
       methParam$ntreeIn <- ntreeIn
       methParam$ntreeOut <- ntreeOut
       methParam$mtryMaxIn <- mtryMaxIn}

      if(method=="ranger")
      {library(ranger)
        if(!missing(svmMethod)){stop('Method is ranger. There should not be svmMethod')}
        if(missing(rfMethod))
        {methParam$rfMethod<-"ranger"
        }else if(rfMethod !="randomForest"&rfMethod!="ranger")
        {stop('other rfMethods not implemented')
        }else if (rfMethod=="ranger")
        {methParam$rfMethod<-rfMethod
        }else if(rfMethod=="randomForest"){
          stop("could only choose one method")
        }
        methParam$ntreeIn <- ntreeIn
        methParam$ntreeOut <- ntreeOut
        methParam$mtryMaxIn <- mtryMaxIn
        }

      if(method=="PLS"){
        if(!missing(svmMethod)){stop('Method is PLS. There should not be svmMethod')}
        if(!missing(rfMethod)){stop('Method is PLS. There should not be rfMethod')}
        #########################
        # Default PLS parameters
        #########################

        methParam$compMax <- compMax
      }


        if(method=="SVM")
        {if(!missing(rfMethod)){stop('Method is SVM. There should not be rfMethod')}
         if(!missing(svmMethod)){
           if(svmMethod!="svm"&svmMethod!="svm"&svmMethod!="svm"){
             stop("This svmMethod not implemented")}
           }
          if(missing(svmMethod)){
            methParam$svmMethod="svm"
  #          methParam$
 #           methParam$
          }else{

          methParam$svmMethod<-svmMethod
#          methParam$
#          methParam$
          }
          }




         }





  ##when oneHot i snot missing but not T& F
  if(!missing(oneHot))
      {if(oneHot!=T&oneHot!=F)
         {stop('oneHot can only be defined as TRUE or FAlSE')}}
  ##when NZV i snot missing but not T& F
  if(!missing(NZV))
      {if(NZV!=T&NZV!=F)
         {stop('NZV can only be defined as TRUE or FAlSE')}}


# Default oneHot values per method
 if(missing(oneHot)) {
    if (methParam$method == 'PLS') {
       oneHot <- TRUE
     } else if(methParam$method == 'RF'||methParam$method=="randomForest"||methParam$method=="ranger") {
       oneHot <- FALSE
     } else if (methParam$method=="SVM")
     {oneHot <- FALSE}
   }
    methParam$oneHot <- oneHot
# Default oneHot values per method
   if(missing(NZV)) {
      if (methParam$method == 'PLS') {
        NZV <- TRUE
      } else if(methParam$method == 'RF'||methParam$method=="randomForest"||methParam$method=="ranger") {
        NZV <- FALSE
      } else if (methParam$method=="SVM")
      {NZV <- FALSE}
     }
     methParam$NZV <- NZV



if (methParam$method=="PLS"&oneHot==F){stop("PLS method must use oneHot encoding. ")}
if (methParam$method=="PLS"&NZV==F){stop("PLS method must use near zero variance. ")}

  #########################
  # Default RF parameters
  # For PLS oneHot and NZV can only be true, for RF, the default value is set as T but can be changed to F
  #########################



  return(methParam)
}
