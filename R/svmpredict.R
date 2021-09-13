#' RF predictions for outer segments and consensus models
#'
#' @param xTrain
#' @param yTrain
#' @param xTest
#' @param yTest
#' @param
#' @param
#' @param kernelmethod
#' @param DA
#'
#' @return  The predicted value of yTest
svmpredict <- function(xTrain,
                   yTrain,
                   xTest,
                   yTest,
                   scale,
                   method=c("svm","ksvm","svmlight"),
                   kernelmethod=c("linear","radical","polynomial","sigmoid"),
                   type=c("C-classification","nu-classification","one-classification","eps-regression","nu-regression"),
                   ) {


  # Allocate return object
  return <- list()
  if(missing(xTrain)|missing(yTrain)){stop("There must be x and y training data")}
  if(missing(xTest)) {
    xTest <- xTrain
    yTest <- yTrain
  }
  dat<-data.frame(x=xTrain,
                  y=yTrain)                           #####one thing to notice is that when using svm, y must be a factor
  # Use "Train" for "Testing" if lacking (for fit-predict)

  if(class(yTest)!=class(yTrain)){stop("yTrain must be the same class as yTest")
  }
  tesdat<-data.frame(x=xTest,
                     y=yTest)
  if(missing(method)){method="svm"}
  if(!method %in% c("svm","ksvm","svmlight")){stop("This method is not applied")}
  ##some notes
  ##for svm, when y is factor, it automatically use C-classification, when y is numeric, it use eps regression

  ##C-classification: the range  is for range 0 ,nu-classification is for range(0,1)
  ##svm has cross validation function build inside,defalut is cross=0, To use that one need to specified cross= in svm()


#####method 1
if(method=="svm"){
library(e1071)
  if (missing(kernelmethod)) { kernelmethod = "radial"}
  if (!kernelmethod %in% c("linear", "radial", "polynomial", "sigmoid"))
  {stop("This method is not included in svm")}

  if (missing(type)) {
    if (class(Y) == "numeric"){type = "C-Classification"}
    if (class(Y) == "numeric" | class(Y) == "integer"){type = "eps-regression"}
  }

  if (!type %in% c("C-classification","nu-classification","one-classification","eps-regression","nu-regression" ))
  { stop("This type is not included in svm")}

  if (type %in% c("C-classification","nu-classification","one-classification")) {
    if (class(Y) != "factor") { stop("dependent variable has to be of factor type for classification mode.")}
  }
  if (type %in% c("C-classification", "nu-classification","one-classification")) {
    if (class(Y) != "numeric" | class(Y) != "integer") { stop("Need numeric dependent variable for regression")}
  }






 suppressWarnings(
  if(kernelmethod=="radial")
   {tunesvm=tune(svm,
                   y~.,
                   data=dat,
                  kernel=kernelmethod,
                   type=type,
                   ranges=list(cost=c(0.0001,0.001,0.01,0.1,1,10,100,1000),
                               gamma=c(0.000001,0.00001,0.0001,0.001,0.01,0.1,1)))

 bestcost<-tunesvm$best.model$cost
 bestgamma<-tunesvm$best.model$gamma ###gamma is only used when kernel method i radial
}   else{tunesvm=tune(svm,
                      y~.,
                      data=dat,
                       kernel=kernelmethod,
                      type=type,
                      ranges=list(cost=c(0.0001,0.001,0.01,0.1,1,10,100,1000))
                     )
bestcost<-tunesvm$best.model$cost
       }
 )
###Here there are 2 methods which should be used.
#(1) use directly best model from tuning
#(2) use parameter selected from tunning to run a new model
 #####################################################
  if(kernelmethod=="radial")
  {return$model <- svm(y~.,
                      data=dat,
                      kernel=kernelmethod,
                      type=type,
                      cost=bestcost,
                      gamma=bestgamma
                      )
  }else {return$model <- svm(y~.,
                            data=dat,
                            kernel=kernelmethod,
                            type=type,
                            cost=bestcost
                            )
  }
################################################

return$model<-tunesvm$best.model

  return$fit<-fitted(return$model)

  return$predicted<-predict(return$model,testdat)
  table(predict=predicted,truth=testdat$y)
}

####for ksvm there are 9 types in total. For simplicity, I will just include the method that is included in svm
####method 2

if(method=="ksvm")
{library(kernlab)
  if(type)
ksvm(y~.,
     data=dat,
     type)

}

  if(method=="svmlight")
  {library(klaR)


  }




  return(return)
}
