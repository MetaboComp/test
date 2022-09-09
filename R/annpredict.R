#' ANN predictions for outer segments and consensus models
#'
#' @param xTrain xTrain
#' @param yTrain yTrain
#' @param xTest xTest
#' @param yTest yTest
#' @param method method
#' @param DA DA
#' @param nodes
#' @param stepmax
#' @param threshold
#' @return  The predicted value of yTest
#' @export
annpredict <- function(xTrain,
                   yTrain,
                   xTest,
                   yTest,
                   DA,
                   nodes,
                   threshold,
                   stepmax,
                   #If set to FALSE, the forest will not be retained in the output object. If xtest is given, defaults to FALSE.
                   method) {

  # Allocate return object
  return <- list()

  # Use "Train" for "Testing" if lacking (for fit-predict)
  if(missing(xTest)) {
    xTest <- xTrain
    yTest <- yTrain
  }



  if(DA==F){
  if(method == 'neuralnet') {
    data=cbind(xTrain,yTrain)
    return$model <- neuralnet(yTrain~.,
                              data=data,
                              hidden=nodes,
                              threshold=threshold,
                              stepmax=stepmax,
                              lifesign="minimal",
                              act.fct="logistic",
                              err.fct="sse"
                              #linear.output = F
                              )
    return$fit <-predict(return$model ,xTrain)
    return$predicted<-predict(return$model ,xTest)
    }


   else if(method == 'nnet') {



  } else {
    stop('other ANN methods not yet implemented')
  }
}
if(DA==T){
  if(method=="neuralnet"){
    data=cbind(xTrain,yTrain)
    ytrain<-as.data.frame(yTrain)
    yTrain<-with(ytrain,
                 data.frame(model.matrix(~yTrain+0)))
    data<-cbind(xTrain,yTrain)
    names_x<-colnames(xTrain)
    names_y<-names(yTrain)

    return$model<-neuralnet(as.formula(
      paste(paste(names_y,collapse="+"),
            "~",
            paste(names_x,collapse="+")
      )),
      data=data,
      threshold=threshold,
      stepmax=stepmax,
      hidden=nodes,
      lifesign="none",
      act.fct="logistic",
      err.fct="sse",
      linear.output = F
    )
    return$fit <-predict(return$model ,xTrain)
    return$predicted<-predict(return$model ,xTest)
  }

}


  return(return)
}
