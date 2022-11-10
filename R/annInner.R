#' ANN model in inner CV loop
#'
#' Called from Wrapper
#'
#' @param xTrain Training data (samples as rows; variables as columns)
#' @param yTrain Training response
#' @param xVal Validation data
#' @param yVal Validation response (for tuning)
#' @param DA Logical for discriminant analysis (classification)
#' @param fitness Fitness function ('MISS', 'AUROC' or 'RMSEP')
#' @param layers Number of layers and how many node in each layed. If multiple layers, it writes as c(n,n,n)
#' @param threshold the Threshold of step size
#' @param stepmax Maximum step in each annInner loop
#'
#' @return An object containing:
#' @return (`miss`, `auc` or `rmsep`) A fitness metric
#' @return `virank` variable importance rankings
#' @export
annInner <- function(xTrain,
                    yTrain,
                    xVal,
                    yVal,
                    nodes,
                    threshold,
                    stepmax,
                    DA,
                    fitness,
                    method) {
  # Allocate return object
  returnIn <- list()
  ###Put it in the main MUVR
if(DA==F){

 if(method=="neuralnet"){


 data=cbind(xTrain,yTrain)

  annModIn<-neuralnet(yTrain~.,
                       data=data,
                       threshold=threshold,
                       stepmax=stepmax,
                       hidden=nodes,
                       lifesign="none",  ##a string specifying how much the function will print during the calculation of the neural network. 'none', '
                       # bact.fct="logistic",
                       #act.fct="linear",
                       err.fct="sse",
                       linear.output = T
                       )

   yValInner<-predict(annModIn,xVal)
   returnIn$virank <- rank(-abs(olden(annModIn,bar_plot=F)$importance))
   names(returnIn$virank) <- rownames(olden(annModIn,bar_plot=F))

 }

   else if (method=="nnet"){##annModIn<-nnet()

     annModIn<-nnet(yTrain~.,
                    x=xTrain,
                    y=yTrain,
                    data=data,
                    #decay=threshold,
                    maxit=stepmax,
                    size=nodes,
                    MaxNWts = 100000,
                    #linout=F

     )

   }
  else{
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

      annModIn<-neuralnet(as.formula(
                            paste(paste(names_y,collapse="+"),
                           "~",
                           paste(names_x,collapse="+")
                           )),
                          data=data,
                          threshold=threshold,
                          stepmax=stepmax,
                          hidden=nodes,
                          lifesign="none",
                          act.fct="logistic",  ##sigmoid is logistic

                          err.fct="ce",
                          linear.output = F
      )
      yValInner<-levels(yVal)[apply(predict(annModIn,xVal),
                        1,
                        which.max)]


  olden_frame<-olden(annModIn,
                     x_names=names_x,
                     y_names=colnames(annModIn$response),
                     out_var=colnames(annModIn$response),
                     bar_plot=F)
  returnIn$virank <- rank(-abs(olden_frame$importance))
  names(returnIn$virank) <- rownames(olden_frame)


  #returnIn$virank <- rank(-garson(annModIn,bar_plot=F)$rel_imp)
  #names(returnIn$virank) <- rownames(garson(annModIn,bar_plot=F))

    }

    if (method=="nnet"){##annModIn<-nnet()


    }
    }




  if (fitness == 'RMSEP') {
    returnIn$rmsep <- sqrt(sum((yVal - yValInner) ^ 2, na.rm = T) / (length(yVal) - sum(is.na(yValInner))))
  }
  if (fitness == 'MISS') {
    # cat(' miss',count)
    if (DA) returnIn$miss <- sum(yValInner != yVal)
    else {
      yClassInner <- ifelse(yValInner > 0, 1, -1)   #####classification binary??????Why?????
      returnIn$miss <- sum(yClassInner != yVal)
    }
  }

  if (fitness == 'BER') {
    returnIn$ber <- getBER(actual = yVal,
                           predicted = yValInner)  ###getBER from MUVR
  }
  ##Balance error rate

  if (fitness == 'AUROC') {
    returnIn$auc <- roc(yVal,
                        yValInner)$auc  ####what is votes
  }

  return(returnIn)
}
