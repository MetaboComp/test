#' Predict MV object using a MUVR class object and a X testing set
#'
#' At present, this function only supports predictions for PLS regression type problems
#'
#' @param MUVRclassobject An 'MUVR' class object obtained from the MVWrap function
#' @param newdata New data for which to predict outcomes
#' @param model What type of model to plot ('min', 'mid' or 'max'). Defaults to 'mid'.
#'
#' @return A pdf with plots of results from multivariate predictions
#' @export
predMV=function(MUVRclassobject,
                newdata,
                model='mid') {
  if (!(class(MUVRclassobject)=='MUVR')) {
    cat('\nWrong object class: Return NULL')
    return(NULL)
  }
  modNum=ifelse(model=='min',1,
                ifelse(model=='mid',
                       2,
                       3))

  method=MUVRclassobject$inData$method
  nRep=MUVRclassobject$inData$nRep
  nOuter=MUVRclassobject$inData$nOuter

 #########################################
  if (method=='PLS') {
    nComps=MUVRclassobject$nCompPerSeg[[modNum]]  ###row is repetition, column is outer segment
  } else {
    library(randomForest)
  }
  #par(mar=c(4,4,0,0)+.5)
  #####################################
####when regression is used
  if (class(MUVRclassobject)[3]=='Regression') {
    yPredPerMod=matrix(ncol=length(MUVRclassobject$outModels),   ###when modReturn is set is true, the length is nRep*nOuter
                       nrow=nrow(newdata),   ##number of observations
                       dimnames=list(paste('observation',1:nrow(newdata),sep=''),
                                     paste('model',1:length(MUVRclassobject$outModels),sep='')))
    n=0
    for(r in 1:nRep) {
      for(i in 1:nOuter) {
        n=n+1
        mod=MUVRclassobject$outModels[[n]][[modNum]]   ###
###when method is PLS
        if (method=='PLS') {
          if ((!colnames(mod$X)%in%colnames(newdata))) {  ###check if training set has variables that is not included in the tesing set
            cat('\nMismatch variable names in model',n,': Return NULL')
            return(NULL)
          } else {

             X=subset(newdata,select=colnames(mod$X))   ###keep testing variables only the variables existing in training set
            nComp=nComps[r,i]   ##r is repetition,i is number of outer segment
            ##yPrePerMod is a matrix of (observation is row, nRep*nOuter)
            yPredPerMod[,n]=predict(mod,newdata=X)$predict[,,nComp]
            ###How this predict function work????
          }
        }

  ###when method is RF
        else {
          if ((!rownames(mod$importance)%in%colnames(newdata))) {  ###variables of higher importance must be included in the testing set
            cat('\nMismatch variable names in model',n,': Return NULL')
            return(NULL)
          } else {
            X=subset(newdata,select=mod$names$X)    ####keep testing variables only the variables existing in training set
            X=subset(newdata,select=rownames(mod$importance))  ###keep testing variables only variable of importance in the training set
            yPredPerMod[,n]=predict(mod,newdata=X)  #
          }
        }
      }
    }
    yPred=apply(yPredPerMod,1,mean)   ##mean because it is a regression problem
    return(list(yPred=yPred,
                yPredPerMod=yPredPerMod))  ###return
  }
 ###########################
##when it is classification problem
  else if (class(MUVRclassobject)[3]=='Classification') {
    yPredPerMod=array(dim=c(nrow(newdata),             ##observations
                            length(levels(MUVRclassobject$inData$Y)),  ###number of levels in Y
                            length(MUVRclassobject$outModels)),###length of list length(nouter*nRep)
                      dimnames=list(paste('observation',1:nrow(newdata),sep=''),
                                    levels(MUVRclassobject$inData$Y),
                                    paste('model',
                                          1:length(MUVRclassobject$outModels),
                                           sep='')))
    n=0
    for(r in 1:nRep) {
      for(i in 1:nOuter) {
        n=n+1
        mod=MUVRclassobject$outModels[[n]][[modNum]]
        ###n =r*i in the end. mod is an output of summary
        ###it is a rfoutmin
        if (method=='PLS') {
          cat('\nNot yet implemented')
          return(NULL)
##########################################################################
###Can pls da be used here in y prediction
###
##############################################################################

        } else {
          if ((!rownames(mod$importance)%in%colnames(newdata))) {
            cat('\nMismatch variable names in model',n,': Return NULL')
            return(NULL)
          } else {
#################################################################################################
###this does not seem work in ranger. The ranger does not have such name as importance
            X=subset(newdata,
                     select=rownames(mod$importance))

            yPredPerMod[,,n]=predict(mod,       ##random forest model vote type
                                     newdata=X,
                                     type='vote')  #
          }
        }
      }
    }
    yPred=apply(yPredPerMod,c(1,2),mean)  ###there is 56 matrix of predicted y (obser,yclass),mean for each value of the matrix

    yClass=levels(MUVRclassobject$inData$Y)[apply(yPred,1,which.max)]  ###which is the position of the max output

    names(yClass)=paste('observation',1:nrow(newdata),sep='')

    return(list(yClass=yClass,
                yPred=yPred,
                yPredPerMod=yPredPerMod))
  } else {
    cat('\nNot yet implemented')
  }
}
