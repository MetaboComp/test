#' Plot stability of selected variables and prediction fitness as a function of number of repetitions 
#'
#' @param MVObject MUVR object
#' @param model 'min' (default), 'mid' or 'max'
#' @param VAll Option of specifying which variables (i.e. names) to consider as reference set. Defaults to variables selected from the `model` of the `MVObject`
#' @param nVarLim Option of specifying upper limit for number of variables
#' @param missLim Option of specifying upper limit for number of misclassifications
#'
#' @return Plot of number of variables, proportion of variables overlapping with reference and prediction accuracy (Q2 for regression; MISS otherwise) as a function of number of repetitions. 
#' @export
plotStability=function(MVObject,model='min',VAll,nVarLim,missLim) {
  regr=any(class(MVObject)=='Regression')
  DA=MVObject$inData$DA
  ML=MVObject$inData$ML
  Y=MVObject$inData$Y
  nModel=ifelse(model=='min',1,ifelse(model=='mid',2,3))
  nVar=round(MVObject$nVar[nModel])
  if(missing(VAll)) VAll=names(sort(MVObject$VIP[,nModel])[1:nVar]) # Final selection of variables
  nRep=MVObject$inData$nRep
  nV=VA=miss=r2=q2=numeric(nRep)
  for (i in 1:nRep) {
    nV[i]=round(mean(MVObject$nVarPerRep[[nModel]][1:i]))
    VA[i]=sum(names(sort(rowMeans(MVObject$VIPPerRep[[nModel]]))[1:nV[i]])%in%VAll)
    if(DA) {
      preds=MVObject$yPredPerRep[[nModel]][,,1:i]
      preds=apply(preds,c(1,2),mean)
      miss[i]=sum(levels(Y)[apply(preds,1,which.max)]!=Y)
    } else {
      preds=MVObject$yPredPerRep[[nModel]][,1:i,drop=F]
      preds=rowMeans(preds)
      PRESS=sum((Y-preds)^2)
      TSS=sum((Y-mean(Y))^2)
      q2[i]=1-(PRESS/TSS)
    }
    if(ML) {
      class=ifelse(preds<0,-1,1)
      miss[i]=sum(class!=Y)
    }
  }
  VA=VA/length(VAll)
  if(missing(nVarLim)) {
    pot=10^floor(log10(max(nV)))
    nVarLim=ceiling(max(nV)/pot)*pot
  }
  par(mar=c(4,8,0,ifelse(ML,8,4))+.5)
  plot(nV,ylim=c(0,nVarLim),type='l',xlab='Number of repetitions',ylab='',axes=F)
  box(bty='u')
  axis(1)
  axis(2,line=4)
  title(ylab='Number of selected variables',line = 6)
  par(new=T)
  plot(VA,type='l',ylim=c(0,1),lty=2,col='red',axes=F,xlab='',ylab='')
  axis(2)
  title(ylab='Proportion of selected variables',line = 2,col.lab="red")
  if(DA | ML) {
    if(missing(missLim)) missLim=length(Y)
    par(new=T)
    plot(miss,ylim=c(0,missLim),type='l',col='blue',lty=3,axes=F,xlab='',ylab='')
    axis(4)
    mtext('Number of misclassifications',side=4,line = 2,col="blue",cex=par()$cex)
  }
  if(regr | ML) {
    par(new=T)
    plot(q2,ylim=c(0,1),type='l',col='green',lty=4,axes=F,xlab='',ylab='')
    axis(4,line=ifelse(regr,NA,4))
    mtext('Q2',side=4,line = ifelse(regr,2,6),col="green",cex=par()$cex)
  }
}
