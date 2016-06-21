#' Plot MV object
#'
#' @param MVObj An `MVobject` obtained from the MVWrap function
#' @param model What type of model to plot ('min', 'mid' or 'max'). Defaults to 'mid'.
#' @return A pdf with plots of results from multivariate predictions
#' @export
plotMV=function(MVObj,model='mid') {
  if (!any(class(MVObj)=='MVObject')) {
    cat('\nWrong object class: Return NULL')
    return(NULL)
  }
  modNum=ifelse(model=='min',1,ifelse(model=='mid',2,3))
  Y=MVObj$inData$Y
  nSamp=length(Y)
  par(mar=c(4,4,0,0)+.5)
  if (class(MVObj)[3]=='Regression') {
    YP=MVObj$yPred[,modNum]
    YPR=MVObj$yPredPerRep[[modNum]]
    matplot(Y,YPR,pch=20,xlab='Original Y',ylab='Predicted Y',col='grey',bty='l',cex=0.5)
    points(Y,YP,pch=20)
    reg=lm(YP~Y)
    abline(reg)
    legend('topleft',legend=c(paste('Model R2 =',signif(MVObj$fitMetric$R2[modNum],3)),paste('Model Q2 =',signif(MVObj$fitMetric$Q2[modNum],3))),bty='n')
  } else if (class(MVObj)[3]=='Classification') {
    YP=MVObj$yPred[[modNum]]
    YPR=MVObj$yPredPerRep[[modNum]]
    classes=1:length(levels(Y))
    classNudge=0.2*((classes-mean(classes))/(mean(classes)-1))
    plot(1:nSamp,Y,type='n',ylim=range(YPR),xlab='Sample number',ylab='Class prediction probability')
    for(cl in classes) {
      matpoints((1:nSamp)+classNudge[cl],YPR[,cl,],pch=20,col=cl+1,cex=0.5)
      points((1:nSamp)+classNudge[cl],YP[,cl],pch=20,col=cl+1)
    }
    for (li in 1:(nSamp+1)) {
      abline(v=li-.5,lty=3,col='grey')
    }
    yClass=MVObj$yClass[,modNum]
    whichWrong=which(yClass!=Y)
    wrongClass=as.numeric(yClass[whichWrong])
    for (w in 1:length(wrongClass)) {
      points(whichWrong[w]+classNudge[wrongClass[w]],YP[whichWrong[w],wrongClass[w]],cex=2)
    }
    legend('topleft',legend=c(levels(Y),'misclassified'),pch=c(rep(16,length(classes)),1),col=c(classes+1,1),cex=0.8,pt.cex=c(rep(0.5,length(classes)),2))
  } else {
    YP=MVObj$yPred[,modNum]
    YPR=MVObj$yPredPerRep[[modNum]]
    matplot(YPR,1:nSamp,pch=20,col='grey',cex=0.5,ylim=c(nSamp,1),ylab='Sample number',xlab='Predicted Y')
    points(YP,1:nSamp,pch=20,col='black')
    abline(h=nSamp/2+0.5,lty=2)
    abline(v=0,lty=2)
  }
}