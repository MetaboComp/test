#' Plot validation metric 
#'
#' Produces a plot of validation metric vs number of variables in model (inner segment)
#' @param MVObj An object of class `MVObject`
#'
#' @return A plot
#' @export
plotVAL=function(MVObj) {
  VAL=MVObj$VAL$VAL
  metric=MVObj$VAL$metric
  count=as.numeric(colnames(VAL))
  nRep=dim(VAL)[3]
  plot(count,count,ylim=range(VAL),xlim=range(count),log='x',type='n',bty='l',ylab=metric,xlab='Number of variables (log scale)')
  for (r in 1:nRep) {
    matlines(count,t(VAL[,,r]),type='l',lty=1,col='grey')
  }
  lines(count,apply(VAL,2,mean))
  for (i in 1:3) {
    abline(v=MVObj$nVar[i],lty=i)
  }
  legend('topright',legend=c('MinModel','MidModel','MaxModel'),lty=1:3,col='black',bty='n')
}
