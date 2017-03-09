#' PLS biplot
#'
#' Makes a biplot of a fitted object (e.g. from a MUVR with PLS core)
#' @param fit A PLS fit (e.g. from MVObject$Fit[[2]])
#' @param comps Which components to plot
#' @param xCol (Optionsal) Continuous vector for grey scale gradient of observation (sample) color (e.g. Y vector in regression analysis)
#' @param labPlSc Boolean to plot observation (sample) names (defaults to TRUE)
#' @param labs (Optional) Labels names
#' @param vars Which variables to plot (names in rownames(loadings))
#' @param labPlLo Boolean to plot variable names (defaults to TRUE)
#' @param pchSc Plotting character for observation scores
#' @param colSc Colors for observation scores (only if xCol omitted)
#' @param colLo Colors for variable loadings (defaults to red)
#'
#' @return A PLS biplot
#' @export
biplotPLS=function(fit,comps=1:2,xCol,labPlSc=TRUE,labs,vars,labPlLo=TRUE,pchSc=16,colSc,colLo=2) {
  scores=fit$variates$X[,comps]
  loads=fit$loadings$X[,comps]
  if(missing(vars)) vars=rownames(loads)
  loads=loads[rownames(loads)%in%vars,]
  nSamp=nrow(scores)
  nVar=nrow(loads)
  if(missing(xCol)) {
    if (missing(colSc)) {
      colSc=rep(1,nSamp) 
      legPlot=FALSE
    } else {
      colScLeg=colSc
      colSc=as.factor(colSc)
      legPlot=TRUE
    }
  } else {
    x.col=10+round(85*((max(xCol)-xCol)/(diff(range(xCol)))))
    colSc=paste("gray",x.col,sep="")
    legPlot=TRUE
  }
  rSc=c(apply(scores,2,min),apply(scores,2,max))
  rLo=c(apply(loads,2,min),apply(loads,2,max))
  scale=max(abs(rSc))/max(abs(rLo))
  lim=1.2*max(scores)
  plot(scores,col=colSc,pch=pchSc,ylim=c(-lim,lim),xlim=c(-lim,lim),bty='l',xlab=paste('Component',comps[1]),ylab=paste('Component',comps[2]),type='n')
  abline(h=0,v=0,lty=2,col='grey')
  if (labPlSc) {
    if (missing(labs)) labs=rownames(scores) 
    # cat(labs)
    text(scores,as.character(labs),pos=3)
  }
  arrows(rep(0,nrow(loads)),rep(0,nrow(loads)),scale*loads[,1],scale*loads[,2],col=colLo)
  if(labPlLo) text(1.1*scale*loads,rownames(loads),col=colLo)
  if (legPlot) {
    if (missing(xCol)) {
      legend('topleft',legend=unique(colScLeg),pch=16,col=unique(colSc),bty='n')
    } else {
      whUnik=!duplicated(xCol)
      unik=xCol[whUnik]
      cols=colSc[whUnik][order(unik)]
      unik=sort(unik)
      if (length(unik)>10) {
        k=(length(unik)-1)/5
        n=1+k*0:5
        cols=cols[n]
        unik=unik[n]
      }
      legend('topleft',legend=signif(unik,3),fill=cols,bty='n')
    }
  }
  points(scores,col=colSc,pch=pchSc)
}