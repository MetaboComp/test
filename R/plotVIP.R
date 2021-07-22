#' Plots variable importance ranking in MUVR object
#'
#' Regardless of MV core method, variables are sorted by rank, where lower is better. `plotVIRank` produces boxplots of variable rankings for all model repetitions.
#'
#' @param MUVRclassObject An MUVR class object
#' @param n Number of top ranking variables to plot (defaults to those selected by MUVR)
#' @param cut Optional value to cut length of variable names to `cut` number of characters
#' @param model Which model to choose ('min', 'mid' {default} or 'max')
#'
#' @return Barplot of variable rankings (lower is better)
#' @export
plotVIRank=function(MUVRclassObject,n,model='mid',cut) {
  nModel=ifelse(model=='min',1,ifelse(model=='mid',2,3))
  nFeat=round(MUVRclassObject$nVar[nModel])
  if(missing(n)) n=nFeat
  VIRank=MUVRclassObject$VIRank[,nModel]
  VIRankRep=MUVRclassObject$VIRankPerRep[[nModel]]
####################################################################################################
  VIRankRep=VIRankRep[order(VIRank),][1:n,]    ###make VIRank a sequence For example, if the first variable is ranking 27th
                                               ###Then it is put in the 27th row of the new dataframe
#I don't know if here is correct. The variable with lowest VIRank has the highest order(VIRank)
#Assume the OTU_28 (variable 29th) has the lowest VIRank,Then it has the higest order Ex 1609
#Here the new VIRankRep choose the original row of variable 1609 put in the new 29th row
# However, variable 1609 has nothing to do with the sequencing?
#
########################################################################################################

  if(n>nFeat) {
    VIRankRep=rbind(VIRankRep[1:nFeat,],rep(NA,ncol(VIRankRep)),VIRankRep[(nFeat+1):n,])
    col=rep(c('yellow','grey'),c(nFeat,(n-nFeat+1)))    ###repeat yellow nFeat times, repear grey n-nFeat+1 times
  } else col=NULL
  ####
  VIRankRep=VIRankRep[nrow(VIRankRep):1,]             ###reverse
  col=rev(col)                                          ###reverse the sequence of col
  boxplot(t(VIRankRep),                            ### row is rep,column is variable
          horizontal=T,                            ###	logical indicating if the boxplots should be horizontal;
          axes=F,                                  ### Do not add any axes
          col=col)                                 ###when outside n it is grey, if inside it is yellow
  axis(1)                                         ###manually add x axis
  labels=rownames(VIRankRep)                    ##add labels

  if(!missing(cut)) labels=substring(labels,1,cut)  ###Extract or replace substrings in a character vector
                                                    ###xxx <- c("asfef", "qwerty", "yuiop[", "b", "stuff.blah.yech")
                                                    ###substr(xxx, 2, 5)
                                                    ###[1] "sfef" "wert" "uiop" ""     "tuff"
  axis(2,                                      ###manually add y axis,
       las=1,                                  ###label is horzontal
       at=1:nrow(VIRankRep),                   ###label position
       labels=labels)                          ###name of label
  if (n>nFeat) abline(h=(n-nFeat+1))           ###add a line that separate the variable inside and outside n
  box(bty='o')                                       ###add a box that has o shape
}
