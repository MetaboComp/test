#' Extract autoselected variables from MUVR model object
#'
#' @param MUVRclassObject An object of MUVR class
#' @param model Which model to use ("min", "mid" (default), or "max")
#'
#' @return Data frame with order, name and average rank of variables (`order`, `name` & `rank`)
#' @export
getVIRank=function(MUVRclassObject,model='mid') {
  nMod=ifelse(model=='min',1,ifelse(model=='mid',2,3))
  nVar=round(MUVRclassObject$nVar[nMod])
  VIRanks=sort(MUVRclassObject$VIP[,nMod])[1:nVar]       ###sequencing them and take the first few of them
  ##In the MUVR class object, it is still called VIP because it is not ranked
  ##only by this getVIRank,it is ranked
  VIRanks=data.frame(order=1:nVar,name=names(VIRanks),rank=VIRanks)
  VIRanks$name=as.character(VIRanks$name)
  return(VIRanks)
}
