#' Make permutations with data and default settings from an actual MUVR object
#'
#' This function will extract data and parameter settings from a MUVR object and run standard permutations.
#' This will fit a standard case of multivariate predictive modelling in either a regression, classification or multilevel case.
#' However, if an analysis has a complex sample dependency which requires constrained permutation of your response vector
#' or if a variable pre-selection is performed for decreased computational burden, then permutaion loops should be constructed manually.
#' In those cases, View(permutations) can be a first start from which to build custom solutions for permutation analysis.
#' @param MUVRclassObject A 'MUVR' class obvject
#' @param nPerm number of permutations to run
#' @param nRep number of repetitions for each permutation (defaults to value of actual model)
#' @param nOuter number of outer validation segments for each permutation (defaults to value of actual model)
#' @param varRatio varRatio for each permutation (defaults to value of actual model)
#' @param parallel whether to run calculations using parallel processing - which requires registered backend (defaults to value of actual model)
#'
#' @return permutation_output: A permutation matrix with permuted fitness statistics (nrow=nPerm and ncol=3 for min/mid/max)
#' @return permutation_type: Either AUROC,Q2 or MISS or BER
#' @export
#'

## library(MUVR)
##  library(doParallel)
##  nCore=detectCores()-1
##  cl=makeCluster(nCore)
##  registerDoParallel(cl)
##  nRep=2*nCore
##  varRatio=.75
##  nOuter=6
##  nPerm=50
## R12ML=MUVR(X=mlr12,ML=T,nRep=nRep,nOuter=nOuter,varRatio=varRatio,method='RF')
##  permR12=permutations(R12ML)
## stopCluster(cl)
## permutationPlot(R12ML,permR12)

permutations=function(MUVRclassObject,
                      nPerm=50,
                      nRep,
                      nOuter,
                      varRatio,
                      parallel,
                      permutation_type= c('AUROC', 'MISS', 'RMSEP',"BER")) {

##substitute() is often paired with deparse(). That function takes the result of substitute(), an expression,
##and turns it into a character vector.

  name=deparse(substitute(MUVRclassObject))   ###the output is "MUVRclassObject"


  X=MUVRclassObject$inData$X
  Y=MUVRclassObject$inData$Y
  ID=MUVRclassObject$inData$ID
  DA=MUVRclassObject$inData$DA
  ML=MUVRclassObject$inData$ML
  scale=MUVRclassObject$inData$scale
  fitness=MUVRclassObject$inData$fitness
  method=MUVRclassObject$inData$method
  methParam=MUVRclassObject$inData$methParam

  if(missing(nRep)) nRep=MUVRclassObject$inData$nRep

  if(missing(nOuter)) {
    nOuter=MUVRclassObject$inData$nOuter
    nInner=MUVRclassObject$inData$nInner
  } else nInner=nOuter-1

  if(missing(varRatio)) varRatio=MUVRclassObject$inData$varRatio
  if(missing(parallel)) parallel=MUVRclassObject$inData$parallel

  if(ML) {
    nSamp=nrow(X)/2
    X=X[1:nSamp,]
    ID=ID[1:nSamp]
  }

###########################################################################################################################3
#I added a new variabler permutation type to include AUROC
  if(!missing(permutation_type)){
    if(permutation_type!="AUROC"&permutation_type!="MISS"&permutation_type!="Q2"&permutation_type!="BER")
    {stop("permutation_type is not correct")}
    if(permutation_type=="Q2"&!any(class(MUVRclassObject)=='Regression'))
    {stop("Classification and Multilevel must use AUROC or MISS for permutation")}
    if(permutation_type!="Q2"&any(class(MUVRclassObject)=='Regression'))
    {stop("Regression must use Q2 for permutation")}
      }
  if(missing(permutation_type)){
    if (any(class(MUVRclassObject)=='Regression')){permutation_type!="Q2"}
        else{permutation_type="MISS"}
  }



  if(permutation_type=="Q2"|permutation_type=="MISS"|permutation_type=="BER"){
  startTime=proc.time()[3]
  permutation_output=matrix(ncol=3,nrow=nPerm)   ##row is permutation number col is 3 (min,mid,max)
  colnames(permutation_output)=c('Min','Mid','Max')
  rownames(permutation_output)=paste("permutation",1:nPerm)
##################################################################################################################################
###I also added permutation test for auc when it is classfication
###hBER could be added too when we want to integrate it into the permutationtest
 for (p in 1:nPerm) {   ####these is to repeat random selection for nPerm times
    cat('\n"',name,'" permutation ',p,' of ',nPerm,'\n',sep = '')

    if (ML) YPerm=sample(c(-1,1),size=nSamp,replace=TRUE)   ##choose the same size sample randomly bootstrapping
      else YPerm=sample(Y)      ##sample(x, size, replace = FALSE, prob = NULL)
                                ## default for size is the number of items inferred from the first argument,

    permMod=MUVR(X=X,
                 Y=YPerm,
                 ID=ID,
                 scale=scale,
                 DA=DA,
                 ML=ML,
                 nRep=nRep,
                 nOuter=nOuter,
                 nInner=nInner,
                 varRatio=varRatio,
                 fitness=fitness,
                 method=method,
                 methParam=methParam,
                 parallel=parallel)

    if (any(class(MUVRclassObject)=='Regression')) permutation_output[p,]=permMod$fitMetric$Q2

    else if (permutation_type=="MISS"){permutation_output[p,]=permMod$miss}
    else {permutation_output[p,]=permMod$ber}

    nowTime=proc.time()[3]

    timePerRep=(nowTime-startTime)/p

    timeLeft=(timePerRep*(nPerm-p))/60

    cat('\nEstimated time left:',timeLeft,'mins\n\n')
 }
  lst<-list(permutation_type,
            permutation_output)
  names(lst)<-c("permutation_type",
                 "permutation_output")
  return(lst)
  }
##############################################################################################333
#When permutation is BER




########################################
#when permutation is AUROC
  if(permutation_type=="AUROC"){
  startTime=proc.time()[3]
  permutation_output=array(NA,dim=c(nPerm,3,ncol(MUVRclassObject$auc)),
             dimnames=list(c(paste("permutation",1:nPerm)),
                           c('Min','Mid','Max'),
                           colnames(MUVRclassObject$auc))) ###if use auc this is a list

  for(s in 1:ncol(MUVRclassObject$auc)){
  for (p in 1:nPerm) {   ####these is to repeat random selection for nPerm times
    cat('\n"',"group",s,name,'" permutation ',p,' of ',nPerm,'\n',sep = '')

    if (ML) YPerm=sample(c(-1,1),size=nSamp,replace=TRUE)   ##choose the same size sample randomly bootstrapping
    else YPerm=sample(Y)      ##sample(x, size, replace = FALSE, prob = NULL)
    ## default for size is the number of items inferred from the first argument,

    permMod=MUVR(X=X,
                 Y=YPerm,
                 ID=ID,
                 scale=scale,
                 DA=DA,
                 ML=ML,
                 nRep=nRep,
                 nOuter=nOuter,
                 nInner=nInner,
                 varRatio=varRatio,
                 fitness=fitness,     ##what if AUROC and MISS
                 method=method,
                 methParam=methParam,
                 parallel=parallel)

    permutation_output[p,,s]=t(permMod$auc[,s])


    nowTime=proc.time()[3]

    timePerRep=(nowTime-startTime)/((s-1)*nPerm+p)

    timeLeft=(timePerRep*(nPerm*ncol(MUVRclassObject$auc)-((s-1)*nPerm+p)))/60

    cat('\nEstimated time left:',timeLeft,'mins\n\n')
  }
  }
  lst<-list(permutation_type,
            permutation_output)
  names(lst)<-c("permutation_type",
                "permutation_output")
  return(lst)

  }


}
