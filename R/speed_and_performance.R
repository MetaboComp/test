#'Speeding and performance test for different models
#'
#'@param X X in MUVR
#'@param Y Y in MUVR
#'@param ID ID in MUVR
#'@param scale scale in MUVR
#'@param nRep nRep in MUVR
#'@param nOuter nOuter in MUVR
#'@param nInner nInner in MUVR
#'@param varRatio = 0.75 varRatio in MUVR
#'@param DA DA in MUVR
#'@param fitness fitness in MUVR
#'@param method method in MUVR
#'@param methParam methParam in MUVR
#'@param ML ML in MUVR
#'@param modReturn modReturn in MUVR
#'@param logg logg in MUVR
#'@param parallel parallel X in MUVR
#'@param repeatMUVR  This is to decide how many MUVR to repeat(since each MUVR may give a differnent value)
#'@param MUVRversion "MUVR" or "MUVRoriginal"
#'@return time   average time for each MUVR
#'@return performance   average miss auc and rmsep
#'@export
#'
#'
speed_and_performance<-function(X,
                               Y,
                               ID,
                               scale = TRUE,
                               nRep = 5,
                               nOuter = 6,
                               nInner,
                               varRatio = 0.75,
                               DA = FALSE,
                               fitness = c('AUROC', 'MISS', 'BER', 'RMSEP'),   ###ber is not integrated in it yet
                               method = c('PLS',' RF'),
                               methParam,
                               ML = FALSE,
                               modReturn = T,
                               logg = FALSE,
                               parallel = TRUE,
                               repeatMUVR=2,
                               MUVRversion=c("MUVR","MUVRoriginal")

)
  {
  performance_miss<-list()

  performance_Q2<-list()
  performance_nVar<-list()

  MUVRclassObject<-list()
  time1<-vector()
  time2<-vector()



  if(MUVRversion=="MUVR"){
  for(i in 1:repeatMUVR){

  cat("The",i,"th time out of",repeatMUVR,"times.", "\n")
  time1[i] <- proc.time()[3]
  nCore=8
  cl=makeCluster(nCore)
  registerDoParallel(cl)
    MUVRclassObject<-MUVR(X,Y,ID,scale,nRep,nOuter, nInner,varRatio,DA ,fitness,method,methParam,
        ML,modReturn,logg ,parallel)
    stopCluster(cl)
  time2[i] <- proc.time()[3]

  if(!is.null(MUVRclassObject$miss)){performance_miss[[i]]<-MUVRclassObject$miss}

  if(!is.null(MUVRclassObject$fitMetric$Q2)){performance_Q2[[i]]<-MUVRclassObject$fitMetric$Q2}
  if(!is.null(MUVRclassObject$nVar)){performance_nVar[[i]]<-MUVRclassObject$nVar}


  }
  }

  if(MUVRversion=="MUVRoriginal"){
    for(i in 1:repeatMUVR){

      cat("The",i,"th time out of",repeatMUVR,"times.", "\n")
      time1[i] <- proc.time()[3]
      nCore=8
      cl=makeCluster(nCore)
      registerDoParallel(cl)
      MUVRclassObject<-MUVRoriginal(X,Y,ID,scale,nRep,nOuter, nInner,varRatio,DA ,fitness,method,methParam,
                            ML,modReturn,logg ,parallel)
      stopCluster(cl)
      time2[i] <- proc.time()[3]

      if(!is.null(MUVRclassObject$miss)){performance_miss[[i]]<-MUVRclassObject$miss}

      if(!is.null(MUVRclassObject$fitMetric$Q2)){performance_Q2[[i]]<-MUVRclassObject$fitMetric$Q2}
      if(!is.null(MUVRclassObject$nVar)){performance_nVar[[i]]<-MUVRclassObject$nVar}


    }
  }



  timeperrun<-time2-time1

if(length(performance_Q2)!=0){

performance_Q2_mean<-matrix(0L,1,3)
 for(s in 1:3){
   for(i in 1:length(performance_Q2))
   {performance_Q2_mean[1,s]<-performance_Q2_mean[1,s]+performance_Q2[[i]][s]}
 }

colnames(performance_Q2_mean)<-c("min","mid","max")
performance_Q2_mean<-performance_Q2_mean/length(performance_Q2)
}

if(length(performance_miss)!=0){
  performance_miss_mean<-matrix(0L,1,3)
  for(s in 1:3){
    for(i in 1:length(performance_miss))
    {performance_miss_mean[1,s]<-performance_miss_mean[1,s]+performance_miss[[i]][s]}
  }

  colnames(performance_miss_mean)<-c("min","mid","max")
  performance_miss_mean<-performance_miss_mean/length(performance_miss)

}

if(length(performance_nVar)!=0){
    performance_nVar_mean<-matrix(0L,1,3)
    for(s in 1:3){
      for(i in 1:length(performance_nVar))
      {performance_nVar_mean[1,s]<-performance_nVar_mean[1,s]+performance_nVar[[i]][s]}
    }

    colnames(performance_nVar_mean)<-c("min","mid","max")
    performance_nVar_mean<-performance_nVar_mean/length(performance_nVar)

  }





  result<-list()
  result$minperrun<-timeperrun/60
  result$mean_minperrun<-mean(timeperrun)/60
  if(exists("performance_miss_mean")){
  result$performance_miss_mean<-performance_miss_mean}

  if(exists("performance_Q2_mean")){
  result$performance_Q2_mean<-performance_Q2_mean}
  if(exists("performance_nVar_mean")){
    result$performance_nVar_mean<-performance_nVar_mean}

  if(exists("performance_miss")){
    result$performance_miss<-performance_miss}

  if(exists("performance_Q2")){
    result$performance_Q2<-performance_Q2}
  if(exists("performance_nVar")){
    result$performance_nVar<-performance_nVar}

  return(result)
}






















