#'Speeding and performance test for different models
#'
#'@param X
#'@param Y
#'@param ID
#'@param scale
#'@param nRep
#'@param nOuter
#'@param nInner
#'@param varRatio = 0.75,
#'@param DA
#'@param fitness
#'@param method
#'@param methParam
#'@param ML
#'@param modReturn
#'@param logg
#'@param parallel
#'@param repeatMUVR  This is to decide how many MUVR to repeat(since each MUVR may give a differnent value)
#'
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
                               modReturn = FALSE,
                               logg = FALSE,
                               parallel = TRUE,
                               repeatMUVR=2

)
  {
  if(method=="RF"){performance_miss<-list()
                   performance_auc<-list()}
    if(method=="PLS"){performance_Q2<-list()}
  MUVRclassObject<-list()
  time1<-vector()
  time2<-vector()
  for(i in 1:repeatMUVR){


  time1[i] <- proc.time()[3]
  cl=makeCluster(nCore)
  registerDoParallel(cl)
    MUVRclassObject<-MUVR(X,Y,ID,scale,nRep,nOuter, nInner,varRatio,DA ,fitness,method,methParam,
        ML,modReturn,logg ,parallel)
    stopCluster(cl)
  time2[i] <- proc.time()[3]

  if(class(MUVRclassObject)[3]=="RF")
  {
  performance_miss[[i]]<-MUVRclassObject$miss
  performance_auc[[i]]<-MUVRclassObject$auc
  }
  else if(class(MUVRclassObject)[3]=="PLS")
  {
  performance_Q2[[i]]<-MUVRclassObject$fitMetric$Q2}


  }


  timeperrun<-time2-time1

if(exists("performance_Q2")){

performance_Q2_mean<-matrix(0L,1,3)
 for(s in 1:3){
   for(i in 1:length(performance_Q2))
   {performance_Q2_mean[1,s]<-performance_Q2_mean[1,s]+performance_Q2[[i]][s]}
 }

colnames(performance_Q2_mean)<-c("min","mid","max")
performance_Q2_mean<-performance_Q2_mean/length(performance_Q2)
}

if(exists("performance_miss")){
  performance_miss_mean<-matrix(0L,1,3)
  for(s in 1:3){
    for(i in 1:length(performance_miss))
    {performance_miss_mean[1,s]<-performance_miss_mean[1,s]+performance_miss[[i]][s]}
  }

  colnames(performance_miss_mean)<-c("min","mid","max")
  performance_miss_mean<-performance_miss_mean/length(performance_miss)

}

if(exists("performance_auc")) {

  if(!is.matrix(performance_auc[[1]])){

  performance_auc_mean<-matrix(0L,1,3)
  for(s in 1:3){
    for(i in 1:length(performance_auc))
    {performance_auc_mean[1,s]<-performance_auc_mean[1,s]+performance_auc[[i]][s]}
  }

  colnames(performance_auc_mean)<-c("min","mid","max")
  performance_auc_mean<-performance_auc_mean/length(performance_auc)

  }else if(is.matrix(performance_auc[[1]]))
          {performance_auc_mean<-matrix(0L,3,dim(performance_auc[[1]])[2])
            for(k in 1:3){
              for(s in 1:dim(performance_auc[[1]])[2]){
                 for(i in 1:length(performance_auc))
              {performance_auc_mean[k,s]<-performance_auc_mean[k,s]+performance_auc[[i]][k,s]}
            }
            }
            rownames(performance_auc_mean)<-c("min","mid","max")
            colnames(performance_auc_mean)<-colnames(performance_auc[[1]])
            performance_auc_mean<-performance_auc_mean/length(performance_auc)}


}

  result<-list()
  result$minperrun<-timeperrun/60
  result$mean_minperrun<-mean(timeperrun)/60
  if(exists("performance_miss_mean")){
  result$performance_miss<-performance_miss_mean}
  if(exists("performance_auc_mean")){
  result$performance_auc<-performance_auc_mean}
  if(exists("performance_Q2_mean")){
  result$performance_Q2<-performance_Q2_mean}

  return(result)
}

