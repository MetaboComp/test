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
  performance_miss<-list()
  performance_auc<-list()
  performance_Q2<-list()
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

  if(!is.null(MUVRclassObject$miss)){performance_miss[[i]]<-MUVRclassObject$miss}
  if(!is.null(MUVRclassObject$auc)){performance_auc[[i]]<-MUVRclassObject$auc}
  if(!is.null(MUVRclassObject$fitMetric$Q2)){performance_Q2[[i]]<-MUVRclassObject$fitMetric$Q2}


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

if(length(performance_auc)!=0) {

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
classperformance1<-suppressWarnings(speed_and_performance(Xotu,Yotu,method="RF",fitness="MISS",methParam = customParamsa(rfMethod="randomForest",NZV=T,oneHot=T)))
classperformance2<-suppressWarnings(speed_and_performance(Xotu,Yotu,method="RF",fitness="MISS",methParam = customParamsa(rfMethod="ranger",NZV=T,oneHot=T)))


regrperformance1<-suppressWarnings(speed_and_performance(XRVIP,YR,ID=IDR,method="RF",fitness="RMSEP",DA=F,methParam = customParamsa(rfMethod="randomForest",NZV=T,oneHot=T)))
regrperformance2<-suppressWarnings(speed_and_performance(XRVIP,YR,ID=IDR,method="RF",fitness="RMSEP",DA=F,methParam = customParamsa(rfMethod="ranger",NZV=T,oneHot=T)))

MLperformance1<-suppressWarnings(speed_and_performance(X=crispEM, ML=TRUE, method="RF",fitness="MISS",methParam = customParamsa(rfMethod="randomForest",NZV=T,oneHot=T)))
MLperformance2<-suppressWarnings(speed_and_performance(X=crispEM, ML=TRUE, method="RF",fitness="MISS",methParam = customParamsa(rfMethod="ranger",NZV=T,oneHot=T)))


################################compare randomForest and ranger in origin and onehotencoding
###classification
factor_variable1<-as.factor(c(rep(33,15),rep(44,7),rep(55,7)))
factor_variable2<-as.factor(c(rep("AB",4),rep("CD",4),rep("EF",6),
                              rep("GH",4),rep("IJ",4),rep("KL",7)))
factor_variable3<-as.factor(c(rep("Tessa",5),rep("Olle",5),rep("Yan",5),
                              rep("Calle",9),rep("Elisa",5)))

character_variable1<-c(rep("one",4),rep("two",4),rep("three",4),
                       rep("four",4),rep("five",4),rep("six",4),rep("seven",5))
character_variable2<-c(rep("yes",8),rep("no",7),
                       rep("yes",7),rep("no",7))
character_variable3<-c(rep("Hahahah",29))

logical_variable1<-c(rep(TRUE,7),rep(FALSE,6),rep(TRUE,7),rep(FALSE,6),rep(TRUE,1),rep(FALSE,2))
logical_variable2<-c(rep(TRUE,8),rep(FALSE,7),rep(TRUE,7),rep(FALSE,7))

z<-data.frame(row.names=1:29)
X<-cbind(z,Xotu,
         factor_variable1,factor_variable2,factor_variable3,
         character_variable1,character_variable2,character_variable3,
         logical_variable1,logical_variable2)

Xone<-onehotencoding(X)
classperformance1_origin<-suppressWarnings(speed_and_performance(X,Yotu,method="RF",fitness="MISS",methParam = customParamsa(rfMethod="randomForest",NZV=T,oneHot=T)))
classperformance1_onehot<-suppressWarnings(speed_and_performance(Xone,Yotu,method="RF",fitness="MISS",methParam = customParamsa(rfMethod="randomForest",NZV=T,oneHot=T)))
classperformance2_origin<-suppressWarnings(speed_and_performance(X,Yotu,method="RF",fitness="MISS",methParam = customParamsa(rfMethod="ranger",NZV=T,oneHot=T)))
classperformance2_onehot<-suppressWarnings(speed_and_performance(Xone,Yotu,method="RF",fitness="MISS",methParam = customParamsa(rfMethod="ranger",NZV=T,oneHot=T)))

#################################################################################
###regression
factor_variable1<-as.factor(c(rep(33,105),rep(44,3),rep(55,4)))
factor_variable2<-as.factor(c(rep("AB",20),rep("CD",10),rep("EF",30),
                              rep("GH",15),rep("IJ",25),rep("KL",12)))
factor_variable3<-as.factor(c(rep("Tessa",25),rep("Olle",30),rep("Yan",12),
                              rep("Calle",25),rep("Elisa",20)))

character_variable1<-c(rep("one",16),rep("two",16),rep("three",16),
                       rep("four",16),rep("five",16),rep("six",16),rep("seven",16))
character_variable2<-c(rep("yes",28),rep("no",28),
                       rep("yes",28),rep("no",28))
character_variable3<-c(rep("Hahahah",112))

logical_variable1<-c(rep(TRUE,16),rep(FALSE,16),rep(TRUE,16),rep(FALSE,16),rep(TRUE,16),rep(FALSE,32))
logical_variable2<-c(rep(TRUE,28),rep(FALSE,28),rep(TRUE,28),rep(FALSE,28))


z<-data.frame(row.names=1:112)
X<-cbind(z,XRVIP,
         factor_variable1,factor_variable2,factor_variable3,
         character_variable1,character_variable2,character_variable3,
         logical_variable1,logical_variable2)
Xone<-onehotencoding(X)
regrperformance1_origin<-suppressWarnings(speed_and_performance(X,YR,ID=IDR,method="RF",fitness="RMSEP",DA=F,methParam = customParamsa(rfMethod="randomForest",NZV=T,oneHot=T)))
regrperformance1_onehot<-suppressWarnings(speed_and_performance(Xone,YR,ID=IDR,method="RF",fitness="RMSEP",DA=F,methParam = customParamsa(rfMethod="randomForest",NZV=T,oneHot=T)))
regrperformance2_origin<-suppressWarnings(speed_and_performance(X,YR,ID=IDR,method="RF",fitness="RMSEP",DA=F,methParam = customParamsa(rfMethod="ranger",NZV=T,oneHot=T)))
regrperformance2_onehot<-suppressWarnings(speed_and_performance(Xone,YR,ID=IDR,method="RF",fitness="RMSEP",DA=F,methParam = customParamsa(rfMethod="ranger",NZV=T,oneHot=T)))

######################################################################
#####Multilevel


factor_variable1<-as.factor(c(rep(33,7),rep(44,7),rep(55,7)))
factor_variable2<-as.factor(c(rep("AB",4),rep("CD",3),rep("EF",3),
                              rep("GH",5),rep("IJ",3),rep("KL",3)))
factor_variable3<-as.factor(c(rep("Tessa",4),rep("Olle",4),rep("Yan",5),
                              rep("Calle",4),rep("Elisa",4)))

character_variable1<-c(rep("one",3),rep("two",3),rep("three",3),
                       rep("four",3),rep("five",3),rep("six",3),rep("seven",3))
character_variable2<-c(rep("yes",6),rep("no",5),
                       rep("yes",5),rep("no",5))
character_variable3<-c(rep("Hahahah",21))

logical_variable1<-c(rep(TRUE,4),rep(FALSE,4),rep(TRUE,4),rep(FALSE,3),rep(TRUE,3),rep(FALSE,3))
logical_variable2<-c(rep(TRUE,6),rep(FALSE,5),rep(TRUE,5),rep(FALSE,5))


z<-data.frame(row.names=1:21)
X<-cbind(z,crispEM,
         factor_variable1,factor_variable2,factor_variable3,
         character_variable1,character_variable2,character_variable3,
         logical_variable1,logical_variable2)




Xone<-onehotencoding(X)
MLperformance1_origin<-suppressWarnings(speed_and_performance(X=X, ML=TRUE, method="RF",fitness="MISS",methParam = customParamsa(rfMethod="randomForest",NZV=T,oneHot=T)))
MLperformance1_onehot<-suppressWarnings(speed_and_performance(X=Xone, ML=TRUE, method="RF",fitness="MISS",methParam = customParamsa(rfMethod="randomForest",NZV=T,oneHot=T)))


MLperformance2_origin<-suppressWarnings(speed_and_performance(X=X, ML=TRUE, method="RF",fitness="MISS",methParam = customParamsa(rfMethod="ranger",NZV=T,oneHot=T)))
MLperformance2_onehot<-suppressWarnings(speed_and_performance(X=Xone, ML=TRUE, method="RF",fitness="MISS",methParam = customParamsa(rfMethod="ranger",NZV=T,oneHot=T)))























