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
#'@param MUVRversion
#'@return time   average time for each MUVR
#'@return result  average miss auc and rmsep
#'@export
#'
#'
test_for_regrandmulti<-function(X,
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
                                  MUVRversion=c("MUVR","MUVRoriginal"))

{
  ###regression data random forest
  test<-list()
  test_nVar<-list()
  test_Q2<-list()
  test_time<-list()

  if(MUVRversion=="MUVR"){
    for(i in 1:repeatMUVR)
    {cat(i,"th time")
      test[[i]]<-MUVR(X,Y,ID,scale,nRep,nOuter, nInner,varRatio,DA ,fitness,method,methParam,
                      ML,modReturn,logg ,parallel)

    }
  }

  if(MUVRversion=="MUVRoriginal"){
    for(i in 1:repeatMUVR)
    {cat(i,"th time", "\n")
      test[[i]]<-MUVRoriginal(X,Y,ID,scale,nRep,nOuter, nInner,varRatio,DA ,fitness,method,methParam,
                              ML,modReturn,logg ,parallel)

    }
  }


  test_Q2_mean_min<-0
  test_nVar_mean_min<-0

  test_Q2_mean_mid<-0
  test_nVar_mean_mid<-0

  test_Q2_mean_max<-0
  test_nVar_mean_max<-0

  test_time_mean<-0

  for(i in 1:repeatMUVR)
  {test_Q2[[i]]<-test[[i]]$fitMetric$Q2
  test_nVar[[i]]<-test[[i]]$nVar
  test_time[[i]]<-test[[i]]$calcMins}



  for(i in 1:repeatMUVR)
  {test_Q2_mean_min<-test_Q2_mean_min+test_Q2[[i]][1]
  test_nVar_mean_min<-test_nVar_mean_min+test_nVar[[i]][1]

  test_Q2_mean_mid<-test_Q2_mean_mid+test_Q2[[i]][2]
  test_nVar_mean_mid<-test_nVar_mean_mid+test_nVar[[i]][2]

  test_Q2_mean_max<-test_Q2_mean_max+test_Q2[[i]][3]
  test_nVar_mean_max<-test_nVar_mean_max+test_nVar[[i]][3]
  test_time_mean<-test_time_mean+test_time[[i]]

  }
  test_Q2_mean<-c(test_Q2_mean_min,test_Q2_mean_mid,test_Q2_mean_max)/repeatMUVR
  test_nVar_mean<-c(test_nVar_mean_min,test_nVar_mean_mid,test_nVar_mean_max)/repeatMUVR
  test_time_mean<-test_time_mean/repeatMUVR
  result<-list()
  result$test_Q2_mean<-test_Q2_mean
  result$test_nVar_mean<-test_nVar_mean
  result$test_time_mean<-test_time_mean
  result$test_Q2<-test_Q2
  result$test_nVar<-test_nVar
  result$test_time<-test_time
  return(result)
}
