

####use the example of XRVIP, Xotu, adding 8 random factor/character/logical variables to test for the speed and out of bag error
####note that only factor and numeric variables can be used in randomforest
library(microbenchmark)
library(randomForest)
library(ranger)
#The data used must be imputed first to remove all the NA value




#'This is to test the speed of randomForest and ranger
#'
#'@param x x is a data frame of many numeric variables
#'@param y
#'@param z is a non- numeric variables matrix
#'@param ntree
#'@param mtry
#'@param times  In microbenchmark
#'@return time. it is a list
#'@export



timesummary<-function(x,
                      y,
                      z,       ###z is the added non numeric variables
                      times=100,
                      mtry=5,
                      ntree=500)
{s<-vector()

for(i in 0:ncol(x)%/%ncol(z))
{if((ncol(x)%/%ncol(z))>=(10^(i))){s[i+1]<-ncol(z)*10^(i)}
  else{s[i+1]=ncol(x)
  break}}
s<-c(0,s)

X<-list()
xone<-list()
timeframe<-list()
time<-array(data=NA,dim=c(4,2,length(s)),
            dimnames=list(c("XrandomForest","Xranger","oneHotrandomForest","oneHotranger"),
                          c("mean","median"),
                          paste0(s,"  varibles in x")))

for(r in 1:length(s)){
  if(r==1){X[[r]]<-data.frame(row.names=1:nrow(x))
  X[[r]]<-cbind(X[[r]],z)
  }
  if(r>1)
  {X[[r]]<-data.frame(row.names=1:nrow(x))
  X[[r]]<-cbind(x[,1:s[r]],z)
  }
}

for(r in 1:length(s))
{randomForest0_1_1<-randomForest(x=X[[r]],
                                 y=y,
                                 ntree=ntree,
                                 mtry=mtry)
pred0_1_1<-predict(randomForest0_1_1,X[[r]])

time0_1_1<-microbenchmark(randomForest0_1_1<-randomForest(x=X[[r]],
                                                          y=y,
                                                          ntree=ntree,
                                                          mtry=mtry),times=times,unit="s")
summarytime0_1_1<-summary(time0_1_1)[,4:5]

ranger0_2_1<-ranger(x=X[[r]],
                    y=y,
                    num.trees = ntree,
                    mtry = mtry,
                    num.threads = 1)
pred0_2_1<-predict(ranger0_2_1,X[[r]])

time0_2_1<-microbenchmark(ranger0_2_1<-ranger(x=X[[r]],
                                              y=y,
                                              num.trees =ntree,
                                              mtry = mtry,
                                              num.threads = 1),times=times,unit="s")
summarytime0_2_1<-summary(time0_2_1)[,4:5]


###5_2 one hot encoding
xone[[r]]=onehotencoding(X[[r]])
randomForest0_1_2<-randomForest(x=xone[[r]],
                                y=y,
                                ntree=ntree,
                                mtry=mtry)
pred0_1_2<-predict(randomForest0_1_2,xone[[r]])

time0_1_2<-microbenchmark(randomFores0_1_2<-randomForest(x=xone[[r]],
                                                         y=y,
                                                         ntree=ntree,
                                                         mtry=mtry),times=times,unit="s")
summarytime0_1_2<-summary(time0_1_2)[,4:5]

ranger0_2_2<-ranger(x=xone[[r]],
                    y=y,
                    num.trees = ntree,
                    mtry = mtry,
                    num.threads = 1)
pred0_2_2<-predict(ranger0_2_2,xone[[r]])

time0_2_2<-microbenchmark(ranger0_2_2<-ranger(x=xone[[r]],
                                              y=y,
                                              num.trees = ntree,
                                              mtry =mtry,
                                              num.threads = 1),times=100,unit="s")
summarytime0_2_2<-summary(time0_2_2)[,4:5]

timeframe[[r]]<-as.data.frame(rbind(summarytime0_1_1,summarytime0_2_1,summarytime0_1_2,summarytime0_2_2))

for(n in 1:4)
{for(m in 1:2){time[n,m,r]=timeframe[[r]][n,m]
}
}
}

return(time)
}





####Use XRVIP and YR and z, z is an 8 non numeric variables data frame that has the same rows as XRVIP and YR
##factor_variable1<-ordered(c(rep(33,105),rep(44,3),rep(55,4)),levels=c("44","55","33"))

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
z<-cbind(z,
         factor_variable1,factor_variable2,factor_variable3,
         character_variable1,character_variable2,character_variable3,
         logical_variable1,logical_variable2)

regrtime<-timesummary(XRVIP,YR,z)



############now test classification for randomForest and ranger
####Use Xotu and Yotu and z, z is an 8 non numeric variables data frame that has the same rows as Xotu and Yotu



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
z<-cbind(z,
         factor_variable1,factor_variable2,factor_variable3,
         character_variable1,character_variable2,character_variable3,
         logical_variable1,logical_variable2)
Xotu_1<-as.matrix(Xotu)
classtime<-timesummary(Xotu_1,Yotu,z)






#################################
####compare1 is made when using mosquito data,using NZV=T,MISS.


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
##Xotu or Xotu_1
Xotu_1<-as.matrix(Xotu)
X<-cbind(z,Xotu_1,
         factor_variable1,factor_variable2,factor_variable3,
         character_variable1,character_variable2,character_variable3,
         logical_variable1,logical_variable2)
X_1<-as.matrix(X)

#####note, in this function, I used customParams() created in my own branch, to differentiate with the one in the master branch. I marked it as customParamsa()
compare1<-microbenchmark(MUVR_randomForest1<-MUVR(X = Xotu_1, Y = Yotu, nRep = nRep, nOuter = nOuter, varRatio = varRatio,
                                                  fitness = "MISS", method = "RF", methParam = customParamsa(method = "randomForest",NZV=T),
                                                  modReturn = T),
                         MUVR_ranger1<-MUVR(X = Xotu_1, Y = Yotu, nRep = nRep, nOuter = nOuter, varRatio = varRatio,
                                            fitness = "MISS", method = "RF", methParam = customParamsa(method = "ranger",NZV=T),
                                            modReturn = T) ,
                         times=2)


compare1_havefactor<-microbenchmark(MUVR_randomForest1<-MUVR(X = X, Y = Yotu, nRep = nRep, nOuter = nOuter, varRatio = varRatio,
                                                  fitness = "MISS", method = "RF", methParam = customParamsa(method = "randomForest",NZV=T),
                                                  modReturn = T),
                         MUVR_ranger1<-MUVR(X = X, Y =Yotu, nRep = nRep, nOuter = nOuter, varRatio = varRatio,
                                            fitness = "MISS", method = "RF", methParam = customParamsa(method = "ranger",NZV=T),
                                            modReturn = T) ,
                         times=2)
compare1_havefactor_onehot<-microbenchmark(MUVR_randomForest1<-MUVR(X = onehotencoding(X), Y = Yotu, nRep = nRep, nOuter = nOuter, varRatio = varRatio,
                                                             fitness = "MISS", method = "RF", methParam = customParamsa(method = "randomForest",NZV=T),
                                                             modReturn = T),
                                    MUVR_ranger1<-MUVR(X = onehotencoding(X), Y =Yotu, nRep = nRep, nOuter = nOuter, varRatio = varRatio,
                                                       fitness = "MISS", method = "RF", methParam = customParamsa(method = "ranger",NZV=T),
                                                       modReturn = T) ,
                                    times=2)

####can be run




##################################################################################################################
####compare 2 is regrModel
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

compare2<-microbenchmark(MUVR_randomForest2<-MUVR(X = XRVIP, Y = YR, nRep = nRep, nOuter = nOuter, varRatio = varRatio,
                                                 fitness = "RMSEP",DA=F, method = "RF", methParam = customParamsa(method = "randomForest",NZV=T),
                                                 modReturn = T),
                         MUVR_ranger2<-MUVR(X = XRVIP, Y = YR, nRep = nRep, nOuter = nOuter, varRatio = varRatio,
                                           fitness = "RMSEP", DA=F,method = "RF", methParam = customParamsa(method = "ranger",NZV=T),
                                           modReturn = T),
                         times=2)

compare2_havefactor<-microbenchmark(MUVR_randomForest2<-MUVR(X = X, Y = YR, nRep = nRep, nOuter = nOuter, varRatio = varRatio,
                                                                    fitness = "RMSEP",DA=F, method = "RF", methParam = customParamsa(method = "randomForest",NZV=T),
                                                                    modReturn = T),
                                           MUVR_ranger2<-MUVR(X = X, Y = YR, nRep = nRep, nOuter = nOuter, varRatio = varRatio,
                                                              fitness = "RMSEP", DA=F,method = "RF", methParam = customParamsa(method = "ranger",NZV=T),
                                                              modReturn = T),
                                           times=2)

compare2_havefactor_onehot<-microbenchmark(MUVR_randomForest2<-MUVR(X = onehotencoding(X), Y = YR, nRep = nRep, nOuter = nOuter, varRatio = varRatio,
                                                                    fitness = "RMSEP",DA=F, method = "RF", methParam = customParamsa(method = "randomForest",NZV=T),
                                                                    modReturn = T),
                                           MUVR_ranger2<-MUVR(X = onehotencoding(X), Y = YR, nRep = nRep, nOuter = nOuter, varRatio = varRatio,
                                                              fitness = "RMSEP", DA=F,method = "RF", methParam = customParamsa(method = "ranger",NZV=T),
                                                              modReturn = T),
                                           times=2)
####################################################################################################################################
##### compare 3 is MLModel
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



compare3<-microbenchmark(MUVR_randomForest3<-MUVR(X = crispEM, ML=T,nRep = nRep, nOuter = nOuter, varRatio = varRatio,
                                                  method = "RF", methParam = customParamsa(method = "randomForest",NZV=T),
                                                 modReturn = T),
                         MUVR_ranger3<-MUVR(X = crispEM, ML=T, nRep = nRep, nOuter = nOuter, varRatio = varRatio,
                                            method = "RF", methParam = customParamsa(method = "ranger",NZV=T),
                                           modReturn = T),
                         times=2)

compare3_havefactor<-microbenchmark(MUVR_randomForest3<-MUVR(X = X, ML=T,nRep = nRep, nOuter = nOuter, varRatio = varRatio,
                                                                    method = "RF", methParam = customParamsa(method = "randomForest",NZV=T),
                                                                    modReturn = T),
                                           MUVR_ranger3<-MUVR(X =X, ML=T, nRep = nRep, nOuter = nOuter, varRatio = varRatio,
                                                              method = "RF", methParam = customParamsa(method = "ranger",NZV=T),
                                                              modReturn = T),
                                           times=2)



compare3_havefactor_onehot<-microbenchmark(MUVR_randomForest3<-MUVR(X = onehotencoding(X), ML=T,nRep = nRep, nOuter = nOuter, varRatio = varRatio,
                                                  method = "RF", methParam = customParamsa(method = "randomForest",NZV=T),
                                                  modReturn = T),
                         MUVR_ranger3<-MUVR(X = onehotencoding(X), ML=T, nRep = nRep, nOuter = nOuter, varRatio = varRatio,
                                            method = "RF", methParam = customParamsa(method = "ranger",NZV=T),
                                            modReturn = T),
                         times=2)



