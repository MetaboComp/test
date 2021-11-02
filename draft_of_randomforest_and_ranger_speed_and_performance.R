
###1
###test with first 0 variables in XRVIP, 8 random factor/character/logical variables in x, y is YR
###Because Y value is not a factor , it runs regression in randomForest and ranger
X<-data.frame(row.names=1:112)
X<-cbind(X,
         factor_variable1,factor_variable2,factor_variable3,
         character_variable1,character_variable2,character_variable3,
         logical_variable1,logical_variable2)
X<-as.data.frame(X)

##1_1 normal data
randomForest1_1_1<-randomForest(x=X,
                                y=YR,
                                ntree=500,
                                mtry=5)
pred1_1_1<-predict(randomForest1_1_1,X)
time1_1_1<-microbenchmark(randomForest1_1_1<-randomForest(x=X,
                                                          y=YR,
                                                          ntree=500,
                                                          mtry=5),times=100,unit="s")
summarytime1_1_1<-summary(time1_1_1)[,4:5]
ranger1_2_1<-ranger(x=X,
                    y=YR,
                    num.trees = 500,
                    mtry = 5,
                    num.threads = 1)
pred1_2_1<-predict(ranger1_2_1,X)
time1_2_1<-microbenchmark(ranger1_2_1<-ranger(x=X,
                                              y=YR,
                                              num.trees = 500,
                                              mtry = 5,
                                              num.threads = 1),times=100,unit="s")
summarytime1_2_1<-summary(time1_2_1)[,4:5]


###1_2 one hot encoding
x=onehotencoding(X)
randomForest1_1_2<-randomForest(x=x,
                                y=YR,
                                ntree=500,
                                mtry=5)
pred1_1_2<-predict(randomForest1_1_2,x)
time1_1_2<-microbenchmark(randomForest1_1_2<-randomForest(x=x,
                                                          y=YR,
                                                          ntree=500,
                                                          mtry=5),times=100,unit="s")
summarytime1_1_2<-summary(time1_1_2)[,4:5]
ranger1_2_2<-ranger(x=x,
                    y=YR,
                    num.trees = 500,
                    mtry = 5,
                    num.threads = 1)
pred1_2_2<-predict(ranger1_2_2,x)
time1_2_2<-microbenchmark(ranger1_2_2<-ranger(x=x,
                                              y=YR,
                                              num.trees = 500,
                                              mtry = 5,
                                              num.threads = 1),times=100,unit="s")
summarytime1_2_2<-summary(time1_2_2)[,4:5]
time1<-rbind(summarytime1_1_1,summarytime1_2_1,summarytime1_1_2,summarytime1_2_2)
rownames(time1)<-c("XrandomForest","Xranger","oneHotrandomForest","oneHot ranger")
time1


###2
###test with first 8 variables in XRVIP, 8 random factor/character/logical variables in x, y is YR
###Because Y value is not a factor , it runs regression in randomForest and ranger
X<-data.frame(row.names=1:112)
X<-cbind(X,XRVIP[,1:8],
         factor_variable1,factor_variable2,factor_variable3,
         character_variable1,character_variable2,character_variable3,
         logical_variable1,logical_variable2)
##2_1 normal data
randomForest2_1_1<-randomForest(x=X,
                                y=YR,
                                ntree=500,
                                mtry=5)
pred2_1_1<-predict(randomForest2_1_1,X)
time2_1_1<-microbenchmark(randomForest2_1_1<-randomForest(x=X,
                                                          y=YR,
                                                          ntree=500,
                                                          mtry=5),times=100,unit="s")
summarytime2_1_1<-summary(time2_1_1)[,4:5]
ranger2_2_1<-ranger(x=X,
                    y=YR,
                    num.trees = 500,
                    mtry = 5,
                    num.threads = 1)
pred2_2_1<-predict(ranger2_2_1,X)
time2_2_1<-microbenchmark(ranger2_2_1<-ranger(x=X,
                                              y=YR,
                                              num.trees = 500,
                                              mtry = 5,
                                              num.threads = 1),times=100,unit="s")
summarytime2_2_1<-summary(time2_2_1)[,4:5]


###2_2 one hot encoding
x=onehotencoding(X)
randomForest2_1_2<-randomForest(x=x,
                                y=YR,
                                ntree=500,
                                mtry=5)
pred2_1_2<-predict(randomForest2_1_2,x)
time2_1_2<-microbenchmark(randomForest2_1_2<-randomForest(x=x,
                                                          y=YR,
                                                          ntree=500,
                                                          mtry=5),times=100,unit="s")
summarytime2_1_2<-summary(time2_1_2)[,4:5]
ranger2_2_2<-ranger(x=x,
                    y=YR,
                    num.trees = 500,
                    mtry = 5,
                    num.threads = 1)
pred2_2_2<-predict(ranger2_2_2,x)
time2_2_2<-microbenchmark(ranger2_2_2<-ranger(x=x,
                                              y=YR,
                                              num.trees = 500,
                                              mtry = 5,
                                              num.threads = 1),times=100,unit="s")
summarytime2_2_2<-summary(time2_2_2)[,4:5]
time2<-rbind(summarytime2_1_1,summarytime2_2_1,summarytime2_1_2,summarytime2_2_2)
rownames(time2)<-c("XrandomForest","Xranger","oneHotrandomForest","oneHot ranger")
time2



###3
###test with first 80 variables in XRVIP, 8 random factor/character/logical variables in x, y is YR
###Because Y value is not a factor , it runs regression in randomForest and ranger
X<-data.frame(row.names=1:112)
X<-cbind(X,XRVIP[,1:80],
         factor_variable1,factor_variable2,factor_variable3,
         character_variable1,character_variable2,character_variable3,
         logical_variable1,logical_variable2)


##3_1 normal data
randomForest3_1_1<-randomForest(x=X,
                                y=YR,
                                ntree=500,
                                mtry=5)
pred3_1_1<-predict(randomForest3_1_1,X)
time3_1_1<-microbenchmark(randomForest3_1_1<-randomForest(x=X,
                                                          y=YR,
                                                          ntree=500,
                                                          mtry=5),times=100,unit="s")
summarytime3_1_1<-summary(time3_1_1)[,4:5]
ranger3_2_1<-ranger(x=X,
                    y=YR,
                    num.trees = 500,
                    mtry = 5,
                    num.threads = 1)
pred3_2_1<-predict(ranger3_2_1,X)
time3_2_1<-microbenchmark(ranger3_2_1<-ranger(x=X,
                                              y=YR,
                                              num.trees = 500,
                                              mtry = 5,
                                              num.threads = 1),times=100,unit="s")
summarytime3_2_1<-summary(time3_2_1)[,4:5]


###3_2 one hot encoding
x=onehotencoding(X)
randomForest3_1_2<-randomForest(x=x,
                                y=YR,
                                ntree=500,
                                mtry=5)
pred3_1_2<-predict(randomForest3_1_2,x)
time3_1_2<-microbenchmark(randomForest3_1_2<-randomForest(x=x,
                                                          y=YR,
                                                          ntree=500,
                                                          mtry=5),times=100,unit="s")
summarytime3_1_2<-summary(time3_1_2)[,4:5]
ranger3_2_2<-ranger(x=x,
                    y=YR,
                    num.trees = 500,
                    mtry = 5,
                    num.threads = 1)
pred3_2_2<-predict(ranger3_2_2,x)
time3_2_2<-microbenchmark(ranger3_2_2<-ranger(x=x,
                                              y=YR,
                                              num.trees = 500,
                                              mtry = 5,
                                              num.threads = 1),times=100,unit="s")
summarytime3_2_2<-summary(time3_2_2)[,4:5]
time3<-rbind(summarytime3_1_1,summarytime3_2_1,summarytime3_1_2,summarytime3_2_2)
rownames(time3)<-c("XrandomForest","Xranger","oneHotrandomForest","oneHot ranger")
time3

###4
###test with first 800 variables in XRVIP, 8 random factor/character/logical variables in x, y is YR
###Because Y value is not a factor , it runs regression in randomForest and ranger

X<-data.frame(row.names=1:112)
X<-cbind(X,XRVIP[,1:800],
         factor_variable1,factor_variable2,factor_variable3,
         character_variable1,character_variable2,character_variable3,
         logical_variable1,logical_variable2)


##4_1 normal data
randomForest4_1_1<-randomForest(x=X,
                                y=YR,
                                ntree=500,
                                mtry=5)
pred4_1_1<-predict(randomForest4_1_1,X)

time4_1_1<-microbenchmark(randomForest4_1_1<-randomForest(x=X,
                                                          y=YR,
                                                          ntree=500,
                                                          mtry=5),times=100,unit="s")
summarytime4_1_1<-summary(time4_1_1)[,4:5]

ranger4_2_1<-ranger(x=X,
                    y=YR,
                    num.trees = 500,
                    mtry = 5,
                    num.threads = 1)
pred4_2_1<-predict(ranger4_2_1,X)

time4_2_1<-microbenchmark(ranger4_2_1<-ranger(x=X,
                                              y=YR,
                                              num.trees = 500,
                                              mtry = 5,
                                              num.threads = 1),times=100,unit="s")
summarytime4_2_1<-summary(time4_2_1)[,4:5]


###4_2 one hot encoding
x=onehotencoding(X)
randomForest4_1_2<-randomForest(x=x,
                                y=YR,
                                ntree=500,
                                mtry=5)
pred4_1_2<-predict(randomForest4_1_2,x)

time4_1_2<-microbenchmark(randomFores4_1_2<-randomForest(x=x,
                                                         y=YR,
                                                         ntree=500,
                                                         mtry=5),times=100,unit="s")
summarytime4_1_2<-summary(time4_1_2)[,4:5]

ranger4_2_2<-ranger(x=x,
                    y=YR,
                    num.trees = 500,
                    mtry = 5,
                    num.threads = 1)
pred4_2_2<-predict(ranger4_2_2,x)

time4_2_2<-microbenchmark(ranger4_2_2<-ranger(x=x,
                                              y=YR,
                                              num.trees = 500,
                                              mtry = 5,
                                              num.threads = 1),times=100,unit="s")
summarytime4_2_2<-summary(time4_2_2)[,4:5]

time4<-rbind(summarytime4_1_1,summarytime4_2_1,summarytime4_1_2,summarytime4_2_2)
rownames(time4)<-c("XrandomForest","Xranger","oneHotrandomForest","oneHot ranger")
time4

###5
###test with all variables in XRVIP, 8 random factor/character/logical variables in x, y is YR
###Because Y value is not a factor , it runs regression in randomForest and ranger
X<-data.frame(row.names=1:112)
X<-cbind(X,XRVIP,
         factor_variable1,factor_variable2,factor_variable3,
         character_variable1,character_variable2,character_variable3,
         logical_variable1,logical_variable2)

##5_1 normal data
randomForest5_1_1<-randomForest(x=X,
                                y=YR,
                                ntree=500,
                                mtry=5)
pred5_1_1<-predict(randomForest5_1_1,X)
time5_1_1<-microbenchmark(randomForest5_1_1<-randomForest(x=X,
                                                          y=YR,
                                                          ntree=500,
                                                          mtry=5),times=100,unit="s")
summarytime5_1_1<-summary(time5_1_1)[,4:5]
ranger5_2_1<-ranger(x=X,
                    y=YR,
                    num.trees = 500,
                    mtry = 5,
                    num.threads = 1)
pred5_2_1<-predict(ranger5_2_1,X)
time5_2_1<-microbenchmark(ranger5_2_1<-ranger(x=X,
                                              y=YR,
                                              num.trees = 500,
                                              mtry = 5,
                                              num.threads = 1),times=100,unit="s")
summarytime5_2_1<-summary(time5_2_1)[,4:5]

###5_2 one hot encoding
x=onehotencoding(X)
randomForest5_1_2<-randomForest(x=x,
                                y=YR,
                                ntree=500,
                                mtry=5)
pred5_1_2<-predict(randomForest5_1_2,x)
time5_1_2<-microbenchmark(randomFores5_1_2<-randomForest(x=x,
                                                         y=YR,
                                                         ntree=500,
                                                         mtry=5),times=100,unit="s")
summarytime5_1_2<-summary(time5_1_2)[,4:5]
ranger5_2_2<-ranger(x=x,
                    y=YR,
                    num.trees = 500,
                    mtry = 5,
                    num.threads = 1)
pred5_2_2<-predict(ranger5_2_2,x)
time5_2_2<-microbenchmark(ranger5_2_2<-ranger(x=x,
                                              y=YR,
                                              num.trees = 500,
                                              mtry = 5,
                                              num.threads = 1),times=100,unit="s")
summarytime5_2_2<-summary(time5_2_2)[,4:5]
time5<-rbind(summarytime5_1_1,summarytime5_2_1,summarytime5_1_2,summarytime5_2_2)
rownames(time5)<-c("XrandomForest","Xranger","oneHotrandomForest","oneHot ranger")
time5
time<-cbind(time1,time2,time3,time4,time5)
colnames(time)<-paste(rep(c("mean","median"),5),c(rep("0",2),rep("8",2),rep("80",2),rep("800",2),rep("1147",2)))

time








#'This is to test the speed and performance of randomForest and ranger
#'
#'
#'
#'@param x
#'@param y
#'@param xtest
#'@param ytest
#'@param
#'@param
#'@param
#'@param times   ###in microbenchmark
#'
#'
#'
#'@return
#'@example data ("airquality");airquality<-as.data.frame(airquality)
#'@example X_training<-airquality[,2:6][1:125,]
#'@example
#'@example X_testing<-airquality[,2:6][126:153,]
#'@example Y_training<-airquality[,1][1:125]
#'@example Y_testing<-airquality[,1][126:153]
#'@example
#'@example
#'
#'@export
#'
#'
#'



testspeed<-function(x,y,xtest,ytest,times)
{library(microbenchmark)
  library(randomForest)
  library(ranger)


}















