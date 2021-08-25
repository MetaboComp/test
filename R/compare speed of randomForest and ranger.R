

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








