#' MUVR: "Multivariate modelling with Unbiased Variable selection"
#'
#' Repeated double cross validation with tuning of variables in the inner loop.
#'
#' @param X Predictor variables. NB: Variables (columns) must have names/unique identifiers. NAs not allowed in data. For multilevel, only the positive half of the difference matrix is specified.
#' @param Y Response vector (Dependent variable). For classification, a factor (or character) variable should be used. For multilevel, Y is calculated automatically.
#' @param ID Subject identifier (for sampling by subject; Assumption of independence if not specified)
#' @param NZV Boolean for whether to filter out near zero variance variables (defaults to TRUE)
#' @param nRep Number of repetitions of double CV. (Defaults to 5)
#' @param nOuter Number of outer CV loop segments. (Defaults to 6)
#' @param nInner Number of inner CV loop segments. (Defaults to nOuter-1)
#' @param DA Boolean for Classification (discriminant analysis) (By default, if Y is numeric -> DA=FALSE. If Y is factor (or character) -> DA=TRUE)
#' @param fitness Fitness function for model tuning (choose either 'AUROC' or 'MISS' (default) for classification; or 'RMSEP' (default) for regression.)
#' @param methParam List with parameter settings for specified MV method (see function code for details)
#' @param ML Boolean for multilevel analysis (defaults to FALSE)
#' @param modReturn Boolean for returning outer segment models (defaults to FALSE). Setting modReturn=TRUE is required for making MUVR predictions using predMV().
#' @param parallel Boolean for whether to perform `foreach` parallel processing (Requires a registered parallel backend; Defaults to `TRUE`)
#' @param alow alpha tuning: lowest value of alpha
#' @param ahigh alpha tuning: highest value of alpha
#' @param astep alpha tuning: number of alphas to try from low to high
#' @param alog alpha tuning: Whether to space tuning of alpha in logarithmic scale (TRUE; default) or normal/arithmetic scale (FALSE)
#' @param keep A group of confounders that you want to manually set as non-zero
#' @param ...
#'
#' @return A MUVR object
#' @export
rdCVnet <- function(X,   ## X should be a dataframe
                    Y,
                    ID,
                    alow=1e-5,  ## lowest value of alpha
                    ahigh=1,    ##
                    astep=11,   ##  number of alphas to try from low to high
                    alog=TRUE,
                    nRep=5,
                    nOuter=6,
                    nInner,
                    NZV=TRUE,
                    DA=FALSE,
                    fitness=c('AUROC', 'MISS', 'BER', 'RMSEP'),
                    methParam,
                    ML=FALSE,
                    modReturn=FALSE,
                    parallel=TRUE,
                    keep=NULL,
                    ...){
  library(glmnet)
  X_original<-X

  # methParams
  if (missing(methParam)){
    methParam <- rdcvNetParams()
    }

  ## If DA, multinomial, if not gaussian


  if(methParam$oneHot==T){
    if(any(class(X)%in% c("data.frame"))) {
      X <- onehotencoding(X)
      cat("X is transformed to a matrix by onehotencoding.", "\n")
    }
  } else {X <- as.matrix(X)}
  # Remove nearZeroVariance variables

  modelReturn <- list(call=match.call())

  if (NZV) {
    nzv <- MUVR::nearZeroVar(X)
    if (length(nzv$Position)>0) {
      modelReturn$nzv=colnames(X)[nzv$Position]
      X <- X[,-nzv$Position]
      cat('\n',length(nzv$Position),'variables with near zero variance detected -> removed from X and stored under $nzv')
    }
  }

  # Number of samples and variables
  nSamp <- nrow(X)
  nVar <- nVar0 <- ncol(X)
  # Sample identifiers
  if (missing(ID)) {
    cat('\nMissing ID -> Assume all unique (i.e. sample independence)')
    ID <- 1:nSamp
  }
  if (ML) {
    X <- rbind(X,-X)
    if (missing(Y)) Y <- rep(-1,nSamp)
    Y <- c(Y,-Y)
    nSamp <- 2 * nSamp
    ID <- c(ID,ID)
    DA <- FALSE
    fitness <- 'MISS'
    cat('\nMultilevel -> Regression on (-1,1) & fitness=MISS')
  }

  if(missing(methParam)){
    methParam=rdcvNetParams()
    }
  # Call in packages
  library(pROC) # for roc and auroc
  library(magrittr) # for pipe operator
  library(foreach) # Parallel processing

  # Set up for parallel/serial
  if (parallel) {"%doVersion%" <- get("%dopar%")
  } else{ "%doVersion%" <- get("%do%")}

  # Initialise modelReturn with function call


  # Start timer
  start.time <- proc.time()[3]

  # Rough check indata
  if (length(dim(X))!=2) {stop('\nWrong format of X matrix.\n')}
  if (is.null(colnames(X))) {stop('\nNo column names in X matrix.\n')}







  # Sort out internal modelling parameters
  if (missing(nInner)) {nInner <- nOuter - 1}

  # DA / Classification
  if (is.character(Y)) {Y <- factor(Y)}
  if (is.factor(Y)) {
    cat('\nY is factor -> Classification (',
        length(unique(Y)),
        ' classes)',
        sep='')
    DA <- TRUE
  }
  if (is.numeric(Y) & DA) {
    Y <- as.factor(Y)
    cat('\nDA=TRUE -> Y as factor -> Classification (',
        length(unique(Y)),
        ' classes)',
        sep='')
  }
  if(DA) {
    methParam$family <- 'multinomial'
  }
  # Check fitness criterion
  # This may not be needed!
  if (missing(fitness)) {
    if (DA) {
      fitness <- 'BER'
      cat('\nMissing fitness -> BER')
    } else {
      fitness <- 'RMSEP'
      cat('\nMissing fitness -> RMSEP')
    }
  }

  # Set up for multilevel analysis

  # No Missingness allowed
  if (any(is.na(X)) | any(is.na(Y))) {
    stop('\nNo missing values allowed in X or Y data.\n')
    }
  if (!is.null(dim(Y))) {
    cat('\nY is not a vector: Return NULL')
    return(NULL)
  }

  # Sanity check
  if (nrow(X)!=length(Y)) {
    cat('\nMust have same nSamp in X and Y: Return NULL')
    return(NULL)
  }


  ## Store indata in list for later model return
  InData <- list(X=X,
                 Y=Y,
                 ID=ID,
                 alow=alow,
                 ahigh=ahigh,
                 astep=astep,
                 nRep=nRep,
                 nOuter=nOuter,
                 nInner=nInner,
                 DA=DA,
                 fitness=fitness,
                 methParam=methParam,
                 ML=ML,
                 parallel=parallel)

  ## Sort sampling based on subjects and not index
  unik <- !duplicated(ID)  # boolean of unique IDs
  unikID <- ID[unik]
  if (DA) {
    if(nOuter > min(table(Y))) {
      warning('\nnOuter is larger than your smallest group size. Consider lowering your nOuter to min(table(Y)).',call.=TRUE)}
    unikY <- Y[unik]  # Counterintuitive, but needed for groupings by Ynames
    Ynames <- sort(unique(Y))  # Find groups
    groups <- length(Ynames) # Number of groups
    groupID <- list()  # Allocate list for indices of groups
    for (g in 1:groups) {
      groupID[[g]] <- unikID[unikY==Ynames[g]]  # Find indices per group
    }
  }

  # Allocate prediction variables
  if (DA) {
    # Predictions per repetition
    yPredSeg <- matrix(nrow=length(Y),
                       ncol=length(levels(Y)),
                       dimnames=list(ID,levels(Y)))
  } else {
    # Predictions per repetition
    yPredSeg <- numeric(length(Y))
  }
##########################################0############################################################################################
  #########################################################################################3
  # tuning values for alpha
  if (alog) {   ## Whether to space tuning of alpha in logarithmic scale (TRUE; default) or normal/arithmetic scale (FALSE)
    avect <- seq(from=log10(alow),
                 to=log10(ahigh),
                 length.out=astep) # logarithmic
    avect <- 10^avect
  } else { avect <- seq(from=alow,
                      to=ahigh,
                      length.out=astep) }# Arithmetic


  ## Choose package/core algorithm according to chosen method
  packs <- c('pROC', 'glmnet', 'magrittr')
  exports <- c('vectSamp', 'uniqDASamp', 'foldVector')

  ############################
  ## Start repetitions
  # reps <- list() # 2+2 lines for pseudomanual troubleshooting
  # for (r in 1:nRep){
  reps <- foreach(r=1:nRep,
                  .packages=packs,
                  .export=exports) %doVersion% {
    # r <- 1
    # r <- r + 1
    if (modReturn) {outMod <- list()}
    cat('\n','   Repetition ',r,' of ',nRep,':',sep='')
    # Sampling into holdout segments
    if (DA & identical(unikID,ID)) {
      allTest <- uniqDASamp(Y, ID, nOuter)

      ##a problem is that in some groups there are only 2 levels
    } else {
      allTest <- vectSamp(unikID,n=nOuter)
    }
    # Allocate result variables per repetition
    if(DA) {
      coefSeg <- array(0,
                       dim=c(nVar0,   ## number of variables
                             length(levels(Y)),  ## number of groups
                             nOuter),
                       dimnames = list(colnames(X),
                                       levels(Y),
                                       paste0('Segment',1:nOuter)))
    } else {
      coefSeg <- matrix(0,
                        nrow=nVar0,
                        ncol=nOuter)
      rownames(coefSeg) <- colnames(X)
      colnames(coefSeg) <- paste0('Segment',1:nOuter)
    }
    alphaSeg <- lambdaSeg <- nonZeroSeg <- numeric(nOuter)

    if(DA){nonZeroSegClass <- matrix(nrow=nOuter,
                                       ncol=length(levels(Y)))}
    VALSeg <- matrix(nrow=nOuter,
                     ncol=astep)   ## number of alphas from low or high
    colnames(VALSeg) <- avect
    rownames(VALSeg) <- paste0('Segment',1:nOuter)
    varSeg <- list()

    if(DA) {varSegClass <- list()}
    ## Perform outer loop segments -> one "majority vote" MV model per segment
    for (i in 1:nOuter) {
      # i <- 1
      # i <- i+1
      cat('\n Segment ',i,' (Inner repeat): ',sep='') # Counter
      ## Draw out test set
      testID <- allTest[[i]] # Draw out segment = holdout set BASED ON UNIQUE ID
      testIndex <- ID%in%testID # Boolean for samples included
      xTest <- X[testIndex,]
      yTest <- Y[testIndex]  ## some of the ytest here may only have 1 level
      inID <- unikID[!unikID%in%testID]  # IDs not in test set
      inIndex <- ID%in%inID
      xIn <- X[inIndex,]
      yIn <- Y[inIndex]
      idIn <- ID[inIndex]
      # Allocate variables for inner repetitions
      VALInner <- matrix(nrow=methParam$nRepInner,
                         ncol=astep)
      # Repeat inner CV several times for improved stability

      for (j in 1:methParam$nRepInner) {  ##Is the same as nInners
        cat(j,'... ')
        if (DA & identical(unikID,ID)) {
          inY <- unikY[!unikID%in%testID]  # Counterintuitive, but needed for grouping by Ynames
          allVal <- uniqDASamp(inY, inID, nInner)
        } else {
          allVal <- vectSamp(inID,n=nInner)
        }
        foldID <- foldVector(foldList = allVal,
                             ID = idIn)  ## identify each ID is in which group
        ## Perform steps with different a
        for (count in 1:astep) {  # Build models with successively fewer variables. Quality metric = number of missclassifications for Validation set
          # count <- 1 # for testing
          # count <- count + 1
          alpha <- avect[count]  ## avect is the alpha vector
          penaltyfactor<-rep(1,ncol(xIn))
          if(!is.null(keep)){
          if(any(!keep%in%colnames(X_original))){
            stop("Could not find that variable(s) in X")}
          if(any(!keep%in%colnames(X))){
              stop("Some variables you want to keep are of near zero varance")}
               filter<-which(colnames(xIn)%in% keep)
               penaltyfactor[filter]<-0
          }
####################################### drop levels if the class not appear
################# safe guard measure  first time
          if(is.factor(yIn)){
            if(length(levels(yIn))!=length(levels(droplevels(yIn)))){

              yIn=droplevels(yIn)
            }
            if(any(table(yIn)==1)){

              remove_samp<-yIn==names(table(yIn))[which(table(yIn)==1)]
              xIn<-xIn[-which(remove_samp),]
              yIn<-yIn[-which(remove_samp)]
            }
            if(length(levels(yIn))!=length(levels(droplevels(yIn)))){

              yIn=droplevels(yIn)
            }
          }

          suppressWarnings(inMod <- cv.glmnet(x = xIn,
                                              y = yIn,
                                              alpha=alpha,
                                              penalty.factor=penaltyfactor,

                                              # exclude=filter, positions, this is which the beta coefficient is forced to Zero
                                              family=methParam$family,
                                              foldid = foldID))
          whichMin <- which.min(inMod$cvm)  ## 	The mean cross-validated error - a vector of length length(lambda).
          VALInner[j,count] <- inMod$cvm[whichMin]
        }  ##each astep loop
        # matplot(avect, t(VALInner), log = 'x', ylab='MSE', main=paste0('Rep ',r,' - Segment ',i), type='b')
        VALSeg[i,] <- colMeans(VALInner)
      } ## where each nRepInner end

      alphaSeg[i] <- avect[which.min(VALSeg[i,])]

      penaltyfactor<-rep(1,ncol(xIn))
      if(!is.null(keep)){
        filter<-which(colnames(xIn)%in% keep)
        penaltyfactor[filter]<-0
      }

################################### drop levels if the class not appear
############## safe guard measure  second time
      if(is.factor(yIn)){
        if(length(levels(yIn))!=length(levels(droplevels(yIn)))){

          yIn=droplevels(yIn)
        }
        if(any(table(yIn)==1)){

          remove_samp<-yIn==names(table(yIn))[which(table(yIn)==1)]
          xIn<-xIn[-which(remove_samp),]
          yIn<-yIn[-which(remove_samp)]
        }
        if(length(levels(yIn))!=length(levels(droplevels(yIn)))){

          yIn=droplevels(yIn)
        }
      }


      suppressWarnings(inMod <- cv.glmnet(x = xIn,
                                          y = yIn,
                                          alpha=alphaSeg[i],
                                          penalty.factor=penaltyfactor,
                                          family=methParam$family,
                                          foldid = foldID))

      if (modReturn){ outMod[[i]] <- inMod}

      whichMin <- which.min(inMod$cvm)
      lambdaSeg[i] <- inMod$lambda.min

      if(DA) {
        coefs <- coef(inMod,
                      s = lambdaSeg[i]) %>% sapply(., as.matrix)

        coefs <- coefs[-1,]   #### remove the intercept
###############################################################
################if one column missing add a new column
        coefs_correct<-matrix(0,
                              nrow(coefs),
                              length(levels(Y)))
        colnames(coefs_correct)<-levels(Y)
        for(m in 1:ncol(coefs_correct)){
          if(colnames(coefs_correct)[m]%in%colnames(coefs)){
            for(n in 1:nrow(coefs_correct))
            {coefs_correct[n,m]<-coefs[n,colnames(coefs_correct)[m]]}
          }
        }

        coefs<-coefs_correct
#########################################################################
        rownames(coefs) <- colnames(X)

        coefSeg[,,i] <- coefs
        nonZeroSeg[i] <- apply(coefs,
                               1,
                               function(x) {any(x!=0)}) %>% sum
        nonZeroSegClass[i,] <- apply(coefs,
                                     2,
                                     function(x) sum(x!=0))
        varSegClass[[i]] <- lapply(as.data.frame(coefs),
                                   function(x, nm=rownames(coefs)) {
          whichNZ <- which(x!=0)
          return(nm[whichNZ])
        })

        ######
        varSeg[[i]] <- unlist(varSegClass[[i]]) %>% unique
        varSeg[[i]] <- varSeg[[i]][order(match(varSeg[[i]],
                                               colnames(X)))] # Resort according to order in X
        yPredSeg[testIndex,] <- predict(inMod,
                                        newx = xTest,
                                        s = lambdaSeg[i],
                                        type = 'response')
      } else {
        nonZeroSeg[i] <- inMod$nzero[whichMin]
        coefSeg[,i] <- as.matrix(coef(inMod,
                                      s=lambdaSeg[i]))[-1,,drop=F]
        ## remove the intercept

        varSeg[[i]] <- names(coefSeg[,i][coefSeg[,i]!=0])
        yPredSeg[testIndex] <- predict(inMod,
                                       newx = xTest,
                                       s = lambdaSeg[i])
      }
      # plot(yTest, yPred)
    }  ## where the outer loop ends

    # matplot(avect, t(VALSeg), log = 'x', ylab='MSE', main=paste0('Rep ',r,' - All segments'), type='b')
    # Per repetition: Average outer loop predictions, VIP ranks and nComp for PLS
    varSegClassReorder <- list()
    if(DA){
    for (i in 1:nOuter) {
      for (j in 1:length(levels(Y))) {
        varSegClassReorder[[levels(Y)[j]]][[i]] <- varSegClass[[i]][levels(Y)[j]]
      }
    }
    }
    parReturn <- list(yPred = yPredSeg,
                      alpha = alphaSeg,
                      lambda = lambdaSeg,
                      nonZero = nonZeroSeg,
                      coef = coefSeg,
                      VAL = VALSeg,
                      vars = varSeg)
    if(DA) {
      parReturn$varsClass <- varSegClassReorder
      parReturn$nonZeroClass <- nonZeroSegClass
    }
    # Return underlying models
    if (modReturn) parReturn$outModel=outMod
    return(parReturn)
    # reps[[r]] <- parReturn
  }   ## Repetition loop end

  ####################################
  # Allocate prediction variables
  if (DA) {
    # Predictions for all repetitions
    yPredRep <- array(dim=c(length(Y),
                            length(levels(Y)),
                            nRep),
                      dimnames=list(ID,
                                    levels(Y),
                                    paste('Rep', 1:nRep,sep='')))
  } else {
    # Predictions for all repetitions
    yPredRep <- matrix(nrow=length(Y),
                       ncol=nRep,
                       dimnames=list(ID, paste('Rep', 1:nRep,sep='')))
  }

  # Allocate validation and parameter variables
  alphaRep <- lambdaRep <- nonZeroRep <-
    matrix(nrow=nRep,
           ncol=nOuter,
           dimnames=list(paste('repetition', 1:nRep, sep=''),
                         paste('segment', 1:nOuter,sep='')))
  if(DA) {coefRep=list()} else {
    coefRep <- array(0,
                     dim=c(nVar0, nOuter, nRep),
                     dimnames=list(colnames(X),
                                   paste0('Segment',1:nOuter),
                                   paste0('Rep',1:nRep)))}



  VALRep <- array(dim=c(nOuter, astep, nRep),
                   dimnames=list(paste('outSeg',
                                       1:nOuter,paste=''),
                                 avect,
                                 paste(rep('rep', nRep),
                                              1:nRep,sep='')))
  varRep <- list()
  if(DA) {
    varsClassRep <- list()
    nonZeroClassRep <- array(dim=c(nOuter,
                                   length(levels(Y)),
                                   nRep),
                             dimnames=list(paste('outSeg', 1:nOuter,paste=''),
                                           levels(Y),
                                           paste(rep('rep', nRep),
                                                 1:nRep,sep='')))
  }

  # Allocate variable for models (for later predictions)
  if (modReturn) {outMods=list()}

  ####################################
  # Aggregate results over predictions
  for (r in 1:nRep) {
    if (DA) {yPredRep[,,r] <- reps[[r]]$yPred
    }else {
      yPredRep[,r] <- reps[[r]]$yPred}
    alphaRep[r,] <- reps[[r]]$alpha
    lambdaRep[r,] <- reps[[r]]$lambda
    nonZeroRep[r,] <- reps[[r]]$nonZero

    if(DA) {
      for (i in 1:length(levels(Y))) {
        varsClassRep[[levels(Y)[i]]][[r]] <- reps[[r]]$varsClass[[levels(Y)[i]]]
      }
      nonZeroClassRep[,,r] <- reps[[r]]$nonZeroClass
    }

    if(DA){ coefRep[[r]] <- reps[[r]]$coef
    }else {
      coefRep[,,r] <- reps[[r]]$coef}

    VALRep[,,r] <- reps[[r]]$VAL
    varRep <- c(varRep,reps[[r]]$vars)

    if (modReturn){
      outMods <- c(outMods,reps[[r]]$outModel)}
  }


  # Average predictions
  if (DA) {
    yPred <- apply(yPredRep,
                   c(1,2),  ## what dimension is not changed by the function
                   mean)
  } else {
    yPred <- apply(yPredRep,
                   1,
                   mean)[1:nrow(X)]
  }
  # matplot(Y,yPredRep, pch=16, col='lightgrey', cex=0.7)
  # points(Y, yPred, pch=16, col='grey20',cex=.9)
  modelReturn$yPred <- yPred

  if (DA) {
    # Classify predictions
    yClass <- apply(yPred,
                    1,
                    function(x) {levels(Y)[which.max(x)]})
    miss <- sum(yClass!=Y)
    auc <- numeric(length(levels(Y)))
    ber<-getBER(actual=Y,
                  predicted=yClass)
    names(auc) <- levels(Y)
    for (cl in 1:length(levels(Y))) {auc[cl] <- roc(Y==(levels(Y)[cl]),
                                                    yPred[,cl],
                                                    quiet = TRUE)$auc}
    # # Report
    modelReturn$yClass <- yClass
    modelReturn$miss <- miss
    modelReturn$ber<-ber
    modelReturn$auc <- auc
  } else if (ML) {
    yClass <- ifelse(yPred>0,1,-1)
    modelReturn$yClass <- yClass
    modelReturn$miss <- sum(yClass!=Y)
    modelReturn$auc <- roc(Y,yPred, quiet=TRUE)$auc
    modelReturn$ber<-getBER(actual=Y,
                            predicted= yClass )
    names(modelReturn$yClass)=paste(1:nSamp,
                                    ID,
                                    sep='_ID')
  }

  # Calculate overall nVar
  VAL <- apply(VALRep,
               2,
               mean)  ### averaged all Repetition and Outer loop to get the overall error
  modelReturn$varTable <- varRep %>% unlist %>% table %>% sort(., decreasing = T)
  modelReturn$alpha <- avect[which.min(VAL)]
  modelReturn$alphaRep <- alphaRep
  modelReturn$lambdaRep <- lambdaRep
  modelReturn$coefRep <- coefRep
  modelReturn$nonZeroRep <- nonZeroRep
  modelReturn$varRep <- varRep
  modelReturn$VAL$metric <- fitness
  modelReturn$VAL$VALRep <- VALRep
  modelReturn$VAL$VAL <- VAL
  modelReturn$VAL$avect <- avect
  if (DA) {
    modelReturn$varClassRep <- sapply(varsClassRep,
                                      function(x) {unlist(x) %>% table %>% sort(., decreasing = T)})
    modelReturn$nonZeroClassRep <- nonZeroClassRep
  }
  if (modReturn) {modelReturn$outModels=outMods}

  modelReturn$yPredPerRep <- yPredRep
  modelReturn$inData=InData

#########################################################
  # Fit-predict model
  penaltyfactor<-rep(1,ncol(X))
  if(!is.null(keep)){
    filter<-which(colnames(X)%in% keep)
    penaltyfactor[filter]<-0
  }

  ####################################### drop levels if the class not appear
  ################# safe guard measure  first time
  if(is.factor(Y)){
    if(length(levels(Y))!=length(levels(droplevels(Y)))){

      Y=droplevels(Y)
    }
    if(any(table(Y)==1)){

      remove_samp<-Y==names(table(Y))[which(table(Y)==1)]
      xIn<-xIn[-which(remove_samp),]
      Y<-Y[-which(remove_samp)]
    }
    if(length(levels(Y))!=length(levels(droplevels(Y)))){

      Y=droplevels(Y)
    }
  }

  ### Use the overall alpha to get the best lambda
  suppressWarnings(cvFit <- cv.glmnet(X,
                                      Y,
                                      penalty.factor=penaltyfactor,
                                      alpha=modelReturn$alpha,
                                      family = methParam$family))
  ### do prediction models
  suppressWarnings(fitPredict <- glmnet(X,
                                        Y,
                                        penalty.factor=penaltyfactor,
                                        alpha=modelReturn$alpha,
                                        lambda = cvFit$lambda.min,
                                        family = methParam$family))
  yFit <- predict(fitPredict,
                  newx = X,
                  type='response')
  # Calculate fit statistics
  if (!DA) {      ## which means for both regression and multilevel
    TSS <- sum((Y-mean(Y))^2)
    RSS <- sum((Y-yFit)^2)
    PRESS=sum((Y-yPred)^2)
    R2=1-(RSS/TSS)
    Q2=1-(PRESS/TSS)
    modelReturn$fitMetric <- data.frame(R2,Q2)
  }


  # Stop timer
  end.time=proc.time()[3]
  modelReturn$calcMins=(end.time-start.time)/60
  cat('\n Elapsed time',(end.time-start.time)/60,'mins \n')
  class(modelReturn) <- c('MUVR',
                          ifelse(DA,
                                 'Classification',
                                 ifelse(ML,'Multilevel','Regression')))
  return(modelReturn)
}
