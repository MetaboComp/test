#' MUVR: "Multivariate modelling with Unbiased Variable selection"
#'
#' Repeated double cross validation with tuning of variables in the inner loop.
#'
#' @param X Predictor variables. NB: Variables (columns) must have names/unique identifiers. NAs not allowed in data. For multilevel, only the positive half of the difference matrix is specified.
#' @param Y Response vector (Dependent variable). For classification, a factor (or character) variable should be used. For multilevel, Y is calculated automatically.
#' @param ID Subject identifier (for sampling by subject; Assumption of independence if not specified)
#' @param scale If TRUE, the predictor variable matrix is scaled to unit variance for PLS modelling.
#' @param nRep Number of repetitions of double CV. (Defaults to 5)
#' @param nOuter Number of outer CV loop segments. (Defaults to 6)
#' @param nInner Number of inner CV loop segments. (Defaults to nOuter-1)
#' @param varRatio Ratio of variables to include in subsequent inner loop iteration. (Defaults to 0.75)
#' @param DA Boolean for Classification (discriminant analysis) (By default, if Y is numeric -> DA = FALSE. If Y is factor (or character) -> DA = TRUE)
#' @param fitness Fitness function for model tuning (choose either 'AUROC' or 'MISS' (default) for classification; or 'RMSEP' (default) for regression.)
#' @param method Multivariate method. Supports 'PLS' and 'RF' (default)
#' @param methParam List with parameter settings for specified MV method (see function code for details)
#' @param ML Boolean for multilevel analysis (defaults to FALSE)
#' @param modReturn Boolean for returning outer segment models (defaults to FALSE). Setting modReturn = TRUE is required for making MUVR predictions using predMV().
#' @param logg Boolean for whether to sink model progressions to `log.txt`
#' @param parallel Boolean for whether to perform `foreach` parallel processing (Requires a registered parallel backend; Defaults to `TRUE`)
#'
#' @return
#' A MUVR object
#' @export
MUVR <- function(X,
                 Y,
                 ID,
                 scale = TRUE,
                 nRep = 5,
                 nOuter = 6,
                 nInner,
                 varRatio = 0.75,
                 DA = FALSE,
                 fitness = c('AUROC', 'MISS', 'BER', 'RMSEP'),
                 method = c('PLS',' RF'),
                 methParam,
                 ML = FALSE,
                 modReturn = FALSE,
                 logg=FALSE,
                 parallel = TRUE,
                 ...) {

  # Start timer
  start.time <- proc.time()[3]

  # Initialise modelReturn with function call
  modelReturn <- list(call = match.call())

  # Default core modelling method
  if (missing(method)) method <- 'RF'

  # Call in relevant package(s)
  library(pROC)
  library(foreach)
  if (method == 'RF') library(randomForest)

  # Parallel processing
  if (parallel) "%doVersion%" <- get("%dopar%") else "%doVersion%" <- get("%do%")

  # Rough check indata
  if (length(dim(X)) != 2) stop('\nWrong format of X matrix.\n')
  if (is.null(colnames(X))) stop('\nNo column names in X matrix.\n')
  X <- as.matrix(X) # PROBLEM: Will not work for factor variables BUT it will work for PLS with one-hot
  # I don't know how much of a problem data frames are from a time perspective. Check matrix vs DF on same dataset for time difference.






  ###########################################################################
  #To test the scenario when X has factor and charactor when using PLS
  #add one factor and one character variable(freelive data X, which originally has 112 numeric samples and 1147 observations)
  #factor varaible has 3 factors(nearzero varianece),character variable has 7 categories
  # factor_variable1<-as.factor(c(rep(33,105),rep(44,3),rep(55,4)))
  #factor_variable2<-as.factor(c(rep("AB",20),rep("CD",10),rep("EF",30),
  #                           rep("GH",15),rep("IJ",25),rep("KL",12)))
  #factor_variable3<-as.factor(c(rep("Tessa",25),rep("Olle",30),rep("Yan",12),
  #                            rep("Calle",25),rep("Elisa",20)))
  # character_variable1<-c(rep("one",16),rep("two",16),rep("three",16),
  #                       rep("four",16),rep("five",16),rep("six",16),rep("seven",16))
  #  character_variable2<-c(rep("yes",28),rep("no",28),
  #                          rep("yes",28),rep("no",28))
  #  X=XRVIP
  #  X<-as.data.frame(X)
  #  X<-cbind(X,
  #        factor_variable1,factor_variable2,factor_variable3,
  #        character_variable1,character_variable2)
  #   check again by using class(X[,1148])
  #

  #################################################################################
  # One-hot expansion of factor variables for PLS,these are the code to add in the model
  #This means that X need to checked manually to make sure each variable class is correct
  #filter the variables that are factors and character first
  X<-as.data.frame(X)
  if (method == "PLS") {
    cat("This is PLS. All variables are transformed to numeric")
    #find factor and character
    X[,which(sapply(X, class) %in% c('factor'))]
    X[,which(sapply(X, class) %in% c('character'))]
    #store names
    X_factor_names<-colnames(X)[which(sapply(X, class) %in% c('factor'))]
    X_character_names<-colnames(X)[which(sapply(X, class) %in% c('character'))]
    #transform character in to factor a

    X_character_frame<-as.data.frame(X[,which(sapply(X, class) %in% c('character'))])

    if(ncol(X_character_frame)!=0)
    {for (c in 1:ncol(X_character_frame))
    {X_character_frame[,c]<-as.factor(X_character_frame[,c])}
    }
    colnames(X_character_frame)<-X_character_names
    ##put factor and previous character in to one data frame and game them previous names
    X_factor_frame<-cbind(as.data.frame(X[,which(sapply(X, class) %in% c('factor'))]),X_character_frame)

    colnames(X_factor_frame)<-c(X_factor_names,X_character_names)
    ##orginal numeric frame
    X_numeric_frame<-as.data.frame(X[,which(sapply(X, class) %in% c('numeric'))])

    ##test how many levels a factor have, if=0 >5}
    if(ncol(X_factor_frame)==0)
    {paste("all variables are numeric")

    }else{


      paste(ncol(X_factor_frame),"non-numeric variables")
      ###a is the list of the names of vairables to be built
      a<-list()
      for(n in 1:ncol(X_factor_frame))
      {if(length(levels(X_factor_frame[,n]))>5)
      {paste(colnames(X_factor_frame)[n],"has",length(levels(X_factor_frame[,n])),"levels, which is too many")}
        a[[n]]<-character()

        for (m in 1:length(levels(X_factor_frame[,n])))
        {a[[n]][m]<-paste0(colnames(X_factor_frame)[n],"_","level","_",
                           ##NOT sort the level by alphabet order
                           levels(factor(X_factor_frame[,n], as.character(unique(X_factor_frame[,n]))))[m])
        }

      }
      ##b is the list of 0/1 matrix of each variable
      b<-list()
      c<-list()
      for(n in 1:ncol(X_factor_frame))
      {b[[n]]<-data.frame(row.names=c(1:nrow(X_factor_frame)))
      c[[n]]<-data.frame(row.names=c(1:nrow(X_factor_frame)))
      for(i in 1:length(levels(X_factor_frame[,n])))
      {b[[n]]<-cbind(b[[n]],X_factor_frame[,n])
      c[[n]]<-cbind(b[[n]],X_factor_frame[,n])
      }
      c[[n]]<-matrix(data=0,
                     nrow=nrow(X_factor_frame),
                     ncol=length(levels(X_factor_frame[,n])))

      for(m in 1:length(levels(X_factor_frame[,n])))
      { for(z in 1:nrow(X_factor_frame))
       {if(b[[n]][z,m]==as.factor(levels(factor(X_factor_frame[,n], as.character(unique(X_factor_frame[,n]))))[m]))
        {c[[n]][z,m]=1
        }
       }
      }
      #Now we can summarize a bit
      # a is the list of new name of all the new variables that will be used
      # b is an intermediate step to use to check for error
      # c is the new 0-1 matrix
      # What needs to be done now is to combine a and c

      rownames(c[[n]])<-rownames(X)
      colnames(c[[n]])<-a[[n]]

      }
      #Now I need to combine the c matrix with original X dataset X_numeric_frame
      new_X_frame<-X_numeric_frame
      for(h in 1:length(c)){new_X_frame<-cbind(new_X_frame,c[[h]])}

        }  ##This is where else end

}   ##this is where method=="PLS" ends
##X_numeric_frame is the data frame that only has original numeirc variable
##new_X_frame is the result of one hot coding, it is a dataframe though
  new_X_matrix<-as.matrix(new_X_frame)
  View(new_X_matrix)
##Now new_X_frame becomes a matrix that contains original numeric variables and transformed non-numeric variable into one hot coding



#########################################################################################################
########################################################################################################
# Remove nearZeroVariance variables for PLS
#  X<-new_X_matrix  This could be used if we decide to use X for all the following code

  if (method == 'PLS') {
    nzv <- MUVR::nearZeroVar(new_X_matrix) # Function borrowed from mixOmics
    if (length(nzv$Position) > 0) {
      modelReturn$nzv <- colnames(new_X_matrix)[nzv$Position]
      new_X_matrix <- new_X_matrix[, -nzv$Position]
      cat('\n',length(nzv$Position),'variables with near zero variance detected -> removed from X and stored under $nzv')
    }
  }


############################################################################################################
############################################################################################################



  # Number of samples and variables
  nSamp <- nrow(X)
  nVar <- nVar0 <- ncol(X)

  # Sample identifiers; Assume independence if ID is not specified
  if (missing(ID)) {
    cat('\nMissing ID -> Assume all unique (i.e. sample independence)')
    ID <- 1:nSamp
  }

  # Figure out number of iterations in inner CV loop
  # Gets tweaked for "keeps"
  var <- numeric()
  cnt <- 0
  while (nVar > 1) {
    cnt <- cnt + 1
    var <- c(var, nVar)
    nVar <- floor(varRatio * nVar)
  }

  # Set up default internal modelling parameters
  if (missing(nInner)) nInner <- nOuter - 1 # Default value for inner segments

  # methParams
  if(any(names(list(...)) == 'nCompMax')) stop('`nCompMax` is deprecated. Use customParams() and the methParam argument in MUVR instead.')
  if (missing(methParam)) {
    if (method == 'PLS') {
      methParam <- list(compMax = 5)
      if (nVar < methParam$compMax) methParam$compMax <- nVar # nCompMax cannot be larger than number of variables!
    } else { # i.e. for RF
      methParam <- list(ntreeIn = 150,
                        ntreeOut = 300,
                        mtryMaxIn = 150)
    }
    methParam$robust <- 0.05
  }

  # Set up for multilevel analysis
  if (ML) {
    X <- rbind(X, -X)
    if (missing(Y)) Y <- rep(-1, nSamp)
    Y <- c(Y, -Y)
    nSamp <- 2 * nSamp
    ID <- c(ID, ID)
    DA <- FALSE
    fitness <- 'MISS'
    cat('\nMultilevel -> Regression on (-1,1) & fitness = MISS')
  }

  # No missingness allowed - Please perform imputations before running MUVR
  if (any(is.na(X)) | any(is.na(Y))) stop('\nNo missing values allowed in X or Y data.\n')
  if (!is.null(dim(Y))) {
    cat('\nY is not a vector: Return NULL')
    return(NULL)
  }

  # DA / Classification
  if (is.character(Y)) Y <- factor(Y)
  if (is.factor(Y)) {
    cat('\nY is factor -> Classification (',length(unique(Y)),' classes)', sep = '')
    DA <- TRUE
  }
  if (is.numeric(Y) & DA) {
    Y <- as.factor(Y)
    cat('\nDA = TRUE -> Y as factor -> Classification (',length(unique(Y)),' classes)', sep = '')
  }

  # Check fitness criterion
  if (missing(fitness)) {
    if (DA) {
      fitness = 'MISS'
      cat('\nMissing fitness -> MISS')
    } else { # I.e. for regression
      fitness = 'RMSEP'
      cat('\nMissing fitness -> RMSEP')
    }
  }

  # Additional sanity check
  if (nrow(X) != length(Y)) stop('\nMust have same nSamp in X and Y: Return NULL')

  # Store indata in list for later model return
  InData <- list(X = X,
                 Y = Y,
                 ID = ID,
                 scale = scale,
                 nRep = nRep,
                 nOuter = nOuter,
                 nInner = nInner,
                 varRatio = varRatio,
                 DA = DA,
                 fitness = fitness,
                 method = method,
                 methParam = methParam,
                 ML = ML,
                 parallel = parallel)

  # Sort sampling based on ID and not index & Allocate prediction
  unik <- !duplicated(ID)  # boolean of unique IDs
  unikID <- ID[unik]       # Actual unique IDs
  if (DA) {
    if(nOuter > min(table(Y))) warning('\nnOuter is larger than your smallest group size. Consider lowering your nOuter to min(table(Y)).',call.=TRUE)
    unikY <- Y[unik]  # Counterintuitive, but needed for groupings by Ynames
    Ynames <- sort(unique(Y))  # Find groups
    groups <- length(Ynames) # Number of groups
    groupID <- list()  # Allocate list for indices of groups
    for (g in 1:groups) {
      groupID[[g]] <- unikID[unikY == Ynames[g]]  # Find indices per group
    }
    # Allocate final predictions for min mid and max models
    yPredMin <- yPredMid <- yPredMax <- array(dim = c(length(Y), # Rows = number of observations
                                                      length(levels(Y)), # Columns = number of classes in Y
                                                      nRep), # 3rd = number of repetitions
                                              dimnames = list(ID, # Observation ID names
                                                              levels(Y), # Names of levels in Y
                                                              paste('Rep', 1:nRep, sep = ''))) # Repetitions
    # Allocate predictions per repetition
    yPredMinR <- yPredMidR <- yPredMaxR <- matrix(nrow = length(Y), # Like above but lacking 3rd dimension (repetitions)
                                                  ncol = length(levels(Y)),
                                                  dimnames = list(ID, levels(Y)))
  } else { # I.e. for regression and ML
    # Allocate final predictions for min mid and max models
    yPredMin <- yPredMid <- yPredMax <- matrix(nrow = length(Y), # Like above but matrix instead of array, since there are not multiple classes
                                               ncol = nRep,
                                               dimnames = list(ID, paste('Rep', 1:nRep, sep = '')))
    # Allocate predictions per repetition
    yPredMinR <- yPredMidR <- yPredMaxR <- numeric(length(Y)) # Like above but lacking columns (repetitions) -> numeric vector
  }

  # Allocate response vectors and matrices for var's, nComp and VIP ranks over repetitions
  missRep <- numeric(nRep)
  names(missRep) <- paste(rep('rep',nRep), 1:nRep, sep = '')
  varRepMin <- varRepMid <- varRepMax <- nCompRepMin <- nCompRepMid <- nCompRepMax <- missRep
  nCompSegMin <- nCompSegMid <- nCompSegMax <- matrix(nrow = nRep,
                                                      ncol = nOuter,
                                                      dimnames = list(paste('repetition', 1:nRep, sep = ''),
                                                                      paste('segment', 1:nOuter, sep='')))
  VIPRepMin <- VIPRepMid <- VIPRepMax <- matrix(data = nVar0,
                                                nrow = nVar0,
                                                ncol = nRep,
                                                dimnames = list(colnames(X),
                                                                paste(rep('rep', nRep), 1:nRep, sep = '')))

  # Allocate array for validation results
  VAL <- array(dim = c(nOuter, cnt, nRep),
               dimnames = list(paste('outSeg', 1:nOuter, paste=''),
                               var,
                               paste(rep('rep',nRep),1:nRep,sep='')))

  # Choose package/core algorithm according to chosen method
  # And prepare for exporting them in 'foreach' (below)
  packs <- c('pROC')
  if(method == 'RF') packs <- c(packs, 'randomForest')
  exports <- 'vectSamp'


  ####################################################
  ## Start repetitions
  ####################################################

  # For manual debugging purposes
  # reps=list()
  # for (r in 1:nRep){

  reps <- foreach(r = 1:nRep, .packages = packs, .export = exports) %doVersion% {

    # For manual debugging purposes
    # r <- 1
    # r <- r + 1

    # Send intermediate outputs to log file: Not a pretty output. For debugging purposes only.
    if (logg) sink('log.txt', append = TRUE)
    # Intermediate info output
    cat('\n','   Repetition ',r,' of ',nRep,':',sep='')

    # Allocate output for models (i.e. for later prediction of external samples)
    if (modReturn) outMod <- list()

    # PERFORM SEGMENTATION INTO OUTER SEGMENTS
    if (DA & identical(unikID, ID)) {
      groupTest <- list()  ## Allocate list for samples within group
      for (gT in 1:groups) {
        groupTest[[gT]] <- vectSamp(groupID[[gT]], n = nOuter)  # Draw random samples within group
      }
      allTest <- groupTest[[1]] # Add 1st groups to 'Master' sample of all groups
      for (gT in 2:groups) {  # Add subsequent groups
        allTest <- allTest[order(sapply(allTest, length))]
        for (aT in 1:nOuter) {
          allTest[[aT]] <- sort(c(allTest[[aT]], groupTest[[gT]][[aT]]))
        }
      }
    } else {
      allTest <- vectSamp(unikID, n=nOuter)
    }

    # Allocate intermediate output objects
    nCompOutMax <- numeric(nOuter)
    names(nCompOutMax) <- paste(rep('outSeg', nOuter), 1:nOuter, sep = '')
    varOutMin <- varOutMid <- varOutMax <- nCompOutMin <- nCompOutMid <- nCompOutMax
    VIPOutMin <- VIPOutMid <- VIPOutMax <- matrix(data = nVar0,
                                                  nrow = nVar0,
                                                  ncol = nOuter,
                                                  dimnames = list(colnames(X),
                                                                  paste(rep('outSeg', nOuter), 1:nOuter, sep = '')))
    VALRep <- matrix(nrow = nOuter, ncol = cnt)

    # Perform outer loop segments -> one "majority vote" MV model per segment
    for (i in 1:nOuter) {

      # For manual debugging purposes
      # i <- 1
      # i <- i + 1

      # Intermediate info output
      cat('\n Segment ',i,' (variables):',sep='') # Counter

      # Draw out test set
      testID <- allTest[[i]] # Draw out segment = holdout set BASED ON UNIQUE ID
      testIndex <- ID%in%testID # Boolean for samples corresponding to unique test IDs
      xTest=X[testIndex,]
      yTest=Y[testIndex]

      # Inner data (not in test)
      inID <- unikID[!unikID%in%testID]  # IDs not in test set
      if (DA & identical(unikID,ID)) inY <- unikY[!unikID%in%testID]  # Counterintuitive, but needed for grouping by Ynames

      # Allocate fitness variables for the inner data
      missIn <- berIn <- aucIn <- rmsepIn <- PRESSIn <- nCompIn <- matrix(nrow = nInner,
                                                                          ncol = cnt,
                                                                          dimnames = list(paste(rep('inSeg', nInner), 1:nInner, sep = ''),
                                                                                          var))
      # Allocate VIP variable for the inner data
      VIPInner <- array(data = nVar0,
                        dim = c(nVar0,cnt,nInner),
                        dimnames = list(colnames(X),
                                        var,
                                        paste(rep('inSeg', nInner), 1:nInner, sep = '')))
      # Set initial inclusion list of variables to all variables
      incVar <- colnames(X)

      # Perform steps with successively fewer variables
      for (count in 1:cnt) {  # Build models with successively fewer variables.

        # For manual debugging purposes
        # count <- 1 # for testing
        # count <- count + 1

        # Intermediate info output
        cat(nVar)

        # Extract the number of variables at the present count (according to the count loop before the foreach loop)
        nVar <- var[count]

        # Tweak method parameters for low number of variables
        if (method == 'PLS') comp <- min(c(nVar, methParam$compMax)) # nComp cannot be > nVar
        if (method=='RF') {
          mtryIn <- ifelse(DA,
                           min(c(methParam$mtryMaxIn,floor(sqrt(nVar)))), # Standard but with upper limit
                           min(c(methParam$mtryMaxIn,floor(nVar/3)))) # Standard but with upper limit
          mtryIn <- max(c(2,mtryIn)) # Lower limit
        }

        # PERFORM SEGMENTATION INTO INNER SEGMENTS
        if (DA & identical(unikID,ID)) {
          groupIDVal <- list()
          for (g in 1:groups) {
            groupIDVal[[g]] <- inID[inY == Ynames[g]]  # Find indices per group
          }
          groupVal <- list()  ## Allocate list for samples within group
          for (gV in 1:groups) {
            groupVal[[gV]] <- vectSamp(groupIDVal[[gV]], n = nInner)  # Draw random samples within group
          }
          allVal = groupVal[[1]] # Add 1st groups to 'Master' sample of all groups
          for (gV in 2:groups) {  # Add subsequent groups
            allVal = allVal[order(sapply(allVal, length))]
            for (aV in 1:nInner) {
              allVal[[aV]] <- sort(c(allVal[[aV]], groupVal[[gV]][[aV]]))
            }
          }
        } else {
          allVal <- vectSamp(inID, n = nInner)
        }

        # Perform inner CV loop
        for (j in 1:nInner) {

          # For manual debugging purposes
          # j <- 1
          # j <- j + 1

          # Intermediate info output
          cat('.') # Counter

          # Extract validation segment
          valID <- allVal[[j]] # Extract IDs
          valIndex <- ID %in% valID # Corresponding observations
          xVal <- X[valIndex,]
          xVal <- subset(xVal, select = incVar) # limit to the current selection of variables (selection below)
          yVal <- Y[valIndex]

          # Extract training data
          trainID <- inID[!inID %in% valID] # The inner but not validation IDs
          trainIndex <- ID %in% trainID # Corresponding observations
          xTrain <- X[trainIndex,]
          xTrain <- subset(xTrain, select = incVar) # limit to the current selection of variables (selection below)
          yTrain <- Y[trainIndex]

          # For debugging
          # sum(trainIndex,valIndex,testIndex)
          # trainIndex|valIndex|testIndex

          # Make inner model
          if (method == 'PLS') {
            inMod <- MUVR::plsInner(xTrain,
                                    yTrain,
                                    xVal,
                                    yVal,
                                    DA,
                                    fitness,
                                    comp,
                                    scale = scale)
            nCompIn[j, count] <- inMod$nComp
          } else {
            inMod <- MUVR::rfInner(xTrain,
                                   yTrain,
                                   xVal,
                                   yVal,
                                   DA,
                                   fitness,
                                   ntree = methParam$ntreeIn,
                                   mtry = mtryIn)
          }

          # Store fitness metric
          if (fitness == 'MISS') {
            missIn[j, count] <- inMod$miss
          } else if (fitness == 'BER') {
            berIn[j, count] <- inMod$ber
          } else if (fitness == 'AUROC') {
            aucIn[j, count] <- inMod$auc
          } else {
            rmsepIn[j, count] <- inMod$rmsep
            PRESSIn[j, count] <- (inMod$rmsep ^ 2) * length(yVal) # Is this really correct???
          }

          # Store VIs
          VIPInner[match(names(inMod$vi), rownames(VIPInner)), count, j] <- inMod$vi
        }

        # Average inner VIP ranks before variable elimination - Tweak for `keeps`
        VIPInAve <- apply(VIPInner[, count, ], 1, mean)
        if (count < cnt) {
          incVar <- names(VIPInAve[order(VIPInAve)])[1:var[count + 1]] # Extract the names of the variables kept for the next iteration
        }
      }

      # PER OUTER SEGMENT:

      # Store fitness curves and construct fitness rank curves for the different "counts"
      if (fitness == 'AUROC') {
        fitRank <- colMeans(-aucIn)
        VALRep[i,] <- colMeans(aucIn)
      } else if (fitness == 'MISS') {
        fitRank <- VALRep[i,] <- colSums(missIn)
      } else if (fitness == 'BER') {
        fitRank <- VALRep[i,] <- colMeans(berIn)
      } else {
        fitRank <- colMeans(rmsepIn)
        VALRep[i,] <- sqrt(colSums(PRESSIn) / sum(!testIndex))
      }
      # Rescale fitRank to range 0 (best) - 1 (worst)
      fitRank=(fitRank-min(fitRank))/abs(diff(range(fitRank)))
      # BugCatch: If all VAL have NaN value -> reset all fitRank to 0
      if(all(is.nan(fitRank))) fitRank <- rep(0,cnt)

      # Extract index of min and max models (remember that counts and nVar go in opposite direction!)
      minIndex <- max(which(fitRank <= methParam$robust))
      maxIndex <- min(which(fitRank <= methParam$robust))

      # nVar at min, mid and max
      varOutMin[i] <- var[minIndex]
      varOutMax[i] <- var[maxIndex]
      varOutMid[i] <- round(exp(mean(log(c(var[minIndex], var[maxIndex]))))) # Geometric mean of min and max

      # midIndex is closest to nVarMid
      midIndex <- which.min(abs(var - varOutMid[i]))

      # For PLS, extract number of components
      if (method == 'PLS') {
        nCompOutMin[i] <- round(mean(nCompIn[, minIndex]))
        nCompOutMid[i] <- round(mean(nCompIn[, midIndex]))
        nCompOutMax[i] <- round(mean(nCompIn[, maxIndex]))
      }

      # Average VIP ranks
      VIPOutMin[, i] <- apply(VIPInner[, minIndex,], 1, mean)
      VIPOutMid[, i] <- apply(VIPInner[, midIndex,], 1, mean)
      VIPOutMax[, i] <- apply(VIPInner[, maxIndex,], 1, mean)

      # Build outer model for min, mid and max and predict YTEST
      # Extract all inner data
      xIn <- X[!testIndex, ] # Perform Validation on all samples except holdout set
      yIn <- Y[!testIndex]

      # Determine consensus choice of included variables for min mid and max
      incVarMin <- rownames(VIPOutMin)[rank(VIPOutMin[,i]) <= varOutMin[i]]
      incVarMid <- rownames(VIPOutMid)[rank(VIPOutMid[,i]) <= varOutMid[i]]
      incVarMax <- rownames(VIPOutMax)[rank(VIPOutMax[,i]) <= varOutMax[i]]

      # Consensus models for PLS
      if (method == 'PLS'){
        # Build min model
        if (DA) {
          plsOutMin <- MUVR::plsda(subset(xIn, select = incVarMin),
                                   yIn,
                                   ncomp = nCompOutMin[i],
                                   near.zero.var = TRUE,
                                   scale=scale)
        } else {
          plsOutMin <- MUVR::pls(subset(xIn,select = incVarMin),
                                 yIn,
                                 ncomp = nCompOutMin[i],
                                 near.zero.var = TRUE,
                                 scale = scale)
        }
        # Extract test data with correct variable selection
        xTestMin <- subset(xTest, select = incVarMin)
        # Extract predictions
        yPredMinR[testIndex] <- predict(plsOutMin,
                                        newdata = xTestMin,
                                        scale = scale)$predict[, , nCompOutMin[i]]  #
        # Mid model - Update similar to "Build min model" above
        if (DA) plsOutMid=MUVR::plsda(subset(xIn,select=incVarMid),yIn,ncomp=nCompOutMid[i],near.zero.var=TRUE,scale=scale) else
          plsOutMid=MUVR::pls(subset(xIn,select=incVarMid),yIn,ncomp=nCompOutMid[i],near.zero.var=TRUE,scale=scale)
        xTestMid=subset(xTest,select=incVarMid)
        yPredMidR[testIndex]=predict(plsOutMid,newdata=xTestMid,scale=scale)$predict[,,nCompOutMid[i]]  #
        # Max model - Update similar to "Build min model" above
        if (DA) plsOutMax=MUVR::plsda(subset(xIn,select=incVarMax),yIn,ncomp=nCompOutMax[i],near.zero.var=TRUE,scale=scale) else
          plsOutMax=MUVR::pls(subset(xIn,select=incVarMax),yIn,ncomp=nCompOutMax[i],near.zero.var=TRUE,scale=scale)
        xTestMax=subset(xTest,select=incVarMax)
        yPredMaxR[testIndex]=predict(plsOutMax,newdata=xTestMax,scale=scale)$predict[,,nCompOutMax[i]]  #
        if (modReturn) {
          outMod[[i]]=list(plsOutMin,plsOutMid,plsOutMax)
        }
      } else { # This is for RF - needs to be expanded for other cores; expand as for PLS above
        rfOutMin=randomForest(x=subset(xIn,select=incVarMin),y=yIn,xtest=subset(xTest,select=incVarMin),ytest=yTest,ntree=methParam$ntreeOut,keep.forest=TRUE)
        if (DA) {
          yPredMinR[testIndex,]=rfOutMin$test$votes
        } else {
          yPredMinR[testIndex]=rfOutMin$test$predicted
        }
        rfOutMid=randomForest(x=subset(xIn,select=incVarMid),y=yIn,xtest=subset(xTest,select=incVarMid),ytest=yTest,ntree=methParam$ntreeOut,keep.forest=TRUE)
        if (DA) {
          yPredMidR[testIndex,]=rfOutMid$test$votes
        } else {
          yPredMidR[testIndex]=rfOutMid$test$predicted
        }
        rfOutMax=randomForest(x=subset(xIn,select=incVarMax),y=yIn,xtest=subset(xTest,select=incVarMax),ytest=yTest,ntree=methParam$ntreeOut,keep.forest=TRUE)
        if (DA) {
          yPredMaxR[testIndex,]=rfOutMax$test$votes
        } else {
          yPredMaxR[testIndex]=rfOutMax$test$predicted
        }
        if (modReturn) {
          outMod[[i]]=list(rfOutMin,rfOutMid,rfOutMax)
        }
      }
    }

    # PER REPETITION:

    # Average outer loop predictions, VIP ranks and nComp for PLS
    parReturn=list(yPredMin=yPredMinR,yPredMid=yPredMidR,yPredMax=yPredMaxR)
    parReturn$VIPRepMin=apply(VIPOutMin,1,mean)
    parReturn$VIPRepMid=apply(VIPOutMid,1,mean)
    parReturn$VIPRepMax=apply(VIPOutMax,1,mean)
    if (method=='PLS'){
      parReturn$nCompRepMin=round(mean(nCompOutMin))
      parReturn$nCompRepMid=round(mean(nCompOutMid))
      parReturn$nCompRepMax=round(mean(nCompOutMax))
      parReturn$nCompSegMin=nCompOutMin
      parReturn$nCompSegMid=nCompOutMid
      parReturn$nCompSegMax=nCompOutMax
    }
    # Validation curves
    parReturn$VAL=VALRep
    # Calculate nVar per repetition
    fitRankRep=colSums(VALRep)
    if(fitness=='AUROC') fitRankRep=-fitRankRep
    fitRankRep=(fitRankRep-min(fitRankRep))/abs(diff(range(fitRankRep)))
    if(all(is.nan(fitRankRep))) fitRankRep=rep(0,cnt) # If all VAL have same value -> reset all fitRankRep to 0
    minIndex=max(which(fitRankRep<=methParam$robust))
    maxIndex=min(which(fitRankRep<=methParam$robust))
    parReturn$varRepMin=var[minIndex]
    parReturn$varRepMid=round(exp(mean(log(c(var[minIndex],var[maxIndex])))))
    parReturn$varRepMax=var[maxIndex]
    # Return underlying models
    if (modReturn) parReturn$outModel=outMod
    if (logg) sink()
    return(parReturn)
    # reps[[r]]=parReturn
  }


  # Unpack the different repetitions
  if (modReturn) outMods=list()
  for (r in 1:nRep) {
    if (DA) yPredMin[,,r]=reps[[r]]$yPredMin else
      yPredMin[,r]=reps[[r]]$yPredMin
    if (DA) yPredMid[,,r]=reps[[r]]$yPredMid else
      yPredMid[,r]=reps[[r]]$yPredMid
    if (DA) yPredMax[,,r]=reps[[r]]$yPredMax else
      yPredMax[,r]=reps[[r]]$yPredMax
    varRepMin[r]=reps[[r]]$varRepMin
    varRepMid[r]=reps[[r]]$varRepMid
    varRepMax[r]=reps[[r]]$varRepMax
    VIPRepMin[,r]=reps[[r]]$VIPRepMin
    VIPRepMid[,r]=reps[[r]]$VIPRepMid
    VIPRepMax[,r]=reps[[r]]$VIPRepMax
    if (method=='PLS') {
      nCompRepMin[r]=reps[[r]]$nCompRepMin
      nCompRepMid[r]=reps[[r]]$nCompRepMid
      nCompRepMax[r]=reps[[r]]$nCompRepMax
      nCompSegMin[r,]=reps[[r]]$nCompSegMin
      nCompSegMid[r,]=reps[[r]]$nCompSegMid
      nCompSegMax[r,]=reps[[r]]$nCompSegMax
    }
    VAL[,,r]=reps[[r]]$VAL
    if (modReturn) outMods=c(outMods,reps[[r]]$outModel)
  }
  # Average predictions
  if (DA) {
    yPred=list()
    yPred[['min']]=apply(yPredMin,c(1,2),mean)
    yPred[['mid']]=apply(yPredMid,c(1,2),mean)
    yPred[['max']]=apply(yPredMax,c(1,2),mean)
  } else {
    yPred=apply(yPredMin,1,mean)[1:nrow(X)]
    yPred=cbind(yPred,apply(yPredMid,1,mean)[1:nrow(X)])
    yPred=cbind(yPred,apply(yPredMax,1,mean)[1:nrow(X)])
    colnames(yPred)=c('min','mid','max')
    rownames(yPred)=paste(1:nSamp,ID,sep='_ID')
  }
  modelReturn$yPred=yPred
  if (DA) {
    auc=matrix(nrow=3,ncol=length(levels(Y)),dimnames=list(c('min','mid','max'),levels(Y)))
    for (cl in 1:length(levels(Y))) {
      auc[1,cl]=roc(Y==(levels(Y)[cl]),yPred[['min']][,cl], quiet = TRUE)$auc
      auc[2,cl]=roc(Y==(levels(Y)[cl]),yPred[['mid']][,cl], quiet = TRUE)$auc
      auc[3,cl]=roc(Y==(levels(Y)[cl]),yPred[['max']][,cl], quiet = TRUE)$auc
    }
    # Classify predictions
    miss=numeric(3)
    yClass=data.frame(Y)
    for (mo in 1:3) {
      classPred=factor(apply(yPred[[mo]],1,function(x) levels(Y)[which.max(x)]),levels=levels(Y))
      miss[mo]=sum(classPred!=Y)
      yClass[,mo]=classPred
    }
    names(miss)=colnames(yClass)=c('min','mid','max')
    rownames(yClass)=paste(1:nSamp,ID,sep='_ID')
    # Report
    modelReturn$yClass=yClass
    modelReturn$miss=miss
    modelReturn$auc=auc
  } else if (ML) {
    modelReturn$yClass=apply(yPred,2,function(x) ifelse(x>0,1,-1))
    modelReturn$miss=apply(modelReturn$yClass,2,function(x) sum(x!=Y))
    modelReturn$auc=apply(yPred,2,function(x) roc(Y,x)$auc)
    colnames(modelReturn$yClass)=names(modelReturn$miss)=names(modelReturn$auc)=c('min','mid','max')
    rownames(modelReturn$yClass)=paste(1:nSamp,ID,sep='_ID')
  }
  # Average VIP ranks over repetitions
  VIP=apply(VIPRepMin,1,mean)
  VIP=cbind(VIP,apply(VIPRepMid,1,mean))
  VIP=cbind(VIP,apply(VIPRepMax,1,mean))
  colnames(VIP)=c('min','mid','max')
  modelReturn$VIP=VIP
  modelReturn$VIPPerRep=list(minModel=VIPRepMin,midModel=VIPRepMid,maxModel=VIPRepMax)
  # Calculate overall nVar
  fitRankAll=apply(VAL,2,mean)
  if(fitness=='AUROC') fitRankAll=-fitRankAll
  fitRankAll=(fitRankAll-min(fitRankAll))/abs(diff(range(fitRankAll))) # rescale to 0-1 range
  if(all(is.nan(fitRankAll))) fitRankAll=rep(0,cnt) # If all VAL have same value -> reset all fitRankAll to 0
  minIndex=max(which(fitRankAll<=methParam$robust))
  maxIndex=min(which(fitRankAll<=methParam$robust))
  nVar=c(var[minIndex],round(exp(mean(log(c(var[minIndex],var[maxIndex]))))),var[maxIndex])
  names(nVar)=c('min','mid','max')
  modelReturn$nVar=nVar
  if (method=='PLS') {
    # Average nComp over repetitions
    nComp=c(round(mean(nCompRepMin)),round(mean(nCompRepMid)),round(mean(nCompRepMax)))
    names(nComp)=c('min','mid','max')
    modelReturn$nComp=nComp
  }
  modelReturn$VAL$metric=fitness
  modelReturn$VAL$VAL=VAL
  if (modReturn) modelReturn$outModels=outMods
  modelReturn$yPredPerRep=list(minModel=yPredMin,midModel=yPredMid,maxModel=yPredMax)
  modelReturn$nVarPerRep=list(minModel=varRepMin,midModel=varRepMid,maxModel=varRepMax)
  if (method=='PLS') {
    modelReturn$nCompPerRep=list(minModel=nCompRepMin,midModel=nCompRepMid,maxModel=nCompRepMax)
    modelReturn$nCompPerSeg=list(minModel=nCompSegMin,midModel=nCompSegMid,maxModel=nCompSegMax)
  }
  modelReturn$inData=InData
  ## Build overall "Fit" method for calculating R2 and visualisations
  incVarMin=names(VIP[rank(VIP[,1])<=round(nVar[1]),1])
  incVarMid=names(VIP[rank(VIP[,2])<=round(nVar[2]),2])
  incVarMax=names(VIP[rank(VIP[,3])<=round(nVar[3]),3])
  if (method=='PLS'){
    # Min model
    if (DA) plsFitMin=MUVR::plsda(subset(X,select=incVarMin),Y,ncomp=round(nComp[1]),near.zero.var=TRUE,scale=scale) else
      plsFitMin=MUVR::pls(subset(X,select=incVarMin),Y,ncomp=round(nComp[1]),near.zero.var=TRUE,scale=scale)
    if (length(plsFitMin$nzv$Position)>0) incVarMin=incVarMin[!incVarMin%in%rownames(plsFitMin$nzv$Metrics)]
    yFitMin=predict(plsFitMin,newdata=subset(X,select=incVarMin),scale=scale)$predict[,,nComp[1]]  #
    # Mid model
    if (DA) plsFitMid=MUVR::plsda(subset(X,select=incVarMid),Y,ncomp=round(nComp[2]),near.zero.var=TRUE,scale=scale) else
      plsFitMid=MUVR::pls(subset(X,select=incVarMid),Y,ncomp=round(nComp[2]),near.zero.var=TRUE,scale=scale)
    if (length(plsFitMid$nzv$Position)>0) incVarMid=incVarMid[!incVarMid%in%rownames(plsFitMid$nzv$Metrics)]
    yFitMid=predict(plsFitMid,newdata=subset(X,select=incVarMid),scale=scale)$predict[,,nComp[2]]  #
    # Max model
    if (DA) plsFitMax=MUVR::plsda(subset(X,select=incVarMax),Y,ncomp=round(nComp[3]),near.zero.var=TRUE,scale=scale) else
      plsFitMax=MUVR::pls(subset(X,select=incVarMax),Y,ncomp=round(nComp[3]),near.zero.var=TRUE,scale=scale)
    if (length(plsFitMax$nzv$Position)>0) incVarMax=incVarMax[!incVarMax%in%rownames(plsFitMax$nzv$Metrics)]
    yFitMax=predict(plsFitMax,newdata=subset(X,select=incVarMax),scale=scale)$predict[,,nComp[3]]  #
    yFit=cbind(yFitMin,yFitMid,yFitMax)
    yRep=ncol(yFit)/3
    colnames(yFit)=rep(c('min','mid','max'),each=yRep)
    rownames(yFit)=ID
    modelReturn$Fit=list(yFit=yFit,plsFitMin=plsFitMin,plsFitMid=plsFitMid,plsFitMax=plsFitMax)
  } else {
    rfFitMin=suppressWarnings(randomForest(subset(X,select=incVarMin),Y)) # suppress warnings for ML (regression against fewer than 5 unique values)
    if (DA) {
      yFitMin=rfFitMin$votes
    } else {
      yFitMin=rfFitMin$predicted
    }
    rfFitMid=suppressWarnings(randomForest(subset(X,select=incVarMid),Y)) # suppress warnings for ML (regression against fewer than 5 unique values)
    if (DA) {
      yFitMid=rfFitMid$votes
    } else {
      yFitMid=rfFitMid$predicted
    }
    rfFitMax=suppressWarnings(randomForest(subset(X,select=incVarMax),Y)) # suppress warnings for ML (regression against fewer than 5 unique values)
    if (DA) {
      yFitMax=rfFitMax$votes
    } else {
      yFitMax=rfFitMax$predicted
    }
    yFit=cbind(yFitMin,yFitMid,yFitMax)
    yRep=ncol(yFit)/3
    colnames(yFit)=rep(c('min','mid','max'),each=yRep)
    rownames(yFit)=ID
    modelReturn$Fit=list(yFit=yFit,rfFitMin=rfFitMin,rfFitMid=rfFitMid,rfFitMax=rfFitMax)
  }
  # Calculate fit statistics
  if (!DA) {
    TSS=sum((Y-mean(Y))^2)
    RSSMin=sum((Y-yFitMin)^2)
    RSSMid=sum((Y-yFitMid)^2)
    RSSMax=sum((Y-yFitMax)^2)
    PRESSMin=sum((Y-yPred[,1])^2)
    PRESSMid=sum((Y-yPred[,2])^2)
    PRESSMax=sum((Y-yPred[,3])^2)
    R2Min=1-(RSSMin/TSS)
    R2Mid=1-(RSSMid/TSS)
    R2Max=1-(RSSMax/TSS)
    Q2Min=1-(PRESSMin/TSS)
    Q2Mid=1-(PRESSMid/TSS)
    Q2Max=1-(PRESSMax/TSS)
    modelReturn$fitMetric=list(R2=c(R2Min,R2Mid,R2Max),Q2=c(Q2Min,Q2Mid,Q2Max))
  } else {
    modelReturn$fitMetric <- list(CR = 1 - (miss / length(Y)))
  }
  # Stop timer
  end.time=proc.time()[3]
  modelReturn$calcMins=(end.time-start.time)/60
  cat('\n Elapsed time',(end.time-start.time)/60,'mins \n')
  class(modelReturn)=c('MVObject',method,ifelse(DA,'Classification',ifelse(ML,'Multilevel','Regression')))
  return(modelReturn)
}
