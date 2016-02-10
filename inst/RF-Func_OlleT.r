### Function for random sampling with "correct" proportion between groups
vectSamp=function(vect,n=4,sampLen) { # Pick 'n' random samples within vector 'vect'
	# sampLen is a vector of number of observations within sample
	# If sampLen is not specified it is automatically calculated
	j=length(vect)
	if (missing(sampLen)) { # Calculation of sampLen vector
		sampLen=numeric(n)
		oddAdd=numeric(n)
		if (j%%n!=0) oddAdd[1:(j%%n)]=1
		sampLen=floor(j/n)+oddAdd
	}
	if (length(sampLen)!=n) 
		stop("Mismatch number of samples")
	if (sum(sampLen)>j) 
		stop("Exceeds maximum samples")
	if (sum(sampLen)<j) 
		cat('\n Warning: Undersampling \n')
	vectSamp=sample(vect) # Randomise vector vect
	samp=list()
	vectInd=1 # Index in randomised vector
	for (i in 1:n) { # Divide randomised vector into samples
		samp[[i]]=sort(vectSamp[vectInd:(vectInd+sampLen[i]-1)])  # Make sample
		vectInd=vectInd+sampLen[i]  # Count up index
	}
	return(samp)
}


### parallelized rdCV
otusRF=function(X,Y,nRep=12,nSeg=8,nInner=7,featRatio=0.75) {  # Perform rdCV PLS 
	start.time=proc.time()[3]
	modelReturn=list()
	nFeat=nFeat0=ncol(X)
	YIndex=1:length(Y)
	name=unique(Y)  # Find groups
	groups=length(name) # Number of groups
	groupInd=list()  # Allocate list for indices of groups
	for (g in 1:groups) { 
	  groupInd[[g]]=YIndex[Y==name[g]]  # Find indices per group
	}
	VIPMin=VIPMid=VIPMax=matrix(nrow=nFeat0,ncol=3) # Allocate matrix for variable importance averaged over repetition
	rownames(VIPMin)=rownames(VIPMid)=rownames(VIPMax)=colnames(X)
	colnames(VIPMin)=colnames(VIPMid)=colnames(VIPMax)=c('Mean','SD','CV')
	# Allocate matrix for prediction of outer segments Y per repetition
	yPredMin=yPredMid=yPredMax=array(dim=c(length(Y),length(unique(Y)),nRep))
	yPredMinR=yPredMidR=yPredMaxR=matrix(nrow=length(Y),ncol=length(unique(Y)))
	# Allocate response vectors and matrices for feat's, nComp and VIP ranks over repetitions
	featRepMin=featRepMid=featRepMax=numeric(nRep)
	names(featRepMin)=names(featRepMid)=names(featRepMax)=paste(rep('rep',nRep),1:nRep,sep='')
	VIPRepMin=VIPRepMid=VIPRepMax=matrix(data=nFeat0,nrow=nFeat0,ncol=nRep)
	rownames(VIPRepMin)=rownames(VIPRepMid)=rownames(VIPRepMax)=colnames(X)
	colnames(VIPRepMin)=colnames(VIPRepMid)=colnames(VIPRepMax)=paste(rep('rep',nRep),1:nRep,sep='')
	feat=numeric()
	cnt=0
	while (nFeat>2) {  
		cnt=cnt+1
		feat=c(feat,nFeat)
		nFeat=floor(featRatio*nFeat)
	}
	# Perform rdCV 
	reps=foreach(r=1:nRep, .packages=c('randomForest','pROC'), .export='vectSamp') %dopar% {
	  sink('log.txt',append=TRUE)
		cat('\n','   Repetition ',r,' of ',nRep,':',sep='')
	  groupTest=list()  ## Allocate list for samples within group
	  for (gT in 1:groups) { 
	    groupTest[[gT]]=vectSamp(groupInd[[gT]],n=nSeg)  # Draw random samples within group
	  }
	  allTest=groupTest[[1]] # Add 1st groups to 'Master' sample of all groups
	  for (gT in 2:groups) {  # Add subsequent groups
	    allTest=allTest[order(sapply(allTest,length))]
	    for (aT in 1:nSeg) {
	      allTest[[aT]]=sort(c(allTest[[aT]],groupTest[[gT]][[aT]]))
	    }
	  }
	  featOutMin=featOutMid=featOutMax=numeric(nSeg)
		names(featOutMin)=names(featOutMid)=names(featOutMax)=paste(rep('outSeg',nSeg),1:nSeg,sep='')
		VIPOutMin=VIPOutMid=VIPOutMax=matrix(data=nFeat0,nrow=nFeat0,ncol=nSeg)
		rownames(VIPOutMin)=rownames(VIPOutMid)=rownames(VIPOutMax)=colnames(X)
		colnames(VIPOutMin)=colnames(VIPOutMid)=colnames(VIPOutMax)=paste(rep('outSeg',nSeg),1:nSeg,sep='')
		for (i in 1:nSeg) {   # Create 'nSamp' models within repetition
			# i=1
			# cat('\n Segment ',i,' (features):',sep='') # Counter
			testIndex=allTest[[i]] # Draw out segment = holdout set BASED ON UNIK
			xTest=X[testIndex,]
			yTest=Y[testIndex]
			InnerIndex=YIndex[-testIndex]
			InnerID=Y[-testIndex]
			missIn=matrix(nrow=nInner,ncol=cnt)
			rownames(missIn)=paste(rep('inSeg',nInner),1:nInner,sep='')
			colnames(missIn)=feat
			VIPInner=array(data=nFeat0,dim=c(nFeat0,cnt,nInner))
			rownames(VIPInner)=colnames(X)
			colnames(VIPInner)=feat
			dimnames(VIPInner)[[3]]=paste(rep('inSeg',nInner),1:nInner,sep='')
			nFeat=nFeat0
			incFeat=colnames(X)
			for (count in 1:cnt) {  # Build models with successively fewer feature. Quality metric = number of missclassifications for Validation set
				nFeat=feat[count]
				cat('\n Rep ',r,'; Seg ',i,'; ',nFeat,' features.')
				groupIndVal=list()
				for (g in 1:groups) { 
				  groupIndVal[[g]]=InnerIndex[InnerID==name[g]]  # Find indices per group
				}
				groupVal=list()  ## Allocate list for samples within group
				for (gV in 1:groups) { 
				  groupVal[[gV]]=vectSamp(groupIndVal[[gV]],n=nInner)  # Draw random samples within group
				}
				allVal=groupVal[[1]] # Add 1st groups to 'Master' sample of all groups
				for (gV in 2:groups) {  # Add subsequent groups
				  allVal=allVal[order(sapply(allVal,length))]
				  for (aV in 1:nInner) {
				    allVal[[aV]]=sort(c(allVal[[aV]],groupVal[[gV]][[aV]]))
				  }
				}
				for (j in 1:nInner) {
					valIndex=allVal[[j]] # Draw out segment = validation set
					xVal=X[valIndex,]
					xVal=subset(xVal,select=incFeat)
					yVal=Y[valIndex]
					trainInd=InnerIndex[-match(valIndex,InnerIndex)] # Define Training segment
					xTrain=X[trainInd,]
					xTrain=subset(xTrain,select=incFeat)
					yTrain=Y[trainInd]
					# perform RF training
					rfInner=randomForest(xTrain,yTrain,xVal,yVal,ntree=300)
					missIn[j,count]=sum(rfInner$test$predicted!=yVal)
					imp=rfInner$importance
					VIPInner[match(rownames(imp),rownames(VIPInner)),count,j]=rank(-imp)
				}
				VIPInAve=apply(VIPInner[,count,],1,mean)
				if (count<cnt) {
					incFeat=names(VIPInAve[order(VIPInAve)])[1:feat[count+1]]
				}
			}
			minIndex=max(which(apply(t(apply(missIn,1,rank)),2,mean)==min(apply(t(apply(missIn,1,rank)),2,mean))))
			maxIndex=min(which(apply(t(apply(missIn,1,rank)),2,mean)==min(apply(t(apply(missIn,1,rank)),2,mean))))
			midIndex=which.min(abs(feat-exp(mean(log(c(feat[minIndex],feat[maxIndex]))))))  # Geometric
			# Per outer segment: Average inner loop features, nComp and VIP ranks 
			featOutMin[i]=feat[minIndex]
			featOutMid[i]=feat[midIndex]
			featOutMax[i]=feat[maxIndex]
			VIPOutMin[,i]=apply(VIPInner[,minIndex,],1,mean)
			VIPOutMid[,i]=apply(VIPInner[,midIndex,],1,mean)
			VIPOutMax[,i]=apply(VIPInner[,maxIndex,],1,mean)
			# Build outer model for min and max nComp and predict YTEST
			xIn=X[-testIndex,] # Perform Validation on all samples except holdout set
			yIn=Y[-testIndex]
			incFeatMin=rownames(VIPOutMin)[rank(VIPOutMin[,i])<=featOutMin[i]]
			xTrainMin=subset(xIn,select=incFeatMin)
			xTestMin=subset(xTest,select=incFeatMin)
			rfOutMin=randomForest(xTrainMin,yIn,xTestMin,yTest,ntree=500)
			yPredMinR[testIndex,]=rfOutMin$test$votes
			incFeatMid=rownames(VIPOutMid)[rank(VIPOutMid[,i])<=featOutMid[i]]
			xTrainMid=subset(xIn,select=incFeatMid)
			xTestMid=subset(xTest,select=incFeatMid)
			rfOutMid=randomForest(xTrainMid,yIn,xTestMid,yTest,ntree=500)
			yPredMidR[testIndex,]=rfOutMid$test$votes
			incFeatMax=rownames(VIPOutMax)[rank(VIPOutMax[,i])<=featOutMax[i]]
			xTrainMax=subset(xIn,select=incFeatMax)
			xTestMax=subset(xTest,select=incFeatMax)
			rfOutMax=randomForest(xTrainMax,yIn,xTestMax,yTest,ntree=500)
			yPredMaxR[testIndex,]=rfOutMax$test$votes
		}
		# Per repetition: Average outer loop features, nComp and VIP ranks 
		featRepMinR=round(mean(featOutMin))
		VIPRepMinR=apply(VIPOutMin,1,mean)
		featRepMidR=round(mean(featOutMid))
		VIPRepMidR=apply(VIPOutMid,1,mean)
		featRepMaxR=round(mean(featOutMax))
		VIPRepMaxR=apply(VIPOutMax,1,mean)
		parReturn=list(yPredMin=yPredMinR,featRepMin=featRepMinR,VIPRepMin=VIPRepMinR,
			yPredMid=yPredMidR,featRepMid=featRepMidR,VIPRepMid=VIPRepMidR,
			yPredMax=yPredMaxR,featRepMax=featRepMaxR,VIPRepMax=VIPRepMaxR)
		return(parReturn)
	}
	for (r in 1:nRep) {
		yPredMin[,,r]=reps[[r]]$yPredMin
		featRepMin[r]=reps[[r]]$featRepMin
		VIPRepMin[,r]=reps[[r]]$VIPRepMin
		yPredMid[,,r]=reps[[r]]$yPredMid
		featRepMid[r]=reps[[r]]$featRepMid
		VIPRepMid[,r]=reps[[r]]$VIPRepMid
		yPredMax[,,r]=reps[[r]]$yPredMax
		featRepMax[r]=reps[[r]]$featRepMax
		VIPRepMax[,r]=reps[[r]]$VIPRepMax
	}
	
	yMin=yMid=yMax=array(dim=c(length(Y),length(unique(Y)),3))
	colnames(yMin)=colnames(yMid)=colnames(yMax)=name
	dimnames(yMin)[[3]]=dimnames(yMid)[[3]]=dimnames(yMax)[[3]]=c('Mean','SD','CV')
	yMin[,,1]=apply(yPredMin,1:2,mean)
	yMin[,,2]=apply(yPredMin,1:2,sd)
	yMin[,,3]=abs(100*yMin[,,2]/yMin[,,1])
	yMid[,,1]=apply(yPredMid,1:2,mean)
	yMid[,,2]=apply(yPredMid,1:2,sd)
	yMid[,,3]=abs(100*yMid[,,2]/yMid[,,1])
	yMax[,,1]=apply(yPredMax,1:2,mean)
	yMax[,,2]=apply(yPredMax,1:2,sd)
	yMax[,,3]=abs(100*yMax[,,2]/yMax[,,1])
	# Classify predictions
	YMinC=t(apply(yMin[,,1],1,function(x) ifelse(x==max(x),1,0)))
	YMidC=t(apply(yMid[,,1],1,function(x) ifelse(x==max(x),1,0)))
	YMaxC=t(apply(yMax[,,1],1,function(x) ifelse(x==max(x),1,0)))
	y1=ifelse(Y==unique(Y)[1],1,0)
	y2=ifelse(Y==unique(Y)[2],1,0)
	y3=ifelse(Y==unique(Y)[3],1,0)
	YC=cbind(y1,y2,y3)
	colnames(YC)=unique(Y)
	# Calculate misclassifications
	missMin=nrow(YC)-sum(YC*YMinC)
	missMid=nrow(YC)-sum(YC*YMidC)
	missMax=nrow(YC)-sum(YC*YMaxC)
	# Average VIP ranks over repetitions
	VIPMin[,1]=apply(VIPRepMin,1,mean)
	VIPMin[,2]=apply(VIPRepMin,1,sd)
	VIPMin[,3]=100*abs(VIPMin[,2]/VIPMin[,1])
	VIPMid[,1]=apply(VIPRepMid,1,mean)
	VIPMid[,2]=apply(VIPRepMid,1,sd)
	VIPMid[,3]=100*abs(VIPMid[,2]/VIPMid[,1])
	VIPMax[,1]=apply(VIPRepMax,1,mean)
	VIPMax[,2]=apply(VIPRepMax,1,sd)
	VIPMax[,3]=100*abs(VIPMax[,2]/VIPMax[,1])
	# Average nFeat over repetitions
	nFeatMin=c(mean(featRepMin),sd(featRepMin),abs(100*sd(featRepMin)/mean(featRepMin)))
	nFeatMid=c(mean(featRepMid),sd(featRepMid),abs(100*sd(featRepMid)/mean(featRepMid)))
	nFeatMax=c(mean(featRepMax),sd(featRepMax),abs(100*sd(featRepMax)/mean(featRepMax)))
	names(nFeatMin)=names(nFeatMid)=names(nFeatMax)=c('Mean','SD','CV')
	# Stop timer
	end.time=proc.time()[3]
	# Generate report-list
	modelReturn$nFeat=list(minModel=nFeatMin,midModel=nFeatMid,maxModel=nFeatMax)
	modelReturn$nFeatPerRep=list(minModel=featRepMin,midModel=featRepMid,maxModel=featRepMax)
	modelReturn$yPred=list(minModel=yMin,midModel=yMid,maxModel=yMax)
	modelReturn$yPredPerRep=list(minModel=yPredMin,midModel=yPredMid,maxModel=yPredMax)
	modelReturn$VIPRank=list(minModel=VIPMin,midModel=VIPMid,maxModel=VIPMax)
	modelReturn$yClass=list(minModel=YMinC,midModel=YMidC,maxModel=YMaxC)
	modelReturn$misClass=list(minModel=missMin,midModel=missMid,maxModel=missMax)
	modelReturn$calcMins=(end.time-start.time)/60
	cat('\n Elapsed time',(end.time-start.time)/60,'mins \n')
	return(modelReturn)
}
