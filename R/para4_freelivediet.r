rm(list=ls())
# getwd()
setwd('C:/R/freelivediet')
library(mixOmics)
library(pROC)
library(doParallel)
library(foreach)

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

### non-parallelized rdCV    ## Extra special: some individuals are present more than once - sampling is performed accordingly. That is why there are many different indices.
ryeCV4=function(X,Y,ID,nRep=10,nSeg=7,nInner=6,comps=5,metric=c('miss','auc','rmsep'),DA=FALSE,featRatio=0.75,P=FALSE,XP=NA) {  # Perform rdCV PLS 
	start.time=proc.time()[3]
	if (missing(metric)) metric='rmsep'
	if (P==T) {  # If model is also used to predict XP matrix, allocate a YP response array
		YP=matrix(nrow=nrow(XP),ncol=nSeg*nRep)  # Predictions for all samples (rows) across outer segments (cols) and repetitions (dim3)
	}
	modelReturn=list()
	nFeat=nFeat0=ncol(X)
	comp=ifelse(nFeat<comps,nFeat,comps)
	unik=!duplicated(ID)
	IDunik=ID[unik]
	Yunik=Y[unik]
	YIndex=1:length(Y)
	ind=1:length(Yunik) # Make indices of unique individuals INSTEAD OF OBSERVATIONS!!!
	if (DA==T) {
		name=unique(Y)  # Find groups
		groups=length(name) # Number of groups
		groupInd=list()  # Allocate list for indices of groups
		for (g in 1:groups) { 
			groupInd[[g]]=ind[Yunik==name[g]]  # Find indices per group
		}
	}
	VIPMin=VIPMid=VIPMax=matrix(nrow=nFeat0,ncol=3) # Allocate matrix for variable importance averaged over repetition
	rownames(VIPMin)=rownames(VIPMid)=rownames(VIPMax)=colnames(X)
	colnames(VIPMin)=colnames(VIPMid)=colnames(VIPMax)=c('Mean','SD','CV')
	# Allocate matrix for prediction of outer segments Y per repetition
	yPredMin=yPredMid=yPredMax=matrix(nrow=length(Y),ncol=nRep)
	colnames(yPredMin)=colnames(yPredMid)=colnames(yPredMax)=paste('Rep',1:nRep,sep='')
	yPredMinR=yPredMidR=yPredMaxR=numeric(length(Y))
	# Allocate response vectors and matrices for feat's, nComp and VIP ranks over repetitions
	featRepMin=featRepMid=featRepMax=nCompRepMin=nCompRepMid=nCompRepMax=missRep=numeric(nRep)
	names(featRepMin)=names(featRepMid)=names(featRepMax)=names(nCompRepMin)=names(nCompRepMid)=names(nCompRepMax)=names(missRep)=paste(rep('rep',nRep),1:nRep,sep='')
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
	# 
	reps=list()
	for (r in 1:nRep) {  # Loop for all repetitions
	# reps=foreach(r=1:nRep, .packages=c('mixOmics','pROC'), .export='vectSamp') %dopar% {
		# r=1
		if (P==T) YPR=matrix(nrow=nrow(XP),ncol=nSeg)
		cat('\n','   Repetition ',r,' of ',nRep,':',sep='')
		if (DA==T) {
			groupTest=list()  ## Allocate list for samples within group
			for (gT in 1:groups) { 
				groupTest[[gT]]=vectSamp(groupInd[[gT]],n=nSeg)  # Draw random samples within group
			}
			groupTest[[2]]=rev(groupTest[[2]])  ## Reverse order for 2nd group 
			## In the case of >2 groups another approach would be suitable
			## eg adding subsequent groups to allTest ordered by length
			allTest=groupTest[[1]] # Add 1st groups to 'Master' sample of all groups
			for (gT in 2:groups) {  # Add subsequent groups
				for (aT in 1:nSeg) {
					allTest[[aT]]=sort(c(allTest[[aT]],groupTest[[gT]][[aT]]))
				}
			}
		} else {
			allTest=vectSamp(ind,n=nSeg)
		}
		featOutMin=featOutMid=featOutMax=nCompOutMin=nCompOutMid=nCompOutMax=numeric(nSeg)
		names(featOutMin)=names(featOutMid)=names(featOutMax)=names(nCompOutMin)=names(nCompOutMax)=paste(rep('outSeg',nSeg),1:nSeg,sep='')
		VIPOutMin=VIPOutMid=VIPOutMax=matrix(data=nFeat0,nrow=nFeat0,ncol=nSeg)
		rownames(VIPOutMin)=rownames(VIPOutMid)=rownames(VIPOutMax)=colnames(X)
		colnames(VIPOutMin)=colnames(VIPOutMid)=colnames(VIPOutMax)=paste(rep('outSeg',nSeg),1:nSeg,sep='')
		for (i in 1:nSeg) {   # Create 'nSamp' models within repetition
			# i=1
			cat('\n Segment ',i,' (features):',sep='') # Counter
			testInd=allTest[[i]] # Draw out segment = holdout set BASED ON UNIK
			TestIndex=numeric()  # BASED ON OBSERVATIONS
			for (xt in 1:length(testInd)) {
				# cat(sum(ID==IDunik[testInd[xt]]))
				TestIndex=c(TestIndex,which(ID==IDunik[testInd[xt]]))
			}
			xTest=X[TestIndex,]
			yTest=Y[TestIndex]
			inInd=ind[-testInd]
			InnerIndex=YIndex[-TestIndex]
			yUnikIn=Yunik[inInd]  
			###
			
			# cat('\n',head(testInd))   # Trouble shooting
			# cat('\n',head(yTrain))
			# cat('\n',as.numeric(xTrain[1,1:5]))
			missIn=aucIn=rmsepIn=nCompIn=matrix(nrow=nInner,ncol=cnt)
			rownames(rmsepIn)=rownames(missIn)=rownames(aucIn)=rownames(nCompIn)=paste(rep('inSeg',nInner),1:nInner,sep='')
			colnames(rmsepIn)=colnames(missIn)=colnames(aucIn)=colnames(nCompIn)=feat
			VIPInner=array(data=nFeat0,dim=c(nFeat0,cnt,nInner))
			rownames(VIPInner)=colnames(X)
			colnames(VIPInner)=feat
			dimnames(VIPInner)[[3]]=paste(rep('inSeg',nInner),1:nInner,sep='')
			nFeat=nFeat0
			# inSamp=vectSamp(inInd,n=nInner)  # Draw random samples within current 2CV segment
			incFeat=colnames(X)
			for (count in 1:cnt) {  # Build models with successively fewer feature. Quality metric = number of missclassifications for Validation set
				# count=1 # for testing
				# count=count+1
				nFeat=feat[count]
				cat(nFeat)
				# if (resampInner==TRUE) inSamp=vectSamp(inInd,n=nInner)  # Resample inner segments for each feature elimination step
				comp=ifelse(nFeat<comps,nFeat,comps)
				if (DA==T) {
					groupIndVal=list()
					for (g in 1:groups) { 
						groupIndVal[[g]]=inInd[yUnikIn==name[g]]  # Find indices per group
					}
					groupVal=list()  ## Allocate list for samples within group
					for (gV in 1:groups) { 
						groupVal[[gV]]=vectSamp(groupIndVal[[gV]],n=nInner)  # Draw random samples within group
					}
					groupVal[[2]]=rev(groupVal[[2]])  ## Reverse order for 2nd group 
					## In the case of >2 groups another approach would be suitable
					## eg adding subsequent groups to allVal ordered by length
					allVal=groupVal[[1]] # Add 1st groups to 'Master' sample of all groups
					for (gV in 2:groups) {  # Add subsequent groups
						for (aV in 1:nInner) {
							allVal[[aV]]=sort(c(allVal[[aV]],groupVal[[gV]][[aV]]))
						}
					}
				} else {
					allVal=vectSamp(inInd,n=nInner)
				}
				for (j in 1:nInner) {
					# j=1 # for testing
					# j=j+1
					cat('.') # Counter
					valInd=allVal[[j]] # Draw out segment = validation set
					ValIndex=numeric()  # BASED ON OBSERVATIONS
					for (xv in 1:length(valInd)) {
						# cat(sum(ID==IDunik[valInd[xv]]))
						ValIndex=c(ValIndex,which(ID==IDunik[valInd[xv]]))
					}
					xVal=X[ValIndex,]
					xVal=subset(xVal,select=incFeat)
					yVal=Y[ValIndex]
					trainInd=InnerIndex[-match(ValIndex,InnerIndex)] # Define Training segment
					xTrain=X[trainInd,]
					xTrain=subset(xTrain,select=incFeat)
					yTrain=Y[trainInd]
					# perform PLS up to 'comps' number of components
					plsInner=pls(xTrain,yTrain,ncomp=comp,mode="regression")
					if (length(plsInner$nzv$Position)>0) {
						removeFeat=rownames(plsInner$nzv$Metrics)
						xVal=xVal[,!colnames(xVal)%in%removeFeat]
					}
					yValInner=predict(plsInner,newdata=xVal)$predict[,,]  # Store  prediction estimates per validation segment 
					if (metric=='miss') {
						# cat(' miss',count)
						yClassInner=ifelse(yValInner>0,1,-1)
						misClass=apply((yVal-yClassInner)*yVal/2,2,sum,na.rm=T)
						missIn[j,count]=min(misClass)
						nCompIn[j,count]=which.min(misClass)
					} 
					if (metric=='auc') {
						# cat(' auc',count)
						auc=apply(yValInner,2,function(x) roc(yVal,x)$auc)
						aucIn[j,count]=max(auc)
						nCompIn[j,count]=which.max(auc)
					}
					if (metric=='rmsep') {
						# cat(' rmsep',count)
						rmsep=apply(yValInner,2,function(x) sqrt(sum((yVal-x)^2,na.rm=T)/(length(yValInner[,1])-sum(is.na(yValInner[,1])))))
						rmsepIn[j,count]=min(rmsep)
						nCompIn[j,count]=which.min(rmsep)
					}
					VIPInner[match(names(vip(plsInner)[,nCompIn[j,count]]),rownames(VIPInner)),count,j]=rank(-vip(plsInner)[,nCompIn[j,count]])
					# VIPInner[match(names(vip(plsInner)[,nCompIn[j]]),rownames(VIPInner)),j]=rank(-vip(plsInner)[,nCompIn[j]])				
				}
				## NB!!! Average VIP ranks over inner segments before feature elimination!!!
				VIPInAve=apply(VIPInner[,count,],1,mean)
				if (count<cnt) {
					incFeat=names(VIPInAve[order(VIPInAve)])[1:feat[count+1]]
				}
			}
			if (metric=='auc') {
				minIndex=max(which(apply(t(apply(aucIn,1,rank)),2,mean)==max(apply(t(apply(aucIn,1,rank)),2,mean))))
				maxIndex=min(which(apply(t(apply(aucIn,1,rank)),2,mean)==max(apply(t(apply(aucIn,1,rank)),2,mean))))
				# Find a middle index | Either arithmetic or geometric mean
				# midIndex=which.min(abs(feat-mean(c(feat[minIndex],feat[maxIndex]))))  # Arithmetic
				midIndex=which.min(abs(feat-exp(mean(log(c(feat[minIndex],feat[maxIndex]))))))  # Geometric
			}
			if (metric=='miss') {
				minIndex=max(which(apply(t(apply(missIn,1,rank)),2,mean)==min(apply(t(apply(missIn,1,rank)),2,mean))))
				maxIndex=min(which(apply(t(apply(missIn,1,rank)),2,mean)==min(apply(t(apply(missIn,1,rank)),2,mean))))
				# Find a middle index | Either arithmetic or geometric mean
				# midIndex=which.min(abs(feat-mean(c(feat[minIndex],feat[maxIndex]))))  # Arithmetic
				midIndex=which.min(abs(feat-exp(mean(log(c(feat[minIndex],feat[maxIndex]))))))  # Geometric
			}
			if (metric=='rmsep') {
				minIndex=max(which(apply(t(apply(rmsepIn,1,rank)),2,mean)==min(apply(t(apply(rmsepIn,1,rank)),2,mean))))
				maxIndex=min(which(apply(t(apply(rmsepIn,1,rank)),2,mean)==min(apply(t(apply(rmsepIn,1,rank)),2,mean))))
				# Find a middle index | Either arithmetic or geometric mean
				# midIndex=which.min(abs(feat-mean(c(feat[minIndex],feat[maxIndex]))))  # Arithmetic
				midIndex=which.min(abs(feat-exp(mean(log(c(feat[minIndex],feat[maxIndex]))))))  # Geometric
			}
			# Per outer segment: Average inner loop features, nComp and VIP ranks 
			featOutMin[i]=feat[minIndex]
			featOutMid[i]=feat[midIndex]
			featOutMax[i]=feat[maxIndex]
			nCompOutMin[i]=round(mean(nCompIn[,minIndex]))
			nCompOutMid[i]=round(mean(nCompIn[,midIndex]))
			nCompOutMax[i]=round(mean(nCompIn[,maxIndex]))
			VIPOutMin[,i]=apply(VIPInner[,minIndex,],1,mean)
			VIPOutMid[,i]=apply(VIPInner[,midIndex,],1,mean)
			VIPOutMax[,i]=apply(VIPInner[,maxIndex,],1,mean)
			# Build outer model for min and max nComp and predict YTEST
			xIn=X[-TestIndex,] # Perform Validation on all samples except holdout set
			yIn=Y[-TestIndex]
			incFeatMin=rownames(VIPOutMin)[rank(VIPOutMin[,i])<=featOutMin[i]]
			xTrainMin=subset(xIn,select=incFeatMin)
			plsOutMin=pls(xTrainMin,yIn,ncomp=nCompOutMin[i],mode="classic")
			xTestMin=subset(xTest,select=incFeatMin)
			yPredMinR[TestIndex]=predict(plsOutMin,newdata=xTestMin)$predict[,,nCompOutMin[i]]  # 	
			incFeatMid=rownames(VIPOutMid)[rank(VIPOutMid[,i])<=featOutMid[i]]
			xTrainMid=subset(xIn,select=incFeatMid)
			plsOutMid=pls(xTrainMid,yIn,ncomp=nCompOutMid[i],mode="classic")
			xTestMid=subset(xTest,select=incFeatMid)
			yPredMidR[TestIndex]=predict(plsOutMid,newdata=xTestMid)$predict[,,nCompOutMid[i]]  # 	
			incFeatMax=rownames(VIPOutMax)[rank(VIPOutMax[,i])<=featOutMax[i]]
			xTrainMax=subset(xIn,select=incFeatMax)
			plsOutMax=pls(xTrainMax,yIn,ncomp=nCompOutMax[i],mode="classic")
			xTestMax=subset(xTest,select=incFeatMax)
			if (length(plsOutMax$nzv$Position)>0) {
				removeFeat=rownames(plsOutMax$nzv$Metrics)
				xTestMax=xTestMax[,!colnames(xTestMax)%in%removeFeat]
			}
			yPredMaxR[TestIndex]=predict(plsOutMax,newdata=xTestMax)$predict[,,nCompOutMax[i]]  # 
			if (P==T) {  # Predict extra prediction samples (XP)
				YPR[,i]=predict(plsOutMid,newdata=subset(XP,select=incFeatMid))$predict[,,nCompOutMid[i]]
			}
		}
		# Per repetition: Average outer loop features, nComp and VIP ranks 
		featRepMinR=round(mean(featOutMin))
		nCompRepMinR=round(mean(nCompOutMin))
		VIPRepMinR=apply(VIPOutMin,1,mean)
		featRepMidR=round(mean(featOutMid))
		nCompRepMidR=round(mean(nCompOutMid))
		VIPRepMidR=apply(VIPOutMid,1,mean)
		featRepMaxR=round(mean(featOutMax))
		nCompRepMaxR=round(mean(nCompOutMax))
		VIPRepMaxR=apply(VIPOutMax,1,mean)
		parReturn=list(yPredMin=yPredMinR,featRepMin=featRepMinR,nCompRepMin=nCompRepMinR,VIPRepMin=VIPRepMinR,
			yPredMid=yPredMidR,featRepMid=featRepMidR,nCompRepMid=nCompRepMidR,VIPRepMid=VIPRepMidR,
			yPredMax=yPredMaxR,featRepMax=featRepMaxR,nCompRepMax=nCompRepMaxR,VIPRepMax=VIPRepMaxR)
		if (P==T) parReturn$YP=YPR
		reps[[r]]=parReturn
		# return(parReturn)
	}
	for (r in 1:nRep) {
		yPredMin[,r]=reps[[r]]$yPredMin
		featRepMin[r]=reps[[r]]$featRepMin
		nCompRepMin[r]=reps[[r]]$nCompRepMin
		VIPRepMin[,r]=reps[[r]]$VIPRepMin
		yPredMid[,r]=reps[[r]]$yPredMid
		featRepMid[r]=reps[[r]]$featRepMid
		nCompRepMid[r]=reps[[r]]$nCompRepMid
		VIPRepMid[,r]=reps[[r]]$VIPRepMid
		yPredMax[,r]=reps[[r]]$yPredMax
		featRepMax[r]=reps[[r]]$featRepMax
		nCompRepMax[r]=reps[[r]]$nCompRepMax
		VIPRepMax[,r]=reps[[r]]$VIPRepMax
		if (P==T) YP[,(nSeg*(r-1)+1):(nSeg*r)]=reps[[r]]$YP
	}
	yMin=yMid=yMax=matrix(nrow=nrow(X),ncol=3)
	colnames(yMin)=colnames(yMid)=colnames(yMax)=c('Mean','SD','CV')
	yMin[,1]=apply(yPredMin,1,mean)[1:nrow(X)]
	yMin[,2]=apply(yPredMin,1,sd)[1:nrow(X)]
	yMin[,3]=abs(100*yMin[,2]/yMin[,1])
	yMid[,1]=apply(yPredMid,1,mean)[1:nrow(X)]
	yMid[,2]=apply(yPredMid,1,sd)[1:nrow(X)]
	yMid[,3]=abs(100*yMid[,2]/yMid[,1])
	yMax[,1]=apply(yPredMax,1,mean)[1:nrow(X)]
	yMax[,2]=apply(yPredMax,1,sd)[1:nrow(X)]
	yMax[,3]=abs(100*yMax[,2]/yMax[,1])
	# Make ROCs
	if (DA==T) {
		ROCMin=roc(Y,yMin[,1])
		ROCMid=roc(Y,yMid[,1])
		ROCMax=roc(Y,yMax[,1])
		# Classify predictions
		yClassMin=ifelse(yMin[,1]>coords(ROCMin,'b',ret='t')[1],1,-1)
		yClassMid=ifelse(yMid[,1]>coords(ROCMid,'b',ret='t')[1],1,-1)
		yClassMax=ifelse(yMax[,1]>coords(ROCMax,'b',ret='t')[1],1,-1)
		# Bind together all Y data
		yMat=cbind(yMin[,1],yMid[,1],yMax[,1],yClassMin,yClassMax,Y)
		colnames(yMat)[1:3]=c('yMin','yMid','yMax')
		# Calculate misclassifications
		missMin=sum(abs((Y-yClassMin))/2)
		missMid=sum(abs((Y-yClassMid))/2)
		missMax=sum(abs((Y-yClassMax))/2)
	}
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
	# Average nComp over repetitions
	nCompMin=c(mean(nCompRepMin),sd(nCompRepMin),abs(100*sd(nCompRepMin)/mean(nCompRepMin)))
	nCompMid=c(mean(nCompRepMid),sd(nCompRepMid),abs(100*sd(nCompRepMid)/mean(nCompRepMid)))
	nCompMax=c(mean(nCompRepMax),sd(nCompRepMax),abs(100*sd(nCompRepMax)/mean(nCompRepMax)))
	names(nCompMin)=names(nCompMid)=names(nCompMax)=c('Mean','SD','CV')
	# Average nFeat over repetitions
	nFeatMin=c(mean(featRepMin),sd(featRepMin),abs(100*sd(featRepMin)/mean(featRepMin)))
	nFeatMid=c(mean(featRepMid),sd(featRepMid),abs(100*sd(featRepMid)/mean(featRepMid)))
	nFeatMax=c(mean(featRepMax),sd(featRepMax),abs(100*sd(featRepMax)/mean(featRepMax)))
	names(nFeatMin)=names(nFeatMid)=names(nFeatMax)=c('Mean','SD','CV')
	# Stop timer
	end.time=proc.time()[3]
	# Generate report-list
	modelReturn$nFeat=list(minModel=nFeatMin,midModel=nFeatMid,maxModel=nFeatMax)
	modelReturn$nComp=list(minModel=nCompMin,midModel=nCompMid,maxModel=nCompMax)
	modelReturn$nFeatPerRep=list(minModel=featRepMin,midModel=featRepMid,maxModel=featRepMax)
	modelReturn$nCompPerRep=list(minModel=nCompRepMin,midModel=nCompRepMid,maxModel=nCompRepMax)
	modelReturn$yPred=list(minModel=yMin,midModel=yMid,maxModel=yMax)
	modelReturn$VIPRank=list(minModel=VIPMin,midModel=VIPMid,maxModel=VIPMax)
	if (P==T) modelReturn$YP=YP
	if (DA==T) {
		modelReturn$auc=list(minModel=ROCMin$auc,midModel=ROCMid$auc,maxModel=ROCMax$auc)
		modelReturn$yClass=list(minModel=yClassMin[1:nrow(X)],midModel=yClassMid[1:nrow(X)],maxModel=yClassMax[1:nrow(X)])
		modelReturn$misClass=list(minModel=missMin,midModel=missMid,maxModel=missMax)
	}
	modelReturn$calcMins=(end.time-start.time)/60
	cat('\n Elapsed time',(end.time-start.time)/60,'mins \n')
	return(modelReturn)
}

### parallelized rdCV
ryeCV4Par=function(X,Y,ID,nRep=10,nSeg=7,nInner=6,comps=5,metric=c('miss','auc','rmsep'),DA=FALSE,featRatio=0.75,P=FALSE,XP=NA) {  # Perform rdCV PLS 
	start.time=proc.time()[3]
	if (missing(metric)) metric='rmsep'
	if (P==T) {  # If model is also used to predict XP matrix, allocate a YP response array
		YP=matrix(nrow=nrow(XP),ncol=nSeg*nRep)  # Predictions for all samples (rows) across outer segments (cols) and repetitions (dim3)
	}
	modelReturn=list()
	nFeat=nFeat0=ncol(X)
	comp=ifelse(nFeat<comps,nFeat,comps)
	unik=!duplicated(ID)
	IDunik=ID[unik]
	Yunik=Y[unik]
	YIndex=1:length(Y)
	ind=1:length(Yunik) # Make indices of unique individuals INSTEAD OF OBSERVATIONS!!!
	if (DA==T) {
		name=unique(Y)  # Find groups
		groups=length(name) # Number of groups
		groupInd=list()  # Allocate list for indices of groups
		for (g in 1:groups) { 
			groupInd[[g]]=ind[Yunik==name[g]]  # Find indices per group
		}
	}
	VIPMin=VIPMid=VIPMax=matrix(nrow=nFeat0,ncol=3) # Allocate matrix for variable importance averaged over repetition
	rownames(VIPMin)=rownames(VIPMid)=rownames(VIPMax)=colnames(X)
	colnames(VIPMin)=colnames(VIPMid)=colnames(VIPMax)=c('Mean','SD','CV')
	# Allocate matrix for prediction of outer segments Y per repetition
	yPredMin=yPredMid=yPredMax=matrix(nrow=length(Y),ncol=nRep)
	colnames(yPredMin)=colnames(yPredMid)=colnames(yPredMax)=paste('Rep',1:nRep,sep='')
	yPredMinR=yPredMidR=yPredMaxR=numeric(length(Y))
	# Allocate response vectors and matrices for feat's, nComp and VIP ranks over repetitions
	featRepMin=featRepMid=featRepMax=nCompRepMin=nCompRepMid=nCompRepMax=missRep=numeric(nRep)
	names(featRepMin)=names(featRepMid)=names(featRepMax)=names(nCompRepMin)=names(nCompRepMid)=names(nCompRepMax)=names(missRep)=paste(rep('rep',nRep),1:nRep,sep='')
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
	# 
	# reps=list()
	# for (r in 1:nRep) {  # Loop for all repetitions
	reps=foreach(r=1:nRep, .packages=c('mixOmics','pROC'), .export='vectSamp') %dopar% {
		# r=1
		if (P==T) YPR=matrix(nrow=nrow(XP),ncol=nSeg)
		cat('\n','   Repetition ',r,' of ',nRep,':',sep='')
		if (DA==T) {
			groupTest=list()  ## Allocate list for samples within group
			for (gT in 1:groups) { 
				groupTest[[gT]]=vectSamp(groupInd[[gT]],n=nSeg)  # Draw random samples within group
			}
			groupTest[[2]]=rev(groupTest[[2]])  ## Reverse order for 2nd group 
			## In the case of >2 groups another approach would be suitable
			## eg adding subsequent groups to allTest ordered by length
			allTest=groupTest[[1]] # Add 1st groups to 'Master' sample of all groups
			for (gT in 2:groups) {  # Add subsequent groups
				for (aT in 1:nSeg) {
					allTest[[aT]]=sort(c(allTest[[aT]],groupTest[[gT]][[aT]]))
				}
			}
		} else {
			allTest=vectSamp(ind,n=nSeg)
		}
		featOutMin=featOutMid=featOutMax=nCompOutMin=nCompOutMid=nCompOutMax=numeric(nSeg)
		names(featOutMin)=names(featOutMid)=names(featOutMax)=names(nCompOutMin)=names(nCompOutMax)=paste(rep('outSeg',nSeg),1:nSeg,sep='')
		VIPOutMin=VIPOutMid=VIPOutMax=matrix(data=nFeat0,nrow=nFeat0,ncol=nSeg)
		rownames(VIPOutMin)=rownames(VIPOutMid)=rownames(VIPOutMax)=colnames(X)
		colnames(VIPOutMin)=colnames(VIPOutMid)=colnames(VIPOutMax)=paste(rep('outSeg',nSeg),1:nSeg,sep='')
		for (i in 1:nSeg) {   # Create 'nSamp' models within repetition
			# i=1
			cat('\n Segment ',i,' (features):',sep='') # Counter
			testInd=allTest[[i]] # Draw out segment = holdout set BASED ON UNIK
			TestIndex=numeric()  # BASED ON OBSERVATIONS
			for (xt in 1:length(testInd)) {
				# cat(sum(ID==IDunik[testInd[xt]]))
				TestIndex=c(TestIndex,which(ID==IDunik[testInd[xt]]))
			}
			xTest=X[TestIndex,]
			yTest=Y[TestIndex]
			inInd=ind[-testInd]
			InnerIndex=YIndex[-TestIndex]
			yUnikIn=Yunik[inInd]  
			###
			
			# cat('\n',head(testInd))   # Trouble shooting
			# cat('\n',head(yTrain))
			# cat('\n',as.numeric(xTrain[1,1:5]))
			missIn=aucIn=rmsepIn=nCompIn=matrix(nrow=nInner,ncol=cnt)
			rownames(rmsepIn)=rownames(missIn)=rownames(aucIn)=rownames(nCompIn)=paste(rep('inSeg',nInner),1:nInner,sep='')
			colnames(rmsepIn)=colnames(missIn)=colnames(aucIn)=colnames(nCompIn)=feat
			VIPInner=array(data=nFeat0,dim=c(nFeat0,cnt,nInner))
			rownames(VIPInner)=colnames(X)
			colnames(VIPInner)=feat
			dimnames(VIPInner)[[3]]=paste(rep('inSeg',nInner),1:nInner,sep='')
			nFeat=nFeat0
			# inSamp=vectSamp(inInd,n=nInner)  # Draw random samples within current 2CV segment
			incFeat=colnames(X)
			for (count in 1:cnt) {  # Build models with successively fewer feature. Quality metric = number of missclassifications for Validation set
				# count=1 # for testing
				# count=count+1
				nFeat=feat[count]
				cat(nFeat)
				# if (resampInner==TRUE) inSamp=vectSamp(inInd,n=nInner)  # Resample inner segments for each feature elimination step
				comp=ifelse(nFeat<comps,nFeat,comps)
				if (DA==T) {
					groupIndVal=list()
					for (g in 1:groups) { 
						groupIndVal[[g]]=inInd[yUnikIn==name[g]]  # Find indices per group
					}
					groupVal=list()  ## Allocate list for samples within group
					for (gV in 1:groups) { 
						groupVal[[gV]]=vectSamp(groupIndVal[[gV]],n=nInner)  # Draw random samples within group
					}
					groupVal[[2]]=rev(groupVal[[2]])  ## Reverse order for 2nd group 
					## In the case of >2 groups another approach would be suitable
					## eg adding subsequent groups to allVal ordered by length
					allVal=groupVal[[1]] # Add 1st groups to 'Master' sample of all groups
					for (gV in 2:groups) {  # Add subsequent groups
						for (aV in 1:nInner) {
							allVal[[aV]]=sort(c(allVal[[aV]],groupVal[[gV]][[aV]]))
						}
					}
				} else {
					allVal=vectSamp(inInd,n=nInner)
				}
				for (j in 1:nInner) {
					# j=1 # for testing
					# j=j+1
					cat('.') # Counter
					valInd=allVal[[j]] # Draw out segment = validation set
					ValIndex=numeric()  # BASED ON OBSERVATIONS
					for (xv in 1:length(valInd)) {
						# cat(sum(ID==IDunik[valInd[xv]]))
						ValIndex=c(ValIndex,which(ID==IDunik[valInd[xv]]))
					}
					xVal=X[ValIndex,]
					xVal=subset(xVal,select=incFeat)
					yVal=Y[ValIndex]
					trainInd=InnerIndex[-match(ValIndex,InnerIndex)] # Define Training segment
					xTrain=X[trainInd,]
					xTrain=subset(xTrain,select=incFeat)
					yTrain=Y[trainInd]
					# perform PLS up to 'comps' number of components
					plsInner=pls(xTrain,yTrain,ncomp=comp,mode="regression")
					if (length(plsInner$nzv$Position)>0) {
						removeFeat=rownames(plsInner$nzv$Metrics)
						xVal=xVal[,!colnames(xVal)%in%removeFeat]
					}
					yValInner=predict(plsInner,newdata=xVal)$predict[,,]  # Store  prediction estimates per validation segment 
					if (metric=='miss') {
						# cat(' miss',count)
						yClassInner=ifelse(yValInner>0,1,-1)
						misClass=apply((yVal-yClassInner)*yVal/2,2,sum,na.rm=T)
						missIn[j,count]=min(misClass)
						nCompIn[j,count]=which.min(misClass)
					} 
					if (metric=='auc') {
						# cat(' auc',count)
						auc=apply(yValInner,2,function(x) roc(yVal,x)$auc)
						aucIn[j,count]=max(auc)
						nCompIn[j,count]=which.max(auc)
					}
					if (metric=='rmsep') {
						# cat(' rmsep',count)
						rmsep=apply(yValInner,2,function(x) sqrt(sum((yVal-x)^2,na.rm=T)/(length(yValInner[,1])-sum(is.na(yValInner[,1])))))
						rmsepIn[j,count]=min(rmsep)
						nCompIn[j,count]=which.min(rmsep)
					}
					VIPInner[match(names(vip(plsInner)[,nCompIn[j,count]]),rownames(VIPInner)),count,j]=rank(-vip(plsInner)[,nCompIn[j,count]])
					# VIPInner[match(names(vip(plsInner)[,nCompIn[j]]),rownames(VIPInner)),j]=rank(-vip(plsInner)[,nCompIn[j]])				
				}
				## NB!!! Average VIP ranks over inner segments before feature elimination!!!
				VIPInAve=apply(VIPInner[,count,],1,mean)
				if (count<cnt) {
					incFeat=names(VIPInAve[order(VIPInAve)])[1:feat[count+1]]
				}
			}
			if (metric=='auc') {
				minIndex=max(which(apply(t(apply(aucIn,1,rank)),2,mean)==max(apply(t(apply(aucIn,1,rank)),2,mean))))
				maxIndex=min(which(apply(t(apply(aucIn,1,rank)),2,mean)==max(apply(t(apply(aucIn,1,rank)),2,mean))))
				# Find a middle index | Either arithmetic or geometric mean
				# midIndex=which.min(abs(feat-mean(c(feat[minIndex],feat[maxIndex]))))  # Arithmetic
				midIndex=which.min(abs(feat-exp(mean(log(c(feat[minIndex],feat[maxIndex]))))))  # Geometric
			}
			if (metric=='miss') {
				minIndex=max(which(apply(t(apply(missIn,1,rank)),2,mean)==min(apply(t(apply(missIn,1,rank)),2,mean))))
				maxIndex=min(which(apply(t(apply(missIn,1,rank)),2,mean)==min(apply(t(apply(missIn,1,rank)),2,mean))))
				# Find a middle index | Either arithmetic or geometric mean
				# midIndex=which.min(abs(feat-mean(c(feat[minIndex],feat[maxIndex]))))  # Arithmetic
				midIndex=which.min(abs(feat-exp(mean(log(c(feat[minIndex],feat[maxIndex]))))))  # Geometric
			}
			if (metric=='rmsep') {
				minIndex=max(which(apply(t(apply(rmsepIn,1,rank)),2,mean)==min(apply(t(apply(rmsepIn,1,rank)),2,mean))))
				maxIndex=min(which(apply(t(apply(rmsepIn,1,rank)),2,mean)==min(apply(t(apply(rmsepIn,1,rank)),2,mean))))
				# Find a middle index | Either arithmetic or geometric mean
				# midIndex=which.min(abs(feat-mean(c(feat[minIndex],feat[maxIndex]))))  # Arithmetic
				midIndex=which.min(abs(feat-exp(mean(log(c(feat[minIndex],feat[maxIndex]))))))  # Geometric
			}
			# Per outer segment: Average inner loop features, nComp and VIP ranks 
			featOutMin[i]=feat[minIndex]
			featOutMid[i]=feat[midIndex]
			featOutMax[i]=feat[maxIndex]
			nCompOutMin[i]=round(mean(nCompIn[,minIndex]))
			nCompOutMid[i]=round(mean(nCompIn[,midIndex]))
			nCompOutMax[i]=round(mean(nCompIn[,maxIndex]))
			VIPOutMin[,i]=apply(VIPInner[,minIndex,],1,mean)
			VIPOutMid[,i]=apply(VIPInner[,midIndex,],1,mean)
			VIPOutMax[,i]=apply(VIPInner[,maxIndex,],1,mean)
			# Build outer model for min and max nComp and predict YTEST
			xIn=X[-TestIndex,] # Perform Validation on all samples except holdout set
			yIn=Y[-TestIndex]
			incFeatMin=rownames(VIPOutMin)[rank(VIPOutMin[,i])<=featOutMin[i]]
			xTrainMin=subset(xIn,select=incFeatMin)
			plsOutMin=pls(xTrainMin,yIn,ncomp=nCompOutMin[i],mode="classic")
			xTestMin=subset(xTest,select=incFeatMin)
			yPredMinR[TestIndex]=predict(plsOutMin,newdata=xTestMin)$predict[,,nCompOutMin[i]]  # 	
			incFeatMid=rownames(VIPOutMid)[rank(VIPOutMid[,i])<=featOutMid[i]]
			xTrainMid=subset(xIn,select=incFeatMid)
			plsOutMid=pls(xTrainMid,yIn,ncomp=nCompOutMid[i],mode="classic")
			xTestMid=subset(xTest,select=incFeatMid)
			yPredMidR[TestIndex]=predict(plsOutMid,newdata=xTestMid)$predict[,,nCompOutMid[i]]  # 	
			incFeatMax=rownames(VIPOutMax)[rank(VIPOutMax[,i])<=featOutMax[i]]
			xTrainMax=subset(xIn,select=incFeatMax)
			plsOutMax=pls(xTrainMax,yIn,ncomp=nCompOutMax[i],mode="classic")
			xTestMax=subset(xTest,select=incFeatMax)
			if (length(plsOutMax$nzv$Position)>0) {
				removeFeat=rownames(plsOutMax$nzv$Metrics)
				xTestMax=xTestMax[,!colnames(xTestMax)%in%removeFeat]
			}
			yPredMaxR[TestIndex]=predict(plsOutMax,newdata=xTestMax)$predict[,,nCompOutMax[i]]  # 
			if (P==T) {  # Predict extra prediction samples (XP)
				YPR[,i]=predict(plsOutMid,newdata=subset(XP,select=incFeatMid))$predict[,,nCompOutMid[i]]
			}
		}
		# Per repetition: Average outer loop features, nComp and VIP ranks 
		featRepMinR=round(mean(featOutMin))
		nCompRepMinR=round(mean(nCompOutMin))
		VIPRepMinR=apply(VIPOutMin,1,mean)
		featRepMidR=round(mean(featOutMid))
		nCompRepMidR=round(mean(nCompOutMid))
		VIPRepMidR=apply(VIPOutMid,1,mean)
		featRepMaxR=round(mean(featOutMax))
		nCompRepMaxR=round(mean(nCompOutMax))
		VIPRepMaxR=apply(VIPOutMax,1,mean)
		parReturn=list(yPredMin=yPredMinR,featRepMin=featRepMinR,nCompRepMin=nCompRepMinR,VIPRepMin=VIPRepMinR,
			yPredMid=yPredMidR,featRepMid=featRepMidR,nCompRepMid=nCompRepMidR,VIPRepMid=VIPRepMidR,
			yPredMax=yPredMaxR,featRepMax=featRepMaxR,nCompRepMax=nCompRepMaxR,VIPRepMax=VIPRepMaxR)
		if (P==T) parReturn$YP=YPR
		# reps[[r]]=parReturn
		return(parReturn)
	}
	for (r in 1:nRep) {
		yPredMin[,r]=reps[[r]]$yPredMin
		featRepMin[r]=reps[[r]]$featRepMin
		nCompRepMin[r]=reps[[r]]$nCompRepMin
		VIPRepMin[,r]=reps[[r]]$VIPRepMin
		yPredMid[,r]=reps[[r]]$yPredMid
		featRepMid[r]=reps[[r]]$featRepMid
		nCompRepMid[r]=reps[[r]]$nCompRepMid
		VIPRepMid[,r]=reps[[r]]$VIPRepMid
		yPredMax[,r]=reps[[r]]$yPredMax
		featRepMax[r]=reps[[r]]$featRepMax
		nCompRepMax[r]=reps[[r]]$nCompRepMax
		VIPRepMax[,r]=reps[[r]]$VIPRepMax
		if (P==T) YP[,(nSeg*(r-1)+1):(nSeg*r)]=reps[[r]]$YP
	}
	yMin=yMid=yMax=matrix(nrow=nrow(X),ncol=3)
	colnames(yMin)=colnames(yMid)=colnames(yMax)=c('Mean','SD','CV')
	yMin[,1]=apply(yPredMin,1,mean)[1:nrow(X)]
	yMin[,2]=apply(yPredMin,1,sd)[1:nrow(X)]
	yMin[,3]=abs(100*yMin[,2]/yMin[,1])
	yMid[,1]=apply(yPredMid,1,mean)[1:nrow(X)]
	yMid[,2]=apply(yPredMid,1,sd)[1:nrow(X)]
	yMid[,3]=abs(100*yMid[,2]/yMid[,1])
	yMax[,1]=apply(yPredMax,1,mean)[1:nrow(X)]
	yMax[,2]=apply(yPredMax,1,sd)[1:nrow(X)]
	yMax[,3]=abs(100*yMax[,2]/yMax[,1])
	# Make ROCs
	if (DA==T) {
		ROCMin=roc(Y,yMin[,1])
		ROCMid=roc(Y,yMid[,1])
		ROCMax=roc(Y,yMax[,1])
		# Classify predictions
		yClassMin=ifelse(yMin[,1]>coords(ROCMin,'b',ret='t')[1],1,-1)
		yClassMid=ifelse(yMid[,1]>coords(ROCMid,'b',ret='t')[1],1,-1)
		yClassMax=ifelse(yMax[,1]>coords(ROCMax,'b',ret='t')[1],1,-1)
		# Bind together all Y data
		yMat=cbind(yMin[,1],yMid[,1],yMax[,1],yClassMin,yClassMax,Y)
		colnames(yMat)[1:3]=c('yMin','yMid','yMax')
		# Calculate misclassifications
		missMin=sum(abs((Y-yClassMin))/2)
		missMid=sum(abs((Y-yClassMid))/2)
		missMax=sum(abs((Y-yClassMax))/2)
	}
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
	# Average nComp over repetitions
	nCompMin=c(mean(nCompRepMin),sd(nCompRepMin),abs(100*sd(nCompRepMin)/mean(nCompRepMin)))
	nCompMid=c(mean(nCompRepMid),sd(nCompRepMid),abs(100*sd(nCompRepMid)/mean(nCompRepMid)))
	nCompMax=c(mean(nCompRepMax),sd(nCompRepMax),abs(100*sd(nCompRepMax)/mean(nCompRepMax)))
	names(nCompMin)=names(nCompMid)=names(nCompMax)=c('Mean','SD','CV')
	# Average nFeat over repetitions
	nFeatMin=c(mean(featRepMin),sd(featRepMin),abs(100*sd(featRepMin)/mean(featRepMin)))
	nFeatMid=c(mean(featRepMid),sd(featRepMid),abs(100*sd(featRepMid)/mean(featRepMid)))
	nFeatMax=c(mean(featRepMax),sd(featRepMax),abs(100*sd(featRepMax)/mean(featRepMax)))
	names(nFeatMin)=names(nFeatMid)=names(nFeatMax)=c('Mean','SD','CV')
	# Stop timer
	end.time=proc.time()[3]
	# Generate report-list
	modelReturn$nFeat=list(minModel=nFeatMin,midModel=nFeatMid,maxModel=nFeatMax)
	modelReturn$nComp=list(minModel=nCompMin,midModel=nCompMid,maxModel=nCompMax)
	modelReturn$nFeatPerRep=list(minModel=featRepMin,midModel=featRepMid,maxModel=featRepMax)
	modelReturn$nCompPerRep=list(minModel=nCompRepMin,midModel=nCompRepMid,maxModel=nCompRepMax)
	modelReturn$yPred=list(minModel=yMin,midModel=yMid,maxModel=yMax)
	modelReturn$VIPRank=list(minModel=VIPMin,midModel=VIPMid,maxModel=VIPMax)
	if (P==T) modelReturn$YP=YP
	if (DA==T) {
		modelReturn$auc=list(minModel=ROCMin$auc,midModel=ROCMid$auc,maxModel=ROCMax$auc)
		modelReturn$yClass=list(minModel=yClassMin[1:nrow(X)],midModel=yClassMid[1:nrow(X)],maxModel=yClassMax[1:nrow(X)])
		modelReturn$misClass=list(minModel=missMin,midModel=missMid,maxModel=missMax)
	}
	modelReturn$calcMins=(end.time-start.time)/60
	cat('\n Elapsed time',(end.time-start.time)/60,'mins \n')
	return(modelReturn)
}


### Setup parallel processing
cl=makeCluster(detectCores())
registerDoParallel(cl)

ryeModel2=ryeCVPar2(X,Y,ID,100,comps=3)
save(ryeModel2,file='ryeModel2.Rdata')
load(file='ryeModel2.Rdata')
### Stop parallel processing
stopCluster(cl)

##############################
# LOOP over successively fewer features
# eliminate highest rank
# take out auroc value at each step

library(doParallel)
library(foreach)
# library(mixOmics)  ## These are called within parallel rdCV
# library(pROC)

### Load data
# ryePath=file.choose()
# ryePath="\\\\storage-ua.slu.se\\home$\\brunius\\My documents\\Papers\\WG_metab\\lcms_rye_da.csv"
# ryeDA=read.csv2(file=ryePath)
# save(ryeDA,file='ryeDA.RData')
load(file='ryeDA.RData')
colnames(ryeDA)[4]='ind'
Y=ryeDA[,3]
ID=ryeDA[,4]
X=ryeDA[,-1:-5]

nRep=1
nSeg=7
nInner=6

### Setup parallel processing
cl=makeCluster(detectCores())
registerDoParallel(cl)

### Perform analysis
nFeat=ncol(X)
Xmod=X
AUC=featRemove=numeric(nFeat)
for (feat in nFeat:2) {
	cat('\n',feat)
	ryeModel=ryeCVPar(Xmod,Y,ID,20,comps=3)
	AUC[feat]=ryeModel$yPredAuc
	featRemove[feat]=names(which.max(ryeModel$VIRankAve))
	Xmod=Xmod[,colnames(Xmod)!=featRemove[feat]]
}

### Stop parallel processing
stopCluster(cl)

### Sort out some variables
featRemove[1]=colnames(X)[!colnames(X) %in% featRemove]
ROC=apply(X,2,function(x) roc(Y,x)$auc)
ROC=ROC[order(names(ROC))]
ROC=ROC[rank(featRemove)]
tab=cbind(AUC,ROC)
write.csv(tab,file='tabPar.csv')

### Make nice plot
t=1:30
png(file='auroc_Par.png',width=1024, height=756, pointsize=20)
plot(AUC,ylim=c(0.5,1),axes=F,col=1,xlab='Variable index',ylab='AUROC')
axis(side=1)
axis(side=1,1:30,label=F,tck=-0.01)
axis(side = 2, las = 1)
spl=smooth.spline(t[-1],AUC[-1],spar=.4)
lines(spl)
points(ROC,pch=2,col=2)
legend(x=5,y=.7,legend=c('PLS model','Univariate'),pch=1:2,col=1:2,lty=c(1,F))
dev.off()