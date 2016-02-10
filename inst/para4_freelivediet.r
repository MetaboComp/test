X=matrix(1:6,ncol=2)
Y=1:3

### parallelized rdCV
ryeCV4Par=function(X,Y,ID,nRep=10,nOuter=7,nInner=6,comps=5,metric=c('miss','auc','rmsep'),DA=FALSE,featRatio=0.75,pred=FALSE,XP=NA) {  # Perform rdCV PLS 
	start.time=proc.time()[3]
	if (missing(metric)) metric='rmsep'
	if (pred) {  # If model is also used to predict XP matrix, allocate a YP response array
		YP=matrix(nrow=nrow(XP),ncol=nOuter*nRep)  # Predictions for all samples (rows) across outer segments (cols) and repetitions (dim3)
	}
	modelReturn=list()
	nVar=nVar0=ncol(X)
	comp=ifelse(nVar<comps,nVar,comps)
	
	# This will probably be done differently!
      	unik=!duplicated(ID)
      	unikID=ID[unik]
      	unikY=Y[unik]
      	YIndex=1:length(Y)
      	ind=1:length(unikY) # Make indices of unique individuals INSTEAD OF OBSERVATIONS!!!

  if (DA) {
		Ynames=unique(Y)  # Find groups
		groups=length(Ynames) # Number of groups
		groupInd=list()  # Allocate list for indices of groups
		for (g in 1:groups) { 
			groupInd[[g]]=ind[unikY==Ynames[g]]  # Find indices per group
		}
	}
	VIPMin=VIPMid=VIPMax=matrix(nrow=nVar0,ncol=3) # Allocate matrix for variable importance averaged over repetition
	rownames(VIPMin)=rownames(VIPMid)=rownames(VIPMax)=colnames(X)
	colnames(VIPMin)=colnames(VIPMid)=colnames(VIPMax)=c('Mean','SD','CV')
	# Allocate matrix for prediction of outer segments Y per repetition
	yPredMin=yPredMid=yPredMax=matrix(nrow=length(Y),ncol=nRep)
	colnames(yPredMin)=colnames(yPredMid)=colnames(yPredMax)=paste('Rep',1:nRep,sep='')
	yPredMinR=yPredMidR=yPredMaxR=numeric(length(Y))
	# Allocate response vectors and matrices for feat's, nComp and VIP ranks over repetitions
	featRepMin=featRepMid=featRepMax=nCompRepMin=nCompRepMid=nCompRepMax=missRep=numeric(nRep)
	names(featRepMin)=names(featRepMid)=names(featRepMax)=names(nCompRepMin)=names(nCompRepMid)=names(nCompRepMax)=names(missRep)=paste(rep('rep',nRep),1:nRep,sep='')
	VIPRepMin=VIPRepMid=VIPRepMax=matrix(data=nVar0,nrow=nVar0,ncol=nRep)
	rownames(VIPRepMin)=rownames(VIPRepMid)=rownames(VIPRepMax)=colnames(X)
	colnames(VIPRepMin)=colnames(VIPRepMid)=colnames(VIPRepMax)=paste(rep('rep',nRep),1:nRep,sep='')
	feat=numeric()
	cnt=0
	while (nVar>2) {  
		cnt=cnt+1
		feat=c(feat,nVar)
		nVar=floor(featRatio*nVar)
	}
	# Perform rdCV 
	# 
	# reps=list()
	# for (r in 1:nRep) {  # Loop for all repetitions
	reps=foreach(r=1:nRep, .packages=c('mixOmics','pROC'), .export='vectSamp') %dopar% {
		# r=1
		if (pred) YPR=matrix(nrow=nrow(XP),ncol=nOuter)
		cat('\n','   Repetition ',r,' of ',nRep,':',sep='')
		if (DA) {
			groupTest=list()  ## Allocate list for samples within group
			for (gT in 1:groups) { 
				groupTest[[gT]]=vectSamp(groupInd[[gT]],n=nOuter)  # Draw random samples within group
			}
			groupTest[[2]]=rev(groupTest[[2]])  ## Reverse order for 2nd group 
			## In the case of >2 groups another approach would be suitable
			## eg adding subsequent groups to allTest ordered by length
			allTest=groupTest[[1]] # Add 1st groups to 'Master' sample of all groups
			for (gT in 2:groups) {  # Add subsequent groups
				for (aT in 1:nOuter) {
					allTest[[aT]]=sort(c(allTest[[aT]],groupTest[[gT]][[aT]]))
				}
			}
		} else {
			allTest=vectSamp(ind,n=nOuter)
		}
		featOutMin=featOutMid=featOutMax=nCompOutMin=nCompOutMid=nCompOutMax=numeric(nOuter)
		names(featOutMin)=names(featOutMid)=names(featOutMax)=names(nCompOutMin)=names(nCompOutMax)=paste(rep('outSeg',nOuter),1:nOuter,sep='')
		VIPOutMin=VIPOutMid=VIPOutMax=matrix(data=nVar0,nrow=nVar0,ncol=nOuter)
		rownames(VIPOutMin)=rownames(VIPOutMid)=rownames(VIPOutMax)=colnames(X)
		colnames(VIPOutMin)=colnames(VIPOutMid)=colnames(VIPOutMax)=paste(rep('outSeg',nOuter),1:nOuter,sep='')
		for (i in 1:nOuter) {   # Create 'nSamp' models within repetition
			# i=1
			cat('\n Segment ',i,' (features):',sep='') # Counter
			testInd=allTest[[i]] # Draw out segment = holdout set BASED ON UNIK
			TestIndex=numeric()  # BASED ON OBSERVATIONS
			for (xt in 1:length(testInd)) {
				# cat(sum(ID==unikID[testInd[xt]]))
				TestIndex=c(TestIndex,which(ID==unikID[testInd[xt]]))
			}
			xTest=X[TestIndex,]
			yTest=Y[TestIndex]
			inInd=ind[-testInd]
			InnerIndex=YIndex[-TestIndex]
			yUnikIn=unikY[inInd]  
			###
			
			# cat('\n',head(testInd))   # Trouble shooting
			# cat('\n',head(yTrain))
			# cat('\n',as.numeric(xTrain[1,1:5]))
			missIn=aucIn=rmsepIn=nCompIn=matrix(nrow=nInner,ncol=cnt)
			rownames(rmsepIn)=rownames(missIn)=rownames(aucIn)=rownames(nCompIn)=paste(rep('inSeg',nInner),1:nInner,sep='')
			colnames(rmsepIn)=colnames(missIn)=colnames(aucIn)=colnames(nCompIn)=feat
			VIPInner=array(data=nVar0,dim=c(nVar0,cnt,nInner))
			rownames(VIPInner)=colnames(X)
			colnames(VIPInner)=feat
			dimnames(VIPInner)[[3]]=paste(rep('inSeg',nInner),1:nInner,sep='')
			nVar=nVar0
			# inSamp=vectSamp(inInd,n=nInner)  # Draw random samples within current 2CV segment
			incFeat=colnames(X)
			for (count in 1:cnt) {  # Build models with successively fewer feature. Quality metric = number of missclassifications for Validation set
				# count=1 # for testing
				# count=count+1
				nVar=feat[count]
				cat(nVar)
				# if (resampInner==TRUE) inSamp=vectSamp(inInd,n=nInner)  # Resample inner segments for each feature elimination step
				comp=ifelse(nVar<comps,nVar,comps)
				if (DA) {
					groupIndVal=list()
					for (g in 1:groups) { 
						groupIndVal[[g]]=inInd[yUnikIn==Ynames[g]]  # Find indices per group
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
						# cat(sum(ID==unikID[valInd[xv]]))
						ValIndex=c(ValIndex,which(ID==unikID[valInd[xv]]))
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
			if (pred) {  # Predict extra prediction samples (XP)
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
		if (pred) parReturn$YP=YPR
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
		if (pred) YP[,(nOuter*(r-1)+1):(nOuter*r)]=reps[[r]]$YP
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
	if (DA) {
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
	# Average nVar over repetitions
	nVarMin=c(mean(featRepMin),sd(featRepMin),abs(100*sd(featRepMin)/mean(featRepMin)))
	nVarMid=c(mean(featRepMid),sd(featRepMid),abs(100*sd(featRepMid)/mean(featRepMid)))
	nVarMax=c(mean(featRepMax),sd(featRepMax),abs(100*sd(featRepMax)/mean(featRepMax)))
	names(nVarMin)=names(nVarMid)=names(nVarMax)=c('Mean','SD','CV')
	# Stop timer
	end.time=proc.time()[3]
	# Generate report-list
	modelReturn$nVar=list(minModel=nVarMin,midModel=nVarMid,maxModel=nVarMax)
	modelReturn$nComp=list(minModel=nCompMin,midModel=nCompMid,maxModel=nCompMax)
	modelReturn$nVarPerRep=list(minModel=featRepMin,midModel=featRepMid,maxModel=featRepMax)
	modelReturn$nCompPerRep=list(minModel=nCompRepMin,midModel=nCompRepMid,maxModel=nCompRepMax)
	modelReturn$yPred=list(minModel=yMin,midModel=yMid,maxModel=yMax)
	modelReturn$VIPRank=list(minModel=VIPMin,midModel=VIPMid,maxModel=VIPMax)
	if (pred) modelReturn$YP=YP
	if (DA) {
		modelReturn$auc=list(minModel=ROCMin$auc,midModel=ROCMid$auc,maxModel=ROCMax$auc)
		modelReturn$yClass=list(minModel=yClassMin[1:nrow(X)],midModel=yClassMid[1:nrow(X)],maxModel=yClassMax[1:nrow(X)])
		modelReturn$misClass=list(minModel=missMin,midModel=missMid,maxModel=missMax)
	}
	modelReturn$calcMins=(end.time-start.time)/60
	cat('\n Elapsed time',(end.time-start.time)/60,'mins \n')
	return(modelReturn)
}
