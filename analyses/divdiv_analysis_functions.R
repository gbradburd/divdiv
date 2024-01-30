divdivAnalyses <- function(z,X,y,phyStr,mod,outName,nIter=1e4,parallel=TRUE,nNodes=4,filterKeep=NULL){
    `%d%` <- parallelizing(args <- as.list(environment()))
    export <- c("divdivAnalysis","makeDB","makeMultiDB")
    outs <- foreach::foreach(i = 1:length(X),
    						.export=export,
    						.packages="rstan") %d% {
								divdivAnalysis(z=z,x=X[[i]],y=y,phyStr=phyStr,
												mod=mod,nIter=nIter,filterKeep=filterKeep)
			}
	doParallel::stopImplicitCluster()
	save(outs,file=paste0(outName,"_outs.Robj"))
	return(outs)
}

divdivAnalysis <- function(z,x,y,phyStr,mod,nIter=1e4,filterKeep=NULL){
    if(inherits(mod,"character")){
	    	mod <- stan_model(model_code=mod)
    }
	if(!is.null(filterKeep)){
		z <- z[filterKeep,]
		phyStr <- phyStr[filterKeep,filterKeep]
	}
	if(length(x)==1){
		db <- makeDB(X=z[[x]],Y=z[[y]],phyStr=phyStr)
	} else {
		preds <- lapply(1:length(x),function(i){z[[x[i]]]})
		db <- makeMultiDB(X=preds,Y=z[[y]],phyStr=phyStr)		
	}
	fit <- rstan::sampling(object=mod,
							data=db,
							iter=nIter,
							thin=nIter/500,
							chains=1,
							control=setNames(list(15),"max_treedepth"))
	out <- list("db"=db,"fit"=fit)
	return(out)
}

getPlotDims <- function(N){
	n <- round(sqrt(N))
	if(N <= n^2){
		dims <- c(n,n)
	} else {
		dims <- c(n,n+1)
	}
	return(dims)
}

getSampCols <- function(db,sampCols){
	sp <- row.names(db$relMat)
	sampCols <- sampCols[names(sampCols) %in% sp]
	return(sampCols)
}

vizAllOuts <- function(outs,predNames,sampPhy,outName,sampCols=NULL,multiPred=FALSE){
	pdf(file=sprintf("%s_fits.pdf",outName),width=14,height=12)
		par(mfrow=getPlotDims(length(outs)))
		for(i in 1:length(outs)){
			betaPPS(outs[[i]]$db,outs[[i]]$fit,500,predName=predNames[i],multiPred=multiPred,sampCols=sampCols)
		}
	dev.off()
	#
	pdf(file=sprintf("%s_phyloFit.pdf",outName),width=12,height=12)
		for(i in 1:length(outs)){
			modAdViz(outs[[i]]$db,outs[[i]]$fit,predName=predNames[i],
						nPPS=500,tree=sampPhy,xlim=c(0,0.05),
						valRange=NULL,qnt=0.999,adj=0.5)
		}
	dev.off()
	#
	pdf(file=sprintf("%s_phyloFit_log.pdf",outName),width=12,height=12)
		for(i in 1:length(outs)){
			modAdViz(outs[[i]]$db,outs[[i]]$fit,predName=predNames[i],
						nPPS=500,tree=sampPhy,xlim=c(0,0.05),
						valRange=NULL,qnt=0.999,adj=0.5,logX=TRUE)
		}
	dev.off()
	#
	pdf(file=sprintf("%s_betas.pdf",outName),width=14,height=10)
		postBetaPlot(outs,predNames,reorder=TRUE,cols=NULL,stdize=TRUE,multiPred=multiPred,qnt=0.999)
	dev.off()
}

getMeanLine <- function(beta,gamma,X,nPPS){
	lns <- lapply(1:nPPS,function(i){beta[i]*X+gamma[i]})
	mnLn <- colMeans(Reduce("rbind",lns))
	return(mnLn)
}

getBestLine <- function(lpp,beta,gamma,X){
	bestIter <- which.max(lpp)
	bestLine <- beta[bestIter] * X + gamma[bestIter]
	return(bestLine)
}

invLogit <- function(x){
	return(exp(x)/(1+exp(x)))
}

logit <- function(x){
	return(log(x/(1-x)))
}

betaPPS <- function(db,fit,nPPS,predName,multiPred=FALSE,sampCols=NULL){
	if(!is.null(sampCols)){
		sampCols <- 	getSampCols(db,sampCols)
	} else {
		sampCols <- "red"
	}
	lpp <- get_logposterior(fit,inc_warmup=FALSE)[[1]]
	nIter <- length(lpp)
	sampledIter <- sample(1:nIter,nPPS,replace=TRUE)
	if(multiPred){
		b <- "beta[1]"
	} else {
		b <- "beta"
	}
	beta <- extract(fit,b,permute=FALSE,inc_warmup=FALSE)[,1,1]
	sig <- 1-2*(abs(0.5-ecdf(beta)(0)))
	siglwd <- ifelse(sig<0.05,4,0.5)
	beta <- beta[sampledIter]
	gamma <- extract(fit,"gamma",permute=FALSE,inc_warmup=FALSE)[sampledIter,1,1]
	# alpha <- extract(fit,"alpha",permute=FALSE,inc_warmup=FALSE)[sampledIter,1,1]
	# phi <- extract(fit,"phi",permute=FALSE,inc_warmup=FALSE)[sampledIter,1,1]
	# mvnMean <- lapply(1:nPPS,function(i){c(gamma[i] + beta[i]*db$X)})
	# theta <- lapply(1:nPPS,function(i){MASS::mvrnorm(1,mu=mvnMean[[i]],Sigma=alpha[i]*db$relMat)})
	# mu <- lapply(theta,function(x){invLogit(x)})
	# shape1 <- lapply(1:nPPS,function(i){phi[[i]]*mu[[i]]})
	# shape2 <- lapply(1:nPPS,function(i){phi[[i]]*(1-mu[[i]])})
	shape1 <- extract(fit,"shape1",permute=FALSE)[sampledIter,1,]
	shape2 <- extract(fit,"shape2",permute=FALSE)[sampledIter,1,]
	pps <- lapply(1:nPPS,function(i){rbeta(db$N,shape1[i,],shape2[i,])})
	ppsCIs <- getPPSci(db,pps)
	if(multiPred){
		x <- db$X[1,]
	} else {
		x <- db$X
	}
	plot(x,db$Y,
		ylim=range(c(db$Y,unlist(ppsCIs))),
		xlab="predictor",ylab="response",
		type='n',main=sprintf("%s (p=%s)",predName,round(sig,4)))
		invisible(
			lapply(1:db$N,
				function(n){
					segments(x0=x[n],
							 y0=ppsCIs[[n]][1],
							 x1=x[n],
							 y1=ppsCIs[[n]][2],
							 lwd=0.75)
				})
		)
	points(x,db$Y,col=sampCols,pch=19,cex=1)
	lnx <- seq(min(x),max(x),length.out=100)
	mnLn <- getMeanLine(beta=beta,gamma=gamma,X=lnx,nPPS=nPPS)
	lines(lnx,invLogit(mnLn),lwd=siglwd)
	box(lwd=2,col=ifelse(sig<0.05,"red","black"))
}

getPPSci <- function(db,pps){
	ppsByN <- lapply(1:db$N,
				function(n){
					quantile(unlist(lapply(pps,"[[",n)),c(0.025,0.975))
				})
	return(ppsByN)
}

colFunc <- function (x, cols, nCols, valRange){
    if (is.null(valRange)) {
        valRange <- c(min(x), max(x))
    }
    cols <- (grDevices::colorRampPalette(cols))(nCols)[findInterval(x, 
        seq(valRange[1], valRange[2], length.out = nCols))]
    return(cols)
}

checkSig <- function(beta){
	sig <- 1-2*(abs(0.5-ecdf(beta)(0)))
	sig <- ifelse(sig < 0.05,TRUE,FALSE)
	return(sig)
}

postBetaPlot <- function(outs,predNames,reorder=TRUE,cols=NULL,stdize=FALSE,multiPred=FALSE,qnt=1,...){
#	recover()
	if(multiPred){
		b <- "beta[1]"
	} else {
		b <- "beta"
	}
	betas <- lapply(1:length(outs),
				function(i){
					extract(outs[[i]]$fit,b,inc_warmup=FALSE,permute=FALSE)[,1,1]
				})
	if(stdize){
		betas <- lapply(1:length(betas),
						function(i){
							if(length(unique(outs[[i]]$db$X[1,])) > 3){
								betas[[i]] * sd(outs[[i]]$db$X[1,])
							} else {
								betas[[i]]
							}
						})
	}
	meanBetas <- unlist(lapply(betas,mean))
	nPredictors <- length(betas)
	sig <- unlist(lapply(betas,function(beta){checkSig(beta)}))
	if(any(sig)){
		predNames[which(sig)] <- paste0(predNames[which(sig)]," * ")
	}
	if(is.null(cols)){
		cols <- colFunc(meanBetas,cols=c("blue","red"),nPredictors,valRange=range(meanBetas))	#rep("blue",nPredictors)
	}
	betaDens <- lapply(betas,density)
	if(reorder){
		predOrder <- order(meanBetas)
	} else {
		predOrder <- 1:nPredictors
	}
	plot(0,type='n',xlim=range(betas)+c(-diff(range(betas))/4,0),yaxt='n',xlab="effect size",ylab="",bty='n',ylim=c(0,(nPredictors+1)*1.25))
		text(x=min(unlist(betas))-diff(range(betas))/8,y=(0.4 + 1:nPredictors)*1.25,labels=predNames[predOrder],srt=0,font=2,cex=1.3)
		abline(v=0,lty=2,lwd=0.5,col=adjustcolor(1,0.5))
	invisible(
		lapply(1:nPredictors,
			function(i){
				plotDens(ymin=i+i*0.25,x=betas[[predOrder[i]]],d=betaDens[[predOrder[i]]],col=cols[predOrder[i]],qnt=qnt)
			}))
}

parallelizing <- function(args) {
    if (args[["parallel"]]) {
        if (!foreach::getDoParRegistered()) {
            if (is.null(args[["nNodes"]])) {
                nNodes <- parallel::detectCores() - 1
            }
            else {
                nNodes <- args[["nNodes"]]
            }
            cl <- parallel::makeCluster(nNodes)
            doParallel::registerDoParallel(cl)
            message("\nRegistered doParallel with ", nNodes, 
                " workers\n")
        }
        else {
            message("\nUsing ", foreach::getDoParName(), " with ", 
                foreach::getDoParWorkers(), " workers")
        }
        d <- foreach::`%dopar%`
    }
    else {
        message("\nRunning sequentially with a single worker\n")
        d <- foreach::`%do%`
    }
    return(d)
}

makeDB <- function(X,Y,phyStr){
	md <- which(is.na(X))
	if(any(md)){
		phyStr <- phyStr[-md,-md]
		Y <- Y[-md]
		X <- X[-md]
	}
	db <- list("N" = nrow(phyStr),
			   "Y" = Y,
			   "nX"= 1,
			   "X" = matrix(X,nrow=1),
			   "relMat" = phyStr)
	return(db)
}

makeMultiDB <- function(X,Y,phyStr){
	md <- unique(unlist(lapply(X,function(x){which(is.na(x))})))
	if(any(md)){
		phyStr <- phyStr[-md,-md]
		Y <- Y[-md]
		X <- lapply(X,function(x){x[-md]})
	}
	db <- list("N" = nrow(phyStr),
			   "Y" = Y,
			   "nX"= length(X),
			   "X" = Reduce("rbind",X,init=NULL),
			   "relMat" = phyStr)
	return(db)
}

getTipOrder <- function(tree){
	isTip <- tree$edge[,2] <= length(tree$tip.label)
	tipOrder <- tree$edge[isTip,2]
	return(tree$tip.label[tipOrder])
}

getDens <- function(x,adj=2,logX=FALSE){
#	recover()
	# if(length(unique(x)) > 1){
		# toDrop <- which(x < quantile(x,qmn))
		# toDrop <- unique(c(toDrop,which(x > quantile(x,qmx))))
		# if(length(toDrop) > 0){
			# x <- x[-toDrop]
		# }
	# }
	if(logX){
		x <- log(x)
	}
	d <- density(x,adjust=adj)
	return(d)
}

plotDens <- function(ymin=0,x,d,col,alpha=0.3,xmin=NULL,peakheight=1,qnt=1,logX=FALSE){
	if(logX){
		x <- log(x)
	}
	if(qnt != 1){
		qmn <- (1-qnt)/2
		qmx <- 1-qmn
		qx <- quantile(x,c(qmn,qmx))
		toDrop <- which(d$x < qx[1])
		toDrop <- c(toDrop,which(d$x > qx[2]))
		dx <- d$x[-toDrop]
		dy <- d$y[-toDrop]
	} else {
		dx <- d$x
		dy <- d$y
	}
	if(is.null(xmin)){
		xmin <- min(dx)
	}
	yvec <- dy/(peakheight*max(dy))
	polygon(x=c(xmin,xmin,dx,max(dx),xmin),
				y=c(ymin,ymin,ymin+yvec,ymin,ymin),
				col=adjustcolor(col,alpha))
}

modAdViz <- function(db,fit,predName,nPPS,tree,xlim,valRange=NULL,qnt=0.95,adj=2,logX=FALSE){
	lpp <- get_logposterior(fit,inc_warmup=FALSE)[[1]]
	nIter <- length(lpp)
	sampledIter <- sample(1:nIter,nPPS,replace=TRUE)
	shape1 <- extract(fit,"shape1",permute=FALSE)[sampledIter,1,]
	shape2 <- extract(fit,"shape2",permute=FALSE)[sampledIter,1,]
	pps <- lapply(1:nPPS,function(i){rbeta(db$N,shape1[i,],shape2[i,])})
	pps <- Reduce("rbind",pps,init=NULL)
	ppsDens <- lapply(1:db$N,function(n){getDens(pps[,n],adj=adj,logX=logX)})
	spNames <- row.names(db$relMat)
	tree <- ape::keep.tip(tree,gsub("_"," ",spNames))
	tree <- ape::ladderize(tree,right=FALSE)
	tipOrder <- getTipOrder(tree)
	nSp <- length(tipOrder)
	spOrder <- match(tipOrder,spNames)
	# ppsFits <- unlist(
				# lapply(1:db$N,
					# function(n){
						# #abs(mean(pps[,n])-db$Y[n])
						# 2*(abs(0.5-ecdf(pps[,n])(db$Y[n])))
				# }))
	virCols <- viridis::viridis_pal(alpha=1,begin=0,end=1,direction=1,option="D")(100)
	cols <- colFunc(log(db$Y),cols=virCols,nCols=100,valRange=valRange) #ppsFits
	par(mfrow=c(1,2),oma=c(3,0,0,0))
	ape::plot.phylo(tree,label.offset=50,cex=0.7,no.margin=TRUE,x.lim=c(0,1950),tip.color=cols,edge.width=1.5)
	text(x=450,y=70,labels=sprintf("predictor:\n%s",predName))
	if(logX){
		xlim <- range(log(db$Y)) + c(-4,1)
		y <- log(db$Y)
	} else {
		y <- db$Y
	}
	plot(0,type='n',
			xlim=xlim,
			yaxt='n',xaxt='n',xlab="",ylab="",bty='n',
			ylim=c(1,nSp)) #c(0,(length(cnDens)+1)*1.25))
	axis(side=1,
		at=c(xlim[1],(xlim[2]+xlim[1])/2,xlim[2]),
		labels=round(c(xlim[1],(xlim[2]+xlim[1])/2,xlim[2]),2))
	#abline(h=0:nSp,lty=2)
	#text(x=xtext,y=(0.4 + 1:nSp)*1.25,labels=spNames[spOrder],srt=0)
	invisible(
		lapply(nSp:1,
			function(i){
				plotDens(ymin=i-0.5,
						 x=pps[,spOrder[i]],
						 d=ppsDens[[spOrder[i]]],
						 col=cols[spOrder[i]],
						 alpha=0.9,
						 qnt=qnt,
						 xmin=NULL,
						 peakheight=0.5,
						 logX=logX)
				#text(x=0.95,y=i,labels=spNames[spOrder[i]])
				points(x=y[spOrder[i]],i-0.25,pch=17,col="red")
				abline(h=spOrder[i]-0.5,lty=3,col="gray")	
			}))
}

phyloViz_scatter <- function(loo,predName,valRange=NULL){
	if(is.null(valRange)){
		valRange <- c(0,1)
	}
	y <- lapply(1:loo$db$N,function(i){
					rstan::extract(loo$looReps[[i]]$fit,"y",
									inc_warmup=FALSE,permute=FALSE)[,1,loo$db$N]
				})
	yMean <- unlist(lapply(y,mean))
	fits <- unlist(lapply(1:loo$db$N,
							function(n){
								2*(abs(0.5-ecdf(y[[n]])(loo$looReps[[n]]$looDB$Ystar)))}))
	virCols <- viridis::viridis_pal(alpha=1,begin=0,end=1,direction=1,option="D")(100)
	cols <- colFunc(fits,cols=virCols,nCols=100,valRange=valRange)
	plot(loo$db$Y,yMean,
			xlim=range(c(loo$db$Y,yMean)),
			ylim=range(c(loo$db$Y,yMean)),
			xlab="observed diversity",
			ylab="leave-one-out predicted diversity",type='n',
			main=predName)
		abline(a=0,b=1,col="red",lty=2)
		invisible(
			lapply(1:loo$db$N,
				function(n){
					qnt <- quantile(y[[n]],c(0.025,0.975))
					segments(x0=loo$db$Y[n],x1=loo$db$Y[n],
								y0=qnt[1],y1=qnt[2],lwd=1.5,col=adjustcolor(cols[n],0.5))
			})
		)
	points(loo$db$Y,yMean,pch=18,cex=1.5,col=cols)
}



