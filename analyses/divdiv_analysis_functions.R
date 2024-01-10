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

mvnPPS <- function(db,fit,nPPS,predName){
	nIter <- length(get_logposterior(fit)[[1]])
	sampledIter <- sample(1:(nIter/2),nPPS,replace=TRUE)
	beta <- extract(fit,"beta",permute=FALSE)[,1,1]
	sig <- 1-2*(abs(0.5-ecdf(beta)(0)))
	siglwd <- ifelse(sig<0.05,4,0.5)
	beta <- beta[sampledIter]
	alpha <- extract(fit,"alpha",permute=FALSE,inc_warmup=FALSE)[sampledIter,1,1]
	gamma <- extract(fit,"gamma",permute=FALSE,inc_warmup=FALSE)[sampledIter,1,1]
	mu <- lapply(1:nPPS,function(i){c(gamma[i] + beta[i]*db$X)})
	pps <- lapply(1:nPPS,function(i){MASS::mvrnorm(1,mu=mu[[i]],Sigma=alpha[i]*db$relMat)})
	ppsCIs <- getPPSci(db,pps)
	plot(db$X,db$Y,
		ylim=range(c(db$Y,unlist(ppsCIs))),
		xlab="predictor",ylab="response",
		type='n',main=sprintf("%s (p=%s)",predName,round(sig,4)))
		invisible(
			lapply(1:db$N,
				function(n){
					segments(x0=db$X[n],
							 y0=ppsCIs[[n]][1],
							 x1=db$X[n],
							 y1=ppsCIs[[n]][2],
							 lwd=0.75)
				})
		)
	points(db$X,db$Y,col="red",pch=19,cex=1)
	x <- seq(min(db$X),max(db$X),length.out=100)
	mnLn <- getMeanLine(beta=beta,gamma=gamma,X=x,nPPS=nPPS)
	lines(x,mnLn,lwd=siglwd)
	box(lwd=2,col=ifelse(sig<0.05,"red","black"))
}

invLogit <- function(x){
	return(exp(x)/(1+exp(x)))
}

logit <- function(x){
	return(log(x/(1-x)))
}

betaPPS <- function(db,fit,nPPS,predName,multiPred=FALSE){
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
	points(x,db$Y,col="red",pch=19,cex=1)
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

expPPS <- function(db,fit,nPPS,predName,multiPred=FALSE){
	lpp <- get_logposterior(fit,inc_warmup=FALSE)[[1]]
	nIter <- length(lpp)
	sampledIter <- sample(1:(nIter/2),nPPS,replace=TRUE)
	if(multiPred){
		b <- "beta[1]"
	} else {
		b <- "beta"
	}
	beta <- extract(fit,b,permute=FALSE,inc_warmup=FALSE)[,1,1]
	sig <- 1-2*(abs(0.5-ecdf(beta)(0)))
	siglwd <- ifelse(sig<0.05,4,0.5)
	beta <- beta[sampledIter]
	#alpha <- extract(fit,"alpha",permute=FALSE,inc_warmup=FALSE)[sampledIter,1,1]
	gamma <- extract(fit,"gamma",permute=FALSE,inc_warmup=FALSE)[sampledIter,1,1]
	#mvnMean <- lapply(1:nPPS,function(i){c(gamma[i] + beta[i]*db$X)})
	#theta <- lapply(1:nPPS,function(i){MASS::mvrnorm(1,mu=mvnMean[[i]],Sigma=alpha[i]*db$relMat)})
	#lambda <- lapply(theta,function(x){exp(x)})
	lambda <- extract(fit,"lambda",permute=FALSE)[sampledIter,1,]
	pps <- lapply(1:nPPS,function(i){rexp(db$N,rate=lambda[[i]])})
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
	points(x,db$Y,col="red",pch=19,cex=1)
	lnx <- seq(min(x),max(x),length.out=100)
	mnLn <- getMeanLine(beta=beta,gamma=gamma,X=lnx,nPPS=nPPS)
	lines(lnx,1/exp(mnLn),lwd=siglwd)
#	lines(0:max(db$X),1/exp(mean(beta)*(0:max(db$X))+mean(gamma)),lwd=siglwd)
	box(lwd=2,col=ifelse(sig<0.05,"red","black"))
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

postBetaPlot <- function(out,predNames,reorder=TRUE,cols=NULL,stdize=FALSE,multiPred=FALSE,qnt=1,...){
#	recover()
	if(multiPred){
		b <- "beta[1]"
	} else {
		b <- "beta"
	}
	betas <- lapply(out$fits,
				function(fit){
					extract(fit,b,inc_warmup=FALSE,permute=FALSE)[,1,1]
				})
	if(stdize){
		betas <- lapply(1:length(betas),
						function(i){
							if(length(unique(out$db[[i]]$X[1,])) > 3){
								betas[[i]] * sd(out$db[[i]]$X[1,])
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

makeLooDB <- function(toDrop,db){
	db$Nstar <- db$N
	db$N <- db$N-1
	db$Ystar <- db$Y[toDrop]
	db$Y <- db$Y[-toDrop]
	reord <- c(c(1:db$Nstar)[-toDrop],toDrop)
	db$relMat <- db$relMat[reord,reord]
	if(db$nX == 1){
		db$X <- db$X[reord]
	} else {
		xMeans <- rowMeans(db$X[2:db$nX,-toDrop])
		db$X[2:db$nX,toDrop] <- xMeans
		db$X <- db$X[,reord]
	}
	return(db)
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

looRep <- function(toDrop,sampleName,db,nIter,mod){
	looDB <- makeLooDB(toDrop,db)
	fit <- rstan::sampling(object=mod,
							data=looDB,
							iter=nIter,
							thin=nIter/500,
							chains=1,
							control=setNames(list(15),"max_treedepth"))
	loo <- list("fit" = fit,
				"looDB" = looDB,
				"dropped" = sampleName)
	return(loo)
}

phyloo <- function(db,mod,sampleNames,nIter,parallel=TRUE,nNodes=4,prefix){
	prespecified <- conStruct:::parallel.prespecify.check(args <- as.list(environment()))
    `%d%` <- parallelizing(args <- as.list(environment()))
	looReps <- foreach::foreach(i = 1:db$N,.export=c("looRep","makeLooDB"),.packages="rstan") %d% {
				looRep(toDrop=i,sampleName=sampleNames[i],db=db,nIter=nIter,mod=mod)
			}
	looReps <- setNames(looReps,paste0("looRep",1:db$N))
	loo <- list("db"=db,
				"looReps"=looReps)
	conStruct:::end.parallelization(prespecified)
	save(loo,file=paste0(prefix,"_loo.Robj"))
	return("loo done")
}

getCNmean <- function(masterDB,i,beta,gamma,sig12,sig22inv,link){
	mu1 <- beta*masterDB$X[,i,drop=FALSE] + gamma
	mu2 <- beta*masterDB$X[,-i,drop=FALSE] + gamma
	cnMean <- mu1 + sig12 %*% sig22inv %*% t(link(masterDB$Y[-i])-mu2)
	return(cnMean)
}

getCNvar <- function(sig11,sig12,sig22inv,sig21){
	cnVar <- sig11 - sig12 %*% sig22inv %*% sig21
	return(cnVar)
}

looCN <- function(nCNsamples,sampleNames,masterDB,looRep,link,invLink){
	i <- which(sampleNames==looRep$dropped)
	alpha <- mean(extract(looRep$fit,"alpha",inc_warmup=FALSE,permute=FALSE))
	beta <- mean(extract(looRep$fit,"beta",inc_warmup=FALSE,permute=FALSE))
	gamma <- mean(extract(looRep$fit,"gamma",inc_warmup=FALSE,permute=FALSE))
	sig22inv <- MASS::ginv(alpha*masterDB$relMat[-i,-i])
	sig11 <- alpha*masterDB$relMat[i,i]
	sig12 <- alpha*masterDB$relMat[i,-i,drop=FALSE]
	sig21 <- alpha*masterDB$relMat[-i,i,drop=FALSE]
	cnMean <- getCNmean(masterDB,i,beta,gamma,sig12,sig22inv,link)
	cnVar <- getCNvar(sig11,sig12,sig22inv,sig21)
	cnSamples <- invLink(rnorm(nCNsamples,mean=cnMean,sd=sqrt(cnVar)))
	return(cnSamples)
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

phylooViz <- function(loo,predName,tree,xlim,valRange=NULL,qnt=0.95,adj=2,logX=FALSE){
	if(is.null(valRange)){
		valRange <- c(0,1)
	}
	y <- lapply(1:loo$db$N,function(i){
					rstan::extract(loo$looReps[[i]]$fit,"y",
									inc_warmup=FALSE,permute=FALSE)[,1,loo$db$N]
				})
	yDens <- lapply(y,getDens,adj=adj,logX=logX)
	spNames <- row.names(loo$db$relMat)
	tree <- ape::keep.tip(tree,gsub("_"," ",spNames))
	tree <- ape::ladderize(tree,right=FALSE)
	tipOrder <- getTipOrder(tree)
	nSp <- length(tipOrder)
	spOrder <- match(tipOrder,spNames)
	fits <- unlist(lapply(1:loo$db$N,
							function(n){
								2*(abs(0.5-ecdf(y[[n]])(loo$looReps[[n]]$looDB$Ystar)))}))
	virCols <- viridis::viridis_pal(alpha=1,begin=0,end=1,direction=1,option="D")(100)
	cols <- colFunc(fits,cols=virCols,nCols=100,valRange=valRange)
	par(mfrow=c(1,2),oma=c(3,0,0,0))
	ape::plot.phylo(tree,label.offset=50,cex=0.7,no.margin=TRUE,x.lim=c(0,1950),tip.color=cols,edge.width=1.5)
	text(x=450,y=70,labels=sprintf("predictor:\n%s",predName))
	if(logX){
		xlim <- range(log(loo$db$Y)) + c(-4,1)
		yax <- log(loo$db$Y)
	} else {
		yax <- loo$db$Y
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
						 x=y[[spOrder[i]]],
						 d=yDens[[spOrder[i]]],
						 col=cols[spOrder[i]],
						 alpha=0.9,
						 qnt=qnt,
						 xmin=NULL,
						 peakheight=0.5,
						 logX=logX)
				#text(x=0.95,y=i,labels=spNames[spOrder[i]])
				points(x=yax[spOrder[i]],i-0.25,pch=17,col="red")
				abline(h=spOrder[i]-0.5,lty=3,col="gray")	
			}))
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



