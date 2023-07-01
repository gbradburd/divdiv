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

betaPPS <- function(db,fit,nPPS,predName){
	#recover()
	lpp <- get_logposterior(fit,inc_warmup=FALSE)[[1]]
	nIter <- length(lpp)
	sampledIter <- sample(1:nIter,nPPS,replace=TRUE)
	beta <- extract(fit,"beta",permute=FALSE,inc_warmup=FALSE)[,1,1]
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
	pps <- lapply(1:nPPS,function(i){rbeta(db$N,shape1[[i]],shape2[[i]])})
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
	lines(x,invLogit(mnLn),lwd=siglwd)
	box(lwd=2,col=ifelse(sig<0.05,"red","black"))
}

getPPSci <- function(db,pps){
	ppsByN <- lapply(1:db$N,
				function(n){
					quantile(unlist(lapply(pps,"[[",n)),c(0.025,0.975))
				})
	return(ppsByN)
}

expPPS <- function(db,fit,nPPS,predName){
	nIter <- length(get_logposterior(fit)[[1]])
	sampledIter <- sample(1:(nIter/2),nPPS,replace=TRUE)
	beta <- extract(fit,"beta",permute=FALSE)[,1,1]
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
	sig <- 1-2*(abs(0.5-ecdf(beta)(0)))
	siglwd <- ifelse(sig<0.05,4,0.5)
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
	lines(x,1/exp(mnLn),lwd=siglwd)
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

postBetaPlot <- function(out,predNames,reorder=TRUE,cols=NULL,...){
	betas <- lapply(out$fits,
				function(fit){
					extract(fit,"beta",inc_warmup=FALSE,permute=FALSE)[,1,1]
				})
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
		abline(v=0,lty=2,lwd=0.5)
	invisible(
		lapply(1:nPredictors,
			function(i){
				plotDens(ymin=i+i*0.25,betaDens[[predOrder[i]]],col=cols[predOrder[i]])
			}))
}

makeLooDB <- function(toDrop,db){
	db$N <- db$N-1
	db$Y <- db$Y[-toDrop]
	db$X <- db$X[1,-toDrop,drop=FALSE]
	db$relMat <- db$relMat[-toDrop,-toDrop]
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

getTipOrder <- function(tree){
	isTip <- tree$edge[,2] <= length(tree$tip.label)
	tipOrder <- tree$edge[isTip,2]
	return(tree$tip.label[tipOrder])
}

plotDens <- function(ymin=0,d,col,alpha=0.3,xmin=NULL,peakheight=1){
	if(is.null(xmin)){
		xmin <- min(d$x)
	}
	yvec <- d$y/(peakheight*max(d$y))
	polygon(x=c(xmin,d$x,xmin),
				y=c(ymin,ymin+yvec,ymin),
				col=adjustcolor(col,alpha))
}

getLooDens <- function(x){
	if(length(unique(x)) > 1){
		toDrop <- which(x < quantile(x,0.1))
		toDrop <- unique(c(toDrop,which(x > quantile(x,0.9))))
		x <- x[-toDrop]
	}
	d <- density(x,adjust=2)
	return(d)
}

phylooViz <- function(db,CNsamples,tree,xlim=c(0.95,1.01),valRange=NULL){
	#recover()
	if(is.null(valRange)){
		valRange <- c(0,1)
	}
	spNames <- row.names(db$relMat)
	tree <- ape::keep.tip(tree,gsub("_"," ",spNames))
	tree <- ape::ladderize(tree,right=FALSE)
	tipOrder <- getTipOrder(tree)
	nSp <- length(tipOrder)
	spOrder <- match(tipOrder,spNames)
	CNsamples <- looCNsamples
	cnFits <- unlist(
				lapply(1:db$N,
					function(n){
						abs(mean(CNsamples[[n]])-db$Y[n])
						#2*(abs(0.5-ecdf(CNsamples[[n]])(db$Y[n])))
				}))
	cnDens <- lapply(CNsamples,getLooDens)
	virCols <- viridis::viridis_pal(alpha=1,begin=0,end=1,direction=1,option="D")(100)
	cols <- colFunc(cnFits,cols=virCols,nCols=100,valRange=valRange)
	par(mfrow=c(1,2),oma=c(3,0,0,0))
	ape::plot.phylo(tree,label.offset=50,cex=0.7,no.margin=TRUE,x.lim=c(0,1950),tip.color=cols,edge.width=1.5)
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
						 d=cnDens[[spOrder[i]]],
						 col=cols[spOrder[i]],
						 alpha=0.9,
						 xmin=NULL,
						 peakheight=0.5)
				#text(x=0.95,y=i,labels=spNames[spOrder[i]])
				points(x=db$Y[spOrder[i]],i-0.25,pch=17,col="red")
				abline(h=spOrder[i]-0.5,lty=3,col="gray")	
			}))
}








