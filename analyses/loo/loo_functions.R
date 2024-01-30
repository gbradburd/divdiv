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