divdivAnalyses <- function(z,X,y,phyStr,mod,outName,nIter=1e4,parallel=TRUE,nNodes=4,filterKeep=NULL,runNames=NULL){
    `%d%` <- parallelizing(args <- as.list(environment()))
    export <- c("divdivAnalysis","makeDB","makeMultiDB")
    outs <- foreach::foreach(i = 1:length(X),
    						.export=export,
    						.packages="rstan") %d% {
								divdivAnalysis(z=z,x=X[[i]],y=y,phyStr=phyStr,
												mod=mod,nIter=nIter,filterKeep=filterKeep)
			}
	doParallel::stopImplicitCluster()
	if(is.null(runNames)){
		runNames <- paste0("pred",1:length(X))
	}
	names(outs) <- runNames
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
							control=setNames(list(15),"max_treedepth"),
							seed=123)
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
			pps <- getBetaPPS(db=outs[[i]]$db,fit=outs[[i]]$fit,nPPS=500,multiPred= multiPred)
			plotBetaPPS(db=outs[[i]]$db,ppsOut=pps,sampCols= sampCols,pt.cex=1.5,trp=0.8,main=predNames[i],xlab="predictor",ylab="diversity")
		}
	dev.off()
	#
	# pdf(file=sprintf("%s_phyloFit.pdf",outName),width=12,height=12)
		# for(i in 1:length(outs)){
			# modAdViz(outs[[i]]$db,outs[[i]]$fit,predName=predNames[i],
						# nPPS=500,tree=sampPhy,xlim=c(0,0.05),
						# valRange=NULL,qnt=0.999,adj=0.5)
		# }
	# dev.off()
	# #
	# pdf(file=sprintf("%s_phyloFit_log.pdf",outName),width=12,height=12)
		# for(i in 1:length(outs)){
			# modAdViz(outs[[i]]$db,outs[[i]]$fit,predName=predNames[i],
						# nPPS=500,tree=sampPhy,xlim=c(0,0.05),
						# valRange=NULL,qnt=0.999,adj=0.5,logX=TRUE)
		# }
	# dev.off()
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

getBetaPPS <- function(db,fit,nPPS,multiPred=FALSE){
	lpp <- get_logposterior(fit,inc_warmup=FALSE)[[1]]
	nIter <- length(lpp)
	sampledIter <- sample(1:nIter,nPPS,replace=TRUE)
	if(multiPred){
		b <- "beta[1]"
	} else {
		b <- "beta"
	}
	beta <- extract(fit,b,permute=FALSE,inc_warmup=FALSE)[,1,1]
	beta <- beta[sampledIter]
	gamma <- extract(fit,"gamma",permute=FALSE,inc_warmup=FALSE)[sampledIter,1,1]
	shape1 <- extract(fit,"shape1",permute=FALSE)[sampledIter,1,]
	shape2 <- extract(fit,"shape2",permute=FALSE)[sampledIter,1,]
	pps <- lapply(1:nPPS,function(i){rbeta(db$N,shape1[i,],shape2[i,])})
	ppsCIs <- getPPSci(db,pps)
	if(multiPred){
		x <- db$X[1,]
	} else {
		x <- db$X
	}
	lnx <- seq(min(x),max(x),length.out=100)
	mnLn <- getMeanLine(beta=beta,gamma=gamma,X=lnx,nPPS=nPPS)
	ppsOut <- list("meanLine"=list("lnx"=lnx,"mnLn"=mnLn),
					"pps"=pps,
					"x"=x,
					"ppsCIs"=ppsCIs)
	return(ppsOut)
}

plotBetaPPS <- function(db,ppsOut,sampCols=NULL,pt.cex=1.5,trp=1,...){
	if(!is.null(sampCols)){
		sampCols <- getSampCols(db,sampCols)
	} else {
		sampCols <- "red"
	}
	plot(ppsOut$x,db$Y,
		ylim=range(c(db$Y,unlist(ppsOut$ppsCIs))),
		type='n',...)
		invisible(
			lapply(1:db$N,
				function(n){
					segments(x0=ppsOut$x[n],
							 y0=ppsOut$ppsCIs[[n]][1],
							 x1=ppsOut$x[n],
							 y1=ppsOut$ppsCIs[[n]][2],
							 lwd=0.75,
							 col=adjustcolor("gray",0.7))
				}))
	points(ppsOut$x,db$Y,col=adjustcolor(sampCols,trp),pch=19,cex=pt.cex)
	lines(ppsOut$meanLine$lnx,invLogit(ppsOut$meanLine$mnLn),lwd=1.5)
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
	#recover()
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
		cols <- colFunc(meanBetas,cols=viridis::viridis(n=nPredictors),nPredictors,valRange=range(meanBetas))
	}
	betaDens <- lapply(betas,density)
	if(reorder){
		predOrder <- order(meanBetas)
	} else {
		predOrder <- 1:nPredictors
	}
	yTxtCoords <- 0.6+c(1:nPredictors)+c(1:nPredictors)*0.05
	yDnsCoords <- c(1:nPredictors)+c(1:nPredictors)*0.05
	plot(0,type='n',xlim=range(betas)+c(-diff(range(betas))/4,0),yaxt='n',xlab="Standardized effect size",ylab="",bty='n',ylim=c(0.9,(nPredictors+1)*1.10))
		text(x=min(unlist(betas))-diff(range(betas))/4,y=yTxtCoords,labels=predNames[predOrder],srt=0,font=2,cex=1.3,pos=4)
		abline(v=0,lty=2,lwd=0.5,col=adjustcolor(1,0.5))
	invisible(
		lapply(nPredictors:1,
			function(i){
				plotDens(ymin=yDnsCoords[i],x=betas[[predOrder[i]]],d=betaDens[[predOrder[i]]],col=cols[predOrder[i]],qnt=qnt,peakheight=0.5,alpha=1,border="black")
			}))
}

postBetaPlot_allPreds <- function(outs,predNames,reorder=TRUE,cols=NULL,stdize=FALSE,qnt=1,...){
	recover()
	betas <- lapply(1:length(predNames),
				function(i){
					extract(outs$fit,sprintf("beta[%s]",i),inc_warmup=FALSE,permute=FALSE)[,1,1]
				})
	if(stdize){
		betas <- lapply(1:length(betas),
						function(i){
							if(length(unique(outs$db$X[1,])) > 3){
								betas[[i]] * sd(outs$db$X[1,])
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
		cols <- colFunc(meanBetas,cols=viridis::viridis(n=nPredictors),nPredictors,valRange=range(meanBetas))
	}
	betaDens <- lapply(betas,density)
	if(reorder){
		predOrder <- order(meanBetas)
	} else {
		predOrder <- 1:nPredictors
	}
	yTxtCoords <- 0.6+c(1:nPredictors)+c(1:nPredictors)*0.05
	yDnsCoords <- c(1:nPredictors)+c(1:nPredictors)*0.05
	plot(0,type='n',xlim=range(betas)+c(-diff(range(betas))/4,0),yaxt='n',xlab="Standardized effect size",ylab="",bty='n',ylim=c(0.9,(nPredictors+1)*1.10))
		text(x=min(unlist(betas))-diff(range(betas))/4,y=yTxtCoords,labels=predNames[predOrder],srt=0,font=2,cex=1.3,pos=4)
		abline(v=0,lty=2,lwd=0.5,col=adjustcolor(1,0.5))
	invisible(
		lapply(nPredictors:1,
			function(i){
				plotDens(ymin=yDnsCoords[i],x=betas[[predOrder[i]]],d=betaDens[[predOrder[i]]],col=cols[predOrder[i]],qnt=qnt,peakheight=0.5,alpha=1,border="black")
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

plotDens <- function(ymin=0,x,d,col,alpha=0.3,xmin=NULL,peakheight=1,qnt=1,logX=FALSE,border="black"){
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
				col=adjustcolor(col,alpha),border=border)
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

highlightClade <- function(tree,dSp1,dSp2,col,x0=NULL,x1=NULL,y0=NULL,y1=NULL,rounding=0.03){
	#recover()
	m <- ape::getMRCA(tree,tip=c(dSp1,dSp2))
	offset <- max(nodeHeights(sampPhy))/100
	decs <- getDescendants(tree,m)
	if(is.null(x0)){
		x0 <- phytools::nodeheight(tree,m) - offset
	}
	if(is.null(x1)){
		x1 <- max(unlist(lapply(1:length(decs),function(i){nodeheight(tree,decs[i])})))
	}
	if(is.null(y0)){
		y0 <- min(node.height(tree)[decs]) - 0.5
	}
	if(is.null(y1)){
		y1 <- max(node.height(tree)[decs]) + 0.5
	}
	berryFunctions::roundedRect(xleft=x0,ybottom=y0,
								xright=x1,ytop=y1,rounding=rounding,
								col=col,border=NA,aspcorrect=TRUE)
}

addPhylopics <- function(cladeCols,ysize=4){
#	recover()
	xleft <- 90
	ytop <- 80
	dshift <- 5
	legendshift <- 100
	txt.cex=1.4
	pos=4
	mammal <- rphylopic::get_phylopic(uuid="44ee07ec-f829-49f9-b242-f4b92bb9cf73") #Halichoerus grypus
	vasc_plant <- rphylopic::get_phylopic(uuid="bacc6a59-80d2-461c-92d4-83cabc91cc39") #Rhizophora mangle
	ochrophyta <- rphylopic::get_phylopic(uuid="61e2982f-07d4-44db-ac00-27a2d86ca855") #Fucus serratus
	crustacea <- rphylopic::get_phylopic(uuid="7197c71a-0653-4e82-bcbb-b156c150826a") #Callinectes sapidus
	echinoderms <- rphylopic::get_phylopic(uuid="0adc540a-d1c1-4dd8-bc4d-e258077c5dbd") #Asteriidae
	bony_fishes <- rphylopic::get_phylopic(uuid="ed245bc2-c8dc-4ad4-9d08-efa19fe6854a") #Sebastes diaconus
	chondrichthyes <- rphylopic::get_phylopic(uuid="8d9410d3-8b41-4c4b-ad27-1bf8dd606e47") #Bathyraja griseocauda
	cnidarians <- rphylopic::get_phylopic(uuid="927676f4-7ea9-4aa6-8e0a-3c70d8253f89") #Nematostella vectensis
	sauropsida <- rphylopic::get_phylopic(uuid="86334821-42ec-4da1-bb9d-53f3d6941c77") #Pygoscelis papua
	molluscs1 <- rphylopic::get_phylopic(uuid="90449630-774d-48dc-b16c-912b69825dee") #Pyroteuthis margaritifera
	molluscs2 <- rphylopic::get_phylopic(uuid="6c2e67f0-14e7-4ba0-ba73-2420cacfa9a3") #Tricolia pullus
	rphylopic::add_phylopic_base(img=bony_fishes,x=xleft,y=ytop,ysize=ysize-2,fill=cladeCols[1])
		text(x=xleft+legendshift,y=ytop,labels="Bony fishes",pos=pos,cex=txt.cex)
	rphylopic::add_phylopic_base(img=sauropsida,x=xleft,y=ytop-dshift,ysize=ysize,fill=cladeCols[2])
		text(x=xleft+legendshift,y=ytop-1*dshift,labels="Birds",pos=pos,cex=txt.cex)
	rphylopic::add_phylopic_base(img=mammal,x=xleft,y=ytop-2*dshift,ysize=ysize,fill=cladeCols[3])
		text(x=xleft+legendshift,y=ytop-2*dshift,labels="Mammals",pos=pos,cex=txt.cex)
	rphylopic::add_phylopic_base(img=chondrichthyes,x=xleft,y=ytop-3*dshift,ysize=ysize,fill=cladeCols[4])
		text(x=xleft+legendshift,y=ytop-3*dshift,labels="Cartilaginous Fishes",pos=pos,cex=txt.cex)
	rphylopic::add_phylopic_base(img=echinoderms,x=xleft,y=ytop-4*dshift,ysize=ysize,fill=cladeCols[5])
		text(x=xleft+legendshift,y=ytop-4*dshift,labels="Echinoderms",pos=pos,cex=txt.cex)
	rphylopic::add_phylopic_base(img=molluscs1,x=xleft,y=ytop-5*dshift,ysize=ysize,fill=cladeCols[6])
		text(x=xleft+legendshift,y=ytop-5*dshift,labels="Molluscs",pos=pos,cex=txt.cex)
	rphylopic::add_phylopic_base(img=crustacea,x=xleft,y=ytop-6*dshift,ysize=ysize,fill=cladeCols[7])
		text(x=xleft+legendshift,y=ytop-6*dshift,labels="Crustaceans",pos=pos,cex=txt.cex)
	rphylopic::add_phylopic_base(img=cnidarians,x=xleft,y=ytop-7*dshift,ysize=ysize+1,fill=cladeCols[8])
		text(x=xleft+legendshift,y=ytop-7*dshift,labels="Cnidarians",pos=pos,cex=txt.cex)
	rphylopic::add_phylopic_base(img=vasc_plant,x=xleft,y=ytop-8*dshift,ysize=ysize+1,fill=cladeCols[9])
		text(x=xleft+legendshift,y=ytop-8*dshift,labels="Vascular Plants",pos=pos,cex=txt.cex)
#	rphylopic::add_phylopic_base(img=ochrophyta,x=50,y=5,ysize=ysize,fill=cladeCols[10])
	#rphylopic::add_phylopic_base(img=molluscs2,x=5,y=5,ysize=1,fill=cladeCols[1])
}

gussyUpPhyloFig <- function(tree,cladeCols,rounding){
	cladeCols1 <- adjustcolor(cladeCols,0.5)
	highlightClade(tree=tree,dSp1="Engraulis encrasicolus",dSp2="Sebastes diaconus",col=cladeCols1[1],x0=NULL,x1=NULL,y0=NULL,y1=NULL,rounding=rounding)
	highlightClade(tree=tree,dSp1="Pygoscelis papua",dSp2="Aythya marila",col=cladeCols1[2],x0=NULL,x1=NULL,y0=NULL,y1=NULL,rounding=rounding)
	highlightClade(tree=tree,dSp1="Phocoena sinus",dSp2="Halichoerus grypus atlantica",col=cladeCols1[3],x0=NULL,x1=NULL,y0=NULL,y1=NULL,rounding=rounding)
	highlightClade(tree=tree,dSp1="Bathyraja panthera",dSp2="Sphyrna tiburo",col=cladeCols1[4],x0=NULL,x1=NULL,y0=NULL,y1=NULL,rounding=rounding)
	highlightClade(tree=tree,dSp1="Paracentrotus lividus",dSp2="Pisaster ochraceus",col=cladeCols1[5],x0=NULL,x1=NULL,y0=NULL,y1=NULL,rounding=rounding)
	highlightClade(tree=tree,dSp1="Chlorostoma funebralis",dSp2="Pyroteuthis margaritifera",col=cladeCols1[6],x0=NULL,x1=NULL,y0=NULL,y1=NULL,rounding=rounding)
	highlightClade(tree=tree,dSp1="Callinectes sapidus",dSp2="Penaeus duorarum",col=cladeCols1[7],x0=NULL,x1=NULL,y0=NULL,y1=NULL,rounding=rounding)
	highlightClade(tree=tree,dSp1="Acropora palmata",dSp2="Ectopleura larynx",col=cladeCols1[8],x0=NULL,x1=NULL,y0=NULL,y1=NULL,rounding=rounding)
	highlightClade(tree=tree,dSp1="Rhizophora mangle",dSp2="Laguncularia racemosa",col=cladeCols1[9],x0=NULL,x1=NULL,y0=NULL,y1=NULL,rounding=rounding)
	addPhylopics(cladeCols)
}

addCladeLine <- function(tree,dSp1,dSp2,col,xcoord,lwd=4,lend=2){
	m <- ape::getMRCA(tree,tip=c(dSp1,dSp2))
	decs <- getDescendants(tree,m)
	yrange <- range(ape::node.height(tree)[decs])
	segments(x0=xcoord,y0=yrange[1],
			 x1=xcoord,y1=yrange[2],
			 col=col,lwd=lwd,lend=lend)
}

gussyUpDivBarPlot <- function(tree,xcoord,cladeCols){
	#recover()
	addCladeLine(tree=tree,dSp1="Engraulis encrasicolus",dSp2="Sebastes diaconus",col=cladeCols[1],xcoord,lwd=3,lend=2)
	addCladeLine(tree=tree,dSp1="Pygoscelis papua",dSp2="Aythya marila",col=cladeCols[2],xcoord,lwd=3,lend=2)
	addCladeLine(tree=tree,dSp1="Phocoena sinus",dSp2="Halichoerus grypus atlantica",col=cladeCols[3],xcoord,lwd=3,lend=2)
	addCladeLine(tree=tree,dSp1="Bathyraja panthera",dSp2="Sphyrna tiburo",col=cladeCols[4],xcoord,lwd=3,lend=2)
	addCladeLine(tree=tree,dSp1="Paracentrotus lividus",dSp2="Pisaster ochraceus",col=cladeCols[5],xcoord,lwd=3,lend=2)
	addCladeLine(tree=tree,dSp1="Chlorostoma funebralis",dSp2="Pyroteuthis margaritifera",col=cladeCols[6],xcoord,lwd=3,lend=2)
	addCladeLine(tree=tree,dSp1="Callinectes sapidus",dSp2="Penaeus duorarum",col=cladeCols[7],xcoord,lwd=3,lend=2)
	addCladeLine(tree=tree,dSp1="Acropora palmata",dSp2="Ectopleura larynx",col=cladeCols[8],xcoord,lwd=3,lend=2)
	addCladeLine(tree=tree,dSp1="Rhizophora mangle",dSp2="Laguncularia racemosa",col=cladeCols[9],xcoord,lwd=3,lend=2)	
}

addLegend <- function(top,xl,xr,txtShft,txtCex){
	grayCols <- gray(c(0.2,0.5,0.8),1)
	gradCols <- rev(colorRampPalette(RColorBrewer::brewer.pal(9,"Blues"))(150)[51:150])
	colRast <- as.raster(matrix(gradCols,nrow=length(gradCols),ncol=1))
	rect(xleft=xl,ybottom=top-3,xright=xr,ytop=top-1,col=grayCols[3])
		text(x=xr+txtShft,y=top-2,labels="present",pos=4,cex=txtCex)
	rect(xleft=xl,ybottom=top-6,xright=xr,ytop=top-4,col=grayCols[2])
		text(x=xr+txtShft,y=top-5,labels="intermediate",pos=4,cex=txtCex)
	rect(xleft=xl,ybottom=top-9,xright=xr,ytop=top-7,col=grayCols[1])
		text(x=xr+txtShft,y=top-8,labels="absent",pos=4,cex=txtCex)
	rasterImage(colRast,xleft=xl,ybottom=top-25,xright=xr,ytop=top-10.5,interpolate=TRUE)
	rect(xleft=xl,ybottom=top-25,xright=xr,ytop=top-10.5)
		text(x=xr+txtShft,y=top-11.5,labels="high",pos=4,cex=txtCex)
		text(x=xr+txtShft,y=top-24,labels="low",pos=4,cex=txtCex)
}

phyViz <- function(db,fit,XX,predNames,tipcols,tree,xlim,valRange=NULL,adj=2,logX=FALSE,rounding=0.1){
	#recover()
	spNames <- row.names(db$relMat)
	tree <- ape::keep.tip(tree,gsub("_"," ",spNames))
	tree <- ape::ladderize(tree,right=FALSE)
	tipOrder <- getTipOrder(tree)
	nSp <- length(tipOrder)
	spOrder <- match(tipOrder,spNames)
	XX <- XX[spOrder,]
	par(oma=c(6,0.25,2,0))
	layout(matrix(c(rep(1,30),rep(2,20),rep(3,20)),10,7))
	ape::plot.phylo(tree,cex=0.7,no.margin=TRUE, #label.offset=50
					x.lim=c(0,1320),tip.color=tipcols,edge.width=1, #x.lim=c(0,1950) #1700
					show.tip.label=FALSE)
	cladeCols <- RColorBrewer::brewer.pal(10,"Paired")
	gussyUpPhyloFig(tree,cladeCols=cladeCols,rounding)
	axisPhylo()
		mtext(side=1,text="Time (Mya)",padj=3)
		mtext(side=3,text="Time-calibrated Phylogeny",padj=1,cex=1.4)
	plot(0,type='n',
				xlim=c(0,1),
				yaxt='n',xaxt='n',xlab="",ylab="",bty='n',
				ylim=c(1,nSp))
	vizPredMat(xx=XX,N=nSp,xnames=predNames,xxext=c(0,1),yyext=c(0,nSp))
	if(logX){
		xlim <- range(log(db$Y)) + c(-4,1)
		y <- log(db$Y)
	} else {
		y <- db$Y
	}
		mtext(side=3,text="Biotic/Abiotic Predictors",padj=1,cex=1.4)
	plot(0,type='n',
			xlim=xlim,
			yaxt='n',xaxt='n',xlab="",ylab="",bty='n',
			ylim=c(1,nSp)) #c(0,(length(cnDens)+1)*1.25))
	axis(side=1,
		at=c(xlim[1],(xlim[2]+xlim[1])/2,xlim[2]),
		labels=round(c(xlim[1],(xlim[2]+xlim[1])/2,xlim[2]),2))
	#text(x=xtext,y=(0.4 + 1:nSp)*1.25,labels=spNames[spOrder],srt=0)
	#abline(h=spOrder-0.5,lty=3,col="gray")
	segments(x0=0,
			 y0=c(1:db$N)-0.25,
			 x1=y[spOrder],
			 y1=c(1:db$N)-0.25,
			 col="red3",
			 lwd=5)
	gussyUpDivBarPlot(tree=tree,xcoord=-0.001,cladeCols=cladeCols)
		mtext(side=1,text="Diversity",padj=3)
		mtext(side=3,text="Species-level Genetic Diversity",padj=1,cex=1.4)
	addLegend(top=84,xl=0.02,xr=0.0235,txtShft=7e-4,txtCex=1.4)
}

vizPredMat <- function(xx,N,xnames,xxext,yyext){
#	recover()
	npreds <- ncol(xx)
	xColLefts <- seq(xxext[1],xxext[2]-diff(xxext)/npreds,length.out=npreds)
	xColRights <- seq(xxext[1]+diff(xxext)/npreds,xxext[2],length.out=npreds)
	invisible(
		lapply(1:npreds,function(i){
			vizPredCol(xx[,i],N,yext=c(yyext),xext=c(xColLefts[i],xColRights[i]))
		})
	)
	text(pos=2,srt=90,x=xColRights-0.02,y=par("usr")[3]+1.5,xpd=NA,labels=xnames)
#	text(pos=2,srt=90,x=xColRights,y=yyext[1],labels=xnames)
#	axis(side=1,srt=45,las=2,at=xColRights,labels=xnames)
}

vizPredCol <- function(x,N,yext,xext){
	#recover()
	valRange <- range(x,na.rm=TRUE)
	if(length(unique(x)) > 4){
		allCols <- colorRampPalette(RColorBrewer::brewer.pal(9,"Blues"))(150)[51:150]
		#wesanderson::wes_palette("Zissou1", 100, type = "continuous")
		#colorspace::diverge_hcl(100,c(246,40),c=96)
		#viridis::viridis_pal(alpha=1,begin=0,end=1,direction=1,option="D")(100)
	} else {
		allCols <- gray(c(0.2,0.5,0.8),1)
	}
	cols <- colFunc(x,cols=allCols,nCols=100,valRange=valRange)
	xls <- rep(xext[1],N)
	xrs <- rep(xext[2],N)
	ybs <- seq(yext[1],yext[2]-diff(yext)/N,length.out=N)
	yts <- seq(yext[1]+diff(yext)/N,yext[2],length.out=N)
	invisible(lapply(1:N,function(i){
		rect(xls[i],ybs[i],xrs[i],yts[i],col=cols[i],border=NA)
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

discreteViolPlot <- function(z,nPreds,predName,xAxLabs,logY=FALSE){
	x <- z[[predName]]
	if(nPreds==2){
		vCols <- gray(c(0.2,0.5,0.8),1)[c(1,3)]
		xRange <- c(-0.5,1.5)
	} else if (nPreds==3){
		vCols <- gray(c(0.2,0.5,0.8),1)
		xRange <- c(-0.5,2.5)
	}
	if(logY){
		yRange <- range(log(z$div)) + c(-0.1+0.1)
		y <- log(z$div)
		yLab <- "Log(Genetic Diversity)"
	} else {
		yRange <- range(z$div) + c(-0.003+0.003)
		y <- z$div
		yLab <- "Genetic Diversity"
	}
	plot(0,type='n',xlab="",xaxt='n',ylab=yLab,xlim=xRange,ylim=yRange)
	for(i in 1:nPreds){
		vioplot::vioplot(at=i-1,y[x==(i-1)],add=TRUE,col=vCols[i])
	}
	points(jitter(x),y=y,col="black",bg=z$cladecolor,pch=21,cex=1.7)
	axis(side=1,at=0:(nPreds-1),labels=xAxLabs)
}
