predictors <- rbind(z[["meanlat.gbif"]],					# abiotic
					z[["n_ECOREGIONS.all"]],				# abiotic
					z[["max.95.sea.gbif"]]/max(z[["max.95.sea.gbif"]]),				# abiotic
					z[["Body_Size"]],				# biotic
					z[["Fecundity_EggSize"]],				# biotic
					z[["Generational_Structure"]],				# biotic
					z[["ReturnToSpawningGround"]],				# biotic
					z[["Spawning_mode"]],				# biotic
					z[["Larval_feeding"]],				# biotic
					z[["PLD_point2"]],				# biotic
					z[["isPlanktonic_atanypoint"]],				# biotic
					z[["ratio.sea.95"]],					# nuisance
					z[["n_samples"]],					# nuisance
					z[["mean_raw_read_cnt"]]/1e7,					# nuisance
					z[["read_length"]],					# nuisance
					z[["mean_locus_depth"]],					# nuisance
					z[["n.totalsnps"]]/1e5)					# nuisance



predNames <- c("mean species latitude","number of ecoregions",
			   "range extent","body size",
			   "egg size","generational structure",
			   "philopatry","spawning mode",
			   "larval feeding","pelagic larval duration",
			   "planktonicity","ratio of range sampled : total",
			   "number of samples","mean raw read count",
			   "read length","mean locus depth","number of SNPs")

nIter <- 1e4

################################
# analyze s with one predictor at a time
#	beta model, unscaled s
################################
db <- lapply(1:nrow(predictors),
			function(i){
				makeDB(predictors[i,],1-z$s,phyStr)
		})
names(db) <- predNames

cl <- parallel::makeCluster(12)
doParallel::registerDoParallel(cl)

fits <- foreach::foreach(i = 1:nrow(predictors)) %dopar% {
	message(sprintf("now analyzing predictor %s/%s",i,nrow(predictors)))
	rstan::sampling(object=betaPhyReg,
					data=db[[i]],
					iter=nIter,
					thin=nIter/500,
					chains=1,
					control=setNames(list(15),"max_treedepth"))
}

parallel::stopCluster(cl)

names(fits) <- predNames
out <- list("db"=db,"fits"=fits)
save(out,file="div.Robj")

pdf(file="div.pdf",width=14,height=12)
	par(mfrow=c(4,5)) ; for(i in 1:nrow(predictors)){betaPPS(out$db[[i]],out$fit[[i]],500,predName=names(out$fit)[i])}
dev.off()

pdf(file="div_phyloFit.pdf",width=12,height=12)
	for(i in 1:nrow(predictors)){
		modAdViz(out$db[[i]],out$fit[[i]],predName=predNames[i],
					nPPS=500,tree=sampPhy,xlim=c(0,0.05),
					valRange=NULL,qnt=0.999,adj=0.5)
	}
dev.off()

pdf(file="div_phyloFit_log.pdf",width=12,height=12)
	for(i in 1:nrow(predictors)){
		modAdViz(out$db[[i]],out$fit[[i]],predName=predNames[i],
					nPPS=500,tree=sampPhy,xlim=c(0,0.05),
					valRange=NULL,qnt=0.999,adj=0.5,logX=TRUE)
	}
dev.off()

pdf(file="div_betas.pdf",width=14,height=10)
	postBetaPlot(out=out,predNames=names(out$fits),reorder=TRUE,cols=NULL,stdize=TRUE,qnt=0.99)
dev.off()

################################
# analyze diversity with one biological predictor 
#	and all the "nuisance" parameters
#	beta model, unscaled diversity
################################
bioPreds <- predictors[1:11,]
bioPredNames <- predNames[1:11]
nuisPreds <- predictors[13:16,]
nuisPredNames <- predNames[13:16]

db <- lapply(1:nrow(bioPreds),
			function(i){
				makeMultiDB(list(predictors[i,],
								 nuisPreds[1,],
								 nuisPreds[2,],
								 nuisPreds[3,],
								 nuisPreds[4,]),
							 Y=1-z$s,
							 phyStr=phyStr)
		})
names(db) <- bioPredNames

cl <- parallel::makeCluster(12)
doParallel::registerDoParallel(cl)

fits <- foreach::foreach(i = 1:nrow(bioPreds)) %dopar% {
	message(sprintf("now analyzing multiple-linear regression with main predictor %s/%s",i,nrow(predictors)))
	rstan::sampling(object=betaPhyReg,
					data=db[[i]],
					iter=nIter,
					thin=nIter/500,
					chains=1,
					control=setNames(list(15),"max_treedepth"))
}

parallel::stopCluster(cl)

names(fits) <- bioPredNames
out <- list("db"=db,"fits"=fits)
save(out,file="div_multiPred.Robj")

pdf(file="div_multiPred.pdf",width=12,height=10)
	par(mfrow=c(3,4)) ; for(i in 1:nrow(bioPreds)){betaPPS(out$db[[i]],out$fit[[i]],500,predName=names(out$fit)[i],multiPred=TRUE)}
dev.off()

pdf(file="div_multiPred_phyloFit.pdf",width=12,height=12)
	for(i in 1:nrow(bioPreds)){
		modAdViz(out$db[[i]],out$fit[[i]],predName=bioPredNames[i],
					nPPS=500,tree=sampPhy,xlim=c(0,0.05),
					valRange=NULL,qnt=0.999,adj=0.5)
	}
dev.off()

pdf(file="div_multiPred_phyloFit_log.pdf",width=12,height=12)
	for(i in 1:nrow(bioPreds)){
		modAdViz(out$db[[i]],out$fit[[i]],predName=bioPredNames[i],
					nPPS=500,tree=sampPhy,xlim=c(0,0.05),
					valRange=NULL,qnt=0.999,adj=0.5,logX=TRUE)
	}
dev.off()


pdf(file="div_multiPred_betas.pdf",width=14,height=10)
	postBetaPlot(out=out,predNames=names(out$fits),reorder=TRUE,cols=NULL,stdize=TRUE,multiPred=TRUE)
dev.off()


################################
# analyze diversity with one biological predictor 
#	and all the "nuisance" parameters
#	beta model, unscaled diversity
#	IN SPECIFIC GROUPS:
#		PLD w/ isPlanktonic filter
#		body size just in fish
#		all analyses dropping birds & mammals
################################

# PLD w/ isPlanktonic filter
isPlanktonic <- which(predictors[11,]==1)
db <- makeMultiDB(list(predictors[10,isPlanktonic],
						 nuisPreds[1,isPlanktonic],
						 nuisPreds[2,isPlanktonic],
						 nuisPreds[3,isPlanktonic],
						 nuisPreds[4,isPlanktonic]),
					 Y=1-z$s[isPlanktonic],
					 phyStr=phyStr[isPlanktonic,isPlanktonic])

fit <- rstan::sampling(object=betaPhyReg,
						data=db,
						iter=nIter,
						thin=nIter/500,
						chains=1,
						control=setNames(list(15),"max_treedepth"))

pdf(file="div_PLD_isPlanktonic.pdf",width=12,height=10)
	betaPPS(db,fit,500,predName="pelagic larval dispersal\nin planktonic species",multiPred=TRUE)
dev.off()
out <- list("db"=db,"fit"=fit)
save(out,file="PLD_isPlanktonic.Robj")


# BodySize w/ isFish filter
fishMRCA <- ape::getMRCA(sampPhy,c(which(sampPhy$tip.label=="Engraulis encrasicolus"),which(sampPhy$tip.label=="Sebastes diaconus")))
isFish <- phangorn::Descendants(sampPhy,fishMRCA,type="tips")[[1]]

db <- makeMultiDB(list(predictors[4,isFish],
						 nuisPreds[1,isFish],
						 nuisPreds[2,isFish],
						 nuisPreds[3,isFish],
						 nuisPreds[4,isFish]),
					 Y=1-z$s[isFish],
					 phyStr=phyStr[isFish,isFish])

fit <- rstan::sampling(object=betaPhyReg,
						data=db,
						iter=nIter,
						thin=nIter/500,
						chains=1,
						control=setNames(list(15),"max_treedepth"))

pdf(file="div_bodySize_isFish.pdf",width=12,height=10)
	betaPPS(db,fit,500,predName="body size\nin fishes",multiPred=TRUE)
dev.off()
out <- list("db"=db,"fit"=fit)
save(out,file="bodySize_isFish.Robj")


# all analyses with only primarily marine spp.
ogMarine <- which(z$taxclade %in% c("cnidarians","crustacea","molluscs","echinoderms","chondrichthyes","bony fishes"))

db <- lapply(1:nrow(bioPreds),
			function(i){
				makeMultiDB(list(predictors[i,ogMarine],
								 nuisPreds[1,ogMarine],
								 nuisPreds[2,ogMarine],
								 nuisPreds[3,ogMarine],
								 nuisPreds[4,ogMarine]),
							 Y=1-z$s[ogMarine],
							 phyStr=phyStr[ogMarine,ogMarine])
		})
names(db) <- bioPredNames

cl <- parallel::makeCluster(12)
doParallel::registerDoParallel(cl)

fits <- foreach::foreach(i = 1:nrow(bioPreds)) %dopar% {
	message(sprintf("now analyzing multiple-linear regression with main predictor %s/%s",i,nrow(predictors)))
	rstan::sampling(object=betaPhyReg,
					data=db[[i]],
					iter=nIter,
					thin=nIter/500,
					chains=1,
					control=setNames(list(15),"max_treedepth"))
}

parallel::stopCluster(cl)

names(fits) <- bioPredNames
out <- list("db"=db,"fits"=fits)
save(out,file="div_multiPred_ogMarine.Robj")

pdf(file="div_multiPred_ogMarine.pdf",width=12,height=10)
	par(mfrow=c(3,4)) ; for(i in 1:nrow(bioPreds)){betaPPS(out$db[[i]],out$fit[[i]],500,predName=names(out$fit)[i],multiPred=TRUE)}
dev.off()

pdf(file="div_multiPred_phyloFit_ogMarine.pdf",width=12,height=12)
	for(i in 1:nrow(bioPreds)){
		modAdViz(out$db[[i]],out$fit[[i]],predName=bioPredNames[i],
					nPPS=500,tree=sampPhy,xlim=c(0,0.05),
					valRange=NULL,qnt=0.999,adj=0.5)
	}
dev.off()

pdf(file="div_multiPred_phyloFit_log_ogMarine.pdf",width=12,height=12)
	for(i in 1:nrow(bioPreds)){
		modAdViz(out$db[[i]],out$fit[[i]],predName=bioPredNames[i],
					nPPS=500,tree=sampPhy,xlim=c(0,0.05),
					valRange=NULL,qnt=0.999,adj=0.5,logX=TRUE)
	}
dev.off()

if(FALSE){
################################
# analyze Nbhd with one predictor at a time
#	exp model
################################

db <- lapply(1:nrow(predictors),
			function(i){
				makeDB(predictors[i,],z$nbhd,phyStr)
		})
names(db) <- predNames

cl <- parallel::makeCluster(12)
doParallel::registerDoParallel(cl)

fits <- foreach::foreach(i = 1:nrow(predictors)) %dopar% {
	message(sprintf("now analyzing predictor %s/%s",i,nrow(predictors)))
	rstan::sampling(object=expPhyReg,
					data=db[[i]],
					iter=nIter,
					thin=nIter/500,
					save_warmup=FALSE,
					chains=1,
					control=setNames(list(15),"max_treedepth"))
}

parallel::stopCluster(cl)

names(fits) <- predNames
out <- list("db"=db,"fits"=fits)
save(out,file="nbhd.Robj")

pdf(file="nbhd.pdf",width=12,height=12)
	par(mfrow=c(4,4)) ; for(i in 1:nrow(predictors)){expPPS(out$db[[i]],out$fit[[i]],500,predName=names(out$fit)[i])}
dev.off()

pdf(file="nbhd_betas.pdf",width=14,height=10)
	postBetaPlot(out=out,predNames=names(out$fits),reorder=TRUE,cols=NULL,stdize=TRUE,multiPred=FALSE)
dev.off()

################################
# analyze Nbhd size with one biological predictor 
#	and all the "nuisance" parameters
################################
bioPreds <- predictors[1:11,]
bioPredNames <- predNames[1:11]
nuisPreds <- predictors[13:16,]
nuisPredNames <- predNames[13:16]

db <- lapply(1:nrow(bioPreds),
			function(i){
				makeMultiDB(list(predictors[i,],
								 nuisPreds[1,],
								 nuisPreds[2,],
								 nuisPreds[3,],
								 nuisPreds[4,]),
							 Y=z$nbhd,
							 phyStr=phyStr)
		})
names(db) <- bioPredNames

cl <- parallel::makeCluster(12)
doParallel::registerDoParallel(cl)

fits <- foreach::foreach(i = 1:nrow(bioPreds)) %dopar% {
	message(sprintf("now analyzing multiple-linear regression with main predictor %s/%s",i,nrow(predictors)))
	rstan::sampling(object=expPhyReg,
					data=db[[i]],
					iter=nIter,
					thin=nIter/500,
					chains=1,
					control=setNames(list(15),"max_treedepth"))
}

parallel::stopCluster(cl)

names(fits) <- bioPredNames
out <- list("db"=db,"fits"=fits)
save(out,file="nbhd_multiPred.Robj")

pdf(file="nbhd_multiPred.pdf",width=12,height=10)
	par(mfrow=c(3,4)) ; for(i in 1:nrow(bioPreds)){expPPS(out$db[[i]],out$fit[[i]],500,predName=names(out$fit)[i],multiPred=TRUE)}
dev.off()

pdf(file="nbhd_multiPred_betas.pdf",width=14,height=10)
	postBetaPlot(out=out,predNames=names(out$fits),reorder=TRUE,cols=NULL,stdize=TRUE,multiPred=TRUE)
dev.off()
}

################
# GRAVEYARD
################

if(FALSE){
################################
# analyze s with multiple biological predictors 
#	and all the "nuisance" parameters
#	beta model, unscaled diversity
################################

# for ecoregions & range size, include both

db <- makeMultiDB(list(predictors[2,],
					   predictors[3,],
					   nuisPreds[1,],
					   nuisPreds[2,],
					   nuisPreds[3,],
					   nuisPreds[4,]),
				   Y=1-z$s,
				   phyStr=phyStr)

fit <- rstan::sampling(object=betaPhyReg,
						data=db,
						iter=nIter,
						thin=nIter/500,
						chains=1,
						control=setNames(list(15),"max_treedepth"))

out <- list("db"=db,"fit"=fit)
save(out,file="ecoreg_range_s.Robj")

# betaPPS(out$db,out$fit,500,predName="ecoreg",multiPred=TRUE)
# b1 <- extract(out$fit,"beta[1]",inc_warmup=TRUE,permute=FALSE)
# b2 <- extract(out$fit,"beta[2]",inc_warmup=TRUE,permute=FALSE)
# plot(b1,b2,xlim=c(-0.03,0.03))
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

vizSinglePredictorFit <- function(db,fit,xName,yName,mName=NULL,cols,unscale=NULL,...){
	b <- rstan::extract(fit,"beta",permute=FALSE)
	sig <- 1-2*(abs(0.5-ecdf(b)(0)))
	siglwd <- ifelse(sig<0.05,4,0.5)
	m <- rstan::extract(fit,"mu",permute=FALSE)
	y <- db$Y
	if(is.null(mName)){
		mName <- yName
	}
	if(!is.null(unscale)){
		y <- db$Y*unscale$mx + unscale$mn
		m <- m*unscale$mx + unscale$mn
		b <- b*unscale$mx
	}
	plot(db$X,y,xlab="",ylab=yName,
		main=sprintf("Effect of %s on %s",xName,mName),
		pch=19,col=cols,cex=2,...)
		abline(mean(m),mean(b),lwd=siglwd)
	print(sig)
}

par(mfrow=c(3,4))
lapply(1:nrow(predictors),
	function(i){
		vizSinglePredictorFit(db=dbns[[i]],
							  fit=nfits[[i]],
							  xName=predNames[i],
							  yName="neighborhood size",
							  cols=z$cladecol,unscale=NULL)
	})

# make some cute figures for the poster
pdf(file="../figures/s_larvalFeeding.pdf",width=8,height=6,pointsize=12)
	par(cex.lab=2,cex.main=2,cex.axis=1.5,mar=c(6,6,5,3))
	vizSinglePredictorFit(db=dbss[[9]],
						  fit=sfits[[9]],
						  xName="larval feeding",
						  yName=expression(paste("1 - ",pi)),
						  mName="genetic diversity",
						  cols=adjustcolor(z$cladecol,0.7),
						  unscale=unscale,
						  xaxt='n')
	axis(side=1,at=c(0,1,2),labels=c("not planktonic","lecithotrophic","planktotrophic"),padj=1)
dev.off()

pdf(file="../figures/nbhd_spawnMode.pdf",width=8,height=6,pointsize=12)
	par(cex.lab=2,cex.main=2,cex.axis=1.5,mar=c(6,6,5,3))
	vizSinglePredictorFit(db=dbns[[8]],
						  fit=nfits[[8]],
						  xName="spawning mode\n",
						  yName="neighborhood size",
						  cols=adjustcolor(z$cladecol,0.7),
						  unscale=NULL,
						  xaxt='n')
	axis(side=1,at=c(0,1,2),labels=c("internal\nfertilization","nesting","free-floating"),padj=1)
dev.off()
################################
# run analyses with all predictors
################################

# modeling collecting phase pi
s <- (z$s.wish - min(z$s.wish))/max((z$s.wish - min(z$s.wish)))
db_s <- list("N" = nrow(phyStr),
		   	 "Y" = s,
		   	 "nX"= nrow(predictors),
		   	 "X" = predictors,
		     "relMat" = phyStr)

fit_s <- sampling(object=mvn,
					data=db_s,
					iter=3e3,
					chains=2,
					control=setNames(list(15),"max_treedepth"))

# modeling neighborhood sizes
db_nbhd <- list("N" = nrow(phyStr),
		   	    "Y" = z$nbhd.wish,
		   	    "nX"= nrow(predictors),
		   	    "X" = predictors,
		        "relMat" = phyStr)

fit_nbhd <- sampling(object=mvn,
						data=db_nbhd,
						iter=3e3,
						chains=2,
						control=setNames(list(15),"max_treedepth"))
save(db_s,fit_s,file="mod_s.Robj")
save(db_nbhd,fit_nbhd,file="mod_nbhd.Robj")

################################
# plot output
################################

plotMarginals <- function(b,predNames,mf){
	nPred <- length(predNames)
	beta.dens <- lapply(1:nPred,function(i){density(b[,i])})
	par(mfrow=mf,mar=c(2,2,2,2))
	for(i in 1:nPred){
		plot(0,
			 main = predNames[i],
			 xlim=range(c(0,b[,i])),ylim=range(beta.dens[[i]]$y),
			 xlab="",cex.main=1,cex.axis=1,cex.lab=1,type='n',
			 ylab="Density")
			polygon(x=c(0,beta.dens[[i]]$x,0),
					y=c(0,beta.dens[[i]]$y,0),
					col=adjustcolor("gray",0.5))
			abline(v=0,lwd=4,lty=2,col=1)
			sig <- demarcateCIs(b[,i],beta.dens[[i]])
		if(sig){
			box(lwd=3,col="red")
		}
	}
}

demarcateCIs <- function(z,z.dens){
	CI <- quantile(z,c(0.025,0.975))
	segments(x0=CI[1],x1=CI[1],
			 y0=0,y1=max(z.dens$y)/2,
			 col=4,lwd=4)
	segments(x0=CI[2],x1=CI[2],
			 y0=0,y1=max(z.dens$y)/2,
			 col=4,lwd=4)
	sig <- TRUE
	if(0 > CI[1] & 0 < CI[2]){
		sig <- FALSE
	}
	return(sig)
}

plotRankBetas <- function(b,predNames){
	mn <- colMeans(b)
	CIs <- apply(b,2,function(x){quantile(x,c(0.025,0.975))})
	CImin <- CIs[1,]
	CImax <- CIs[2,]
	rdr <- order(mn)
	plot(mn[rdr],ylab="effect size",xaxt='n',xlab="",ylim=range(c(CImin,CImax)))
		points(CImax[rdr],pch=25,pt.col="black",bg="black")
		points(CImin[rdr],pch=24,pt.col="black",bg="black")
		abline(h=0,lty=2)
		axis(side=1,at=1:length(predNames),labels=predNames[rdr],las=2)
}


b_s <- rstan::extract(fit_s,"beta",permute=FALSE)
mu_s <- rstan::extract(fit_s,"mu",permute=FALSE)
b_nbhd <- rstan::extract(fit_nbhd,"beta",permute=FALSE)
mu_nbhd <- rstan::extract(fit_nbhd,"mu",permute=FALSE)


pdf(file="divdiv_s.pdf",width=4,height=8)
	plotMarginals(b_s[,1,],predNames,c(6,2))
		plot(0,type='n',bty='n',yaxt='n',xaxt='n',xlab="",ylab="")
		legend(x="topright",col=c(4,1,"red"),lty=c(1,2,1),lwd=3,legend=c("95% CI","Zero","p < 0.05"),cex=1.5,bty='n')
dev.off()

pdf(file="divdiv_nbhd.pdf",width=4,height=8)
	plotMarginals(b_nbhd[,1,],predNames,c(6,2))
		plot(0,type='n',bty='n',yaxt='n',xaxt='n',xlab="",ylab="")
		legend(x="topright",col=c(4,1,"red"),lty=c(1,2,1),lwd=3,legend=c("95% CI","Zero","p < 0.05"),cex=1.5,bty='n')
dev.off()

################################
# visualize model fit
################################


parMeans <- extract(fit_nbhd,"meanVec",permute=FALSE,inc_warmup=TRUE)[,1,]
res <- db_nbhd$Y - colMeans(parMeans)
mu <- extract(fit_nbhd,"mu",permute=FALSE,inc_warmup=TRUE)[,1,1]
pmci <- apply(parMeans,2,function(x){quantile(x,c(0.005,0.995))})
par(mfrow=c(3,4))
for(i in 1:11){

plot(db_nbhd$X[i,],res+mean(b_nbhd[,,i])*db_nbhd$X[i,],main=predNames[i],ylim=range(pmci))
	segments(x0=db_nbhd$X[i,],x1=db_nbhd$X[i,],
			 y0=pmci[1,],y1=pmci[2,],col=z$cladecol)
	points(db_nbhd$X[i,],z$nbhd.wish,pch=19,col=z$cladecolor)
	abline(mean(mu),colMeans(b_nbhd[,1,])[i],col="black",lty=2,lwd=2)
}

plot(0,col=z$cladecol,pch=19,cex=2,type='n',xlim=c(0,db_nbhd$N),ylim=range(pmci),ylab="",xlab="")
	segments(x0=1:db_nbhd$N,x1=1:db_nbhd$N,
			 y0=pmci[1,],y1=pmci[2,],col=z$cladecol)
	points(1:db_nbhd$N,db_nbhd$Y,pch=19,cex=2,col=z$cladecol)
#	abline(0,1,col="black",lty=2,lwd=2)

bspt <- colMeans(b_s[,1,])
bnpt <- colMeans(b_nbhd[,1,])

plot(bnpt,type='n') ; text(bnpt,labels=predNames)
	plot(predictors[11,],s,pch=19,col=adjustcolor(1,0.3))
		abline(mean(mu_s),b=bspt[11])

################################
# visualize data
################################

pdf(file="../figures/pi_across_spp.pdf",width=10,height=5,pointsize=12)
par(mar=c(10,5,3,1))
	plot(z$s.wish[order(z$s.wish)],
			xlab="",ylab=expression(paste("1 - deep-time ",pi)),
			xaxt='n',pch=19,cex=2,main="Genetic diversity",
			ylim=range(z$s.wish)+c(-0.001,0.001),col=adjustcolor(z$cladecolor[order(z$s.wish)],0.5))
	text(1:length(z$s.wish),par("usr")[3],labels=gsub("_"," ",z$species[order(z$s.wish)]),srt=50,adj=c(1.1,1.1),xpd=TRUE)
dev.off()

pdf(file="../figures/nbhd_across_spp.pdf",width=10,height=5,pointsize=12)
par(mar=c(10,5,3,1))
	plot(z$nbhd.wish[order(z$nbhd.wish)],
			xlab="",ylab="neighborhood size",
			xaxt='n',pch=19,cex=2,main="Neighborhood size",
			ylim=range(z$nbhd.wish)+c(-3,5),col=adjustcolor(z$cladecolor[order(z$nbhd.wish)],0.5))
	text(1:length(z$nbhd.wish),par("usr")[3],labels=gsub("_"," ",z$species[order(z$nbhd.wish)]),srt=50,adj=c(1.1,1.1),xpd=TRUE)
dev.off()


################################
# analyze s with one predictor at a time
#	MVN model, log(1-homozygosity)
################################

s <- log(1-z$s.wish)
db <- lapply(1:nrow(predictors),
			function(i){
				makeDB(predictors[i,],s,phyStr)
		})
fits <- lapply(1:nrow(predictors),
			function(i){
				sampling(object=mvnPhyReg,
					data=db[[i]],
					iter=nIter,
					thin=nIter/500,
					chains=1,
					control=setNames(list(15),"max_treedepth"))
			})
names(db) <- predNames
names(fits) <- predNames
out <- list("db"=db,"fits"=fits)
save(out,file="mvn_s_logpi.Robj")

pdf(file="mvn_s_logpi.pdf",width=12,height=10)
	par(mfrow=c(3,4)) ; for(i in 1:11){mvnPPS(out$db[[i]],out$fit[[i]],500,predName=names(out$fit)[i])}
dev.off()

################################
# analyze s with one predictor at a time
#	beta model, scaled s
mnS <- min(z$s.wish)
mxS <- max(z$s.wish - mnS)
s <- (z$s.wish - mnS)/mxS
s[which(s==0)] <- 1e-8
s[which(s==1)] <- 1-1e-8

db <- lapply(1:nrow(predictors),
			function(i){
				makeDB(predictors[i,],s,phyStr)
		})
fits <- lapply(1:nrow(predictors),
			function(i){
				sampling(object=betaPhyReg,
					data=db[[i]],
					iter=nIter,
					thin=nIter/500,
					chains=1,
					control=setNames(list(15),"max_treedepth"))
			})
names(db) <- predNames
names(fits) <- predNames
out <- list("db"=db,"fits"=fits)
save(out,file="beta_s_scld.Robj")

pdf(file="beta_s_scld.pdf",width=12,height=10)
	par(mfrow=c(3,4)) ; for(i in 1:11){betaPPS(out$db[[i]],out$fit[[i]],500,predName=names(out$fit)[i])}
dev.off()

load("mvn_s_logpi.Robj")
pdf(file="mvn_s_logpi.pdf",width=12,height=10)
	par(mfrow=c(3,4)) ; for(i in 1:11){mvnPPS(out$db[[i]],out$fit[[i]],500,predName=names(out$fit)[i])}
dev.off()
pdf(file="betas_mvn_s_logpi.pdf",width=12,height=10)
	postBetaPlot(out=out,predNames=names(out$fits),reorder=TRUE,cols=NULL)
dev.off()

load("beta_s_scld.Robj")
pdf(file="beta_s_scld.pdf",width=12,height=10)
	par(mfrow=c(3,4)) ; for(i in 1:11){betaPPS(out$db[[i]],out$fit[[i]],500,predName=names(out$fit)[i])}
dev.off()
pdf(file="betas_beta_s_scld.pdf",width=12,height=10)
	postBetaPlot(out=out,predNames=names(out$fits),reorder=TRUE,cols=NULL)
dev.off()


phylooViz <- function(db,CNsamples,tree,xlim){
	library(dplyr)
	library(ggplot2)
	library(ggridges)
	library(ape)
	tree <- ladderize(tree,right=FALSE)
	tipOrder <- getTipOrder(tree)
	plot(tree)
	load("db_MVN_unscaled_s.Robj")
	db <- db[[4]]
	load("s_mvn_body_size_loo_CNsamples.Robj")
	CNsamples <- looCNsamples
	
	spNames <- row.names(db$relMat)
	cnFits <- unlist(
				lapply(1:db$N,
					function(n){
						2*(abs(0.5-ecdf(CNsamples[[n]])(db$Y[n])))
				}))
	df <- Reduce("rbind",
			lapply(1:length(CNsamples[[1]]),
				function(i){data.frame("species"=spNames,
									   "cns"=unlist(lapply(CNsamples,"[[",i)),
									   "cnf"=cnFits,
									   "s"=db$Y)}))
	df <- df %>% arrange(species) %>% group_by(species) %>% mutate(o = cur_group_id()) %>% arrange(tipOrder)

	tp <- ggtree(tree) + xlim(c(0,1.5e3)) + geom_tiplab(size=2)
	df %>% ggplot() + 
  			geom_density_ridges(aes(x=cns,y=species,fill=cnf), scale=2,rel_min_height=0.05) +
  			geom_segment(aes(x=s, xend=s, y=o, yend=o+1.5), colour="red") + 
  			xlim(c(0.973,1.001)) +
  			theme_bw() +
  			theme(panel.grid = element_blank()) +
  			labs(y = "Species")
	tp + 
	geom_facet(panel="CV fit") + 

	
	p1 <- ggtree(tree)
	p2 <- facet_plot(p1,panel="CV fit",data=df,geom_joy,mapping=aes(x=))

	ggtree(tree) + 
		xlim(c(0,2e3)) + 
		geom_tiplab(offset=500,hjust=0,align=FALSE,size=1.8) + 
		geom_facet(panel="CV fit",data=df,geom=geom_col,
					aes(x=df$cnFits),orientation="y",width=0.6)
	
}








