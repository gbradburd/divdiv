################################################################
################################################################
# Running DivDiv analyses
################################################################
################################################################


################################
# load libraries and define 
#	useful functions
################################
library(rstan)
library(ape)

colFunc <- function (x, cols, nCols, valRange){
    if (is.null(valRange)) {
        valRange <- c(min(x), max(x))
    }
    cols <- (grDevices::colorRampPalette(cols))(nCols)[findInterval(x, 
        seq(valRange[1], valRange[2], length.out = nCols))]
    return(cols)
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
################################
# compile rstan models
################################
source("trait_mod_stan_block.R")
mvn <- stan_model(model_code=mvn)

################################
# get phylo structure,
# load master dataframe,
#	create dataBlocks
################################
z <- read.csv("../data/master_df.csv",header=TRUE,stringsAsFactors=FALSE)
z <- z[-which(z$species=="Pocillopora_damicornis"),]
z <- z[-which(is.na(z$s.wish)),]

load("../data/phylo/divdiv_phy_from_timetreebeta5.Robj")
sampPhy <- phy
sampPhy <- ape::keep.tip(sampPhy,gsub("_"," ",z$species))
phyStr <- ape::vcv(sampPhy,corr=TRUE)

predictors <- rbind(z[["meanlat.gbif"]],
					z[["n_ECOREGIONS.all"]],
					z[["max.95.sea.gbif"]]/max(z[["max.95.sea.gbif"]]),
					z[["Body_Size"]],
					z[["Fecundity_EggSize"]],
					z[["Generational_Structure"]],
					z[["ReturnToSpawningGround"]],
					z[["Spawning_mode"]],
					z[["Larval_feeding"]],
					z[["PLD_point2"]],
					z[["isPlanktonic_atanypoint"]])

tmp <- predictors
# fill in missing data w/ grand mean for that predictor
md <- which(is.na(predictors),arr.ind=TRUE)
for(i in 1:nrow(md)){
	predictors[md[i,1],md[i,2]] <- mean(tmp[md[i,1],],na.rm=TRUE)
}

predNames <- c("mean species latitude","number of ecoregions",
			   "range extent","body size",
			   "egg size","generational structure",
			   "philopatry","spawning mode",
			   "larval feeding","pelagic larval duration",
			   "planktonicity")

################################
# run analyses with one predictor 
#	at a time
################################

s <- (z$s.wish - min(z$s.wish))/max((z$s.wish - min(z$s.wish)))
dbss <- lapply(1:nrow(predictors),
			function(i){
				list("N" = nrow(phyStr),
				   	 "Y" = s,
				   	 "nX"= 1,
				   	 "X" = predictors[i,,drop=FALSE],
				     "relMat" = phyStr)
		})
sfits <- lapply(1:nrow(predictors),
			function(i){
				sampling(object=mvn,
					data=dbss[[i]],
					iter=3e3,
					chains=1,
					control=setNames(list(15),"max_treedepth"))
			})

unscale <- list("mn" = min(z$s.wish),"mx" = max((z$s.wish - min(z$s.wish))))
par(mfrow=c(3,4))
lapply(1:nrow(predictors),
	function(i){
		vizSinglePredictorFit(db=dbss[[i]],
							  fit=sfits[[i]],
							  xName=predNames[i],
							  yName="deep-time diversity",
							  cols=z$cladecol,unscale=unscale)
	})

dbns <- lapply(1:nrow(predictors),
			function(i){
				list("N" = nrow(phyStr),
				   	 "Y" = z$nbhd.wish,
				   	 "nX"= 1,
				   	 "X" = predictors[i,,drop=FALSE],
				     "relMat" = phyStr)
		})
nfits <- lapply(1:nrow(predictors),
			function(i){
				sampling(object=mvn,
					data=dbns[[i]],
					iter=3e3,
					chains=1,
					control=setNames(list(15),"max_treedepth"))
			})
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
