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

load("../data/phylo/divdiv_phy_from_timetreebeta5.Robj")
sampPhy <- phy
sampPhy <- ape::keep.tip(sampPhy,gsub("_"," ",z$species))
phyStr <- ape::vcv(sampPhy,corr=TRUE)

predictors <- rbind(z[["meanlat.gbif"]],
					z[["n_ECOREGIONS.all"]],
					z[["maxgbif.sea_km"]],
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

################################
# run analyses
################################

# modeling collecting phase pi
s <- (z$s - min(z$s))/max((z$s - min(z$s)))
db_s <- list("N" = nrow(phyStr),
		   	 "Y" = s,
		   	 "nX"= nrow(predictors),
		   	 "X" = predictors,
		     "relMat" = phyStr)

# modeling neighborhood sizes
db_nbhd <- list("N" = nrow(phyStr),
		   	    "Y" = z$nbhd,
		   	    "nX"= nrow(predictors),
		   	    "X" = predictors,
		        "relMat" = phyStr)

fit_s <- sampling(object=mvn,
					data=db_s,
					iter=3e3,
					chains=2,
					control=setNames(list(15),"max_treedepth"))

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
	par(mfrow=mf,mar=c(3,3,3,3))
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

predNames <- c("meanlat.gbif","n_ECOREGIONS.all",
			   "maxgbif.sea_km","Body_Size",
			   "Fecundity_EggSize","Generational_Structure",
			   "ReturnToSpawningGround","Spawning_mode",
			   "Larval_feeding","PLD_point2",
			   "isPlanktonic_atanypoint")

b_s <- rstan::extract(fit_s,"beta",permute=FALSE)
b_nbhd <- rstan::extract(fit_nbhd,"beta",permute=FALSE)

pdf(file="divdiv_s.pdf",width=12,height=8)
	plotMarginals(b_s[,1,],predNames,c(3,4))
		plot(0,type='n',bty='n',yaxt='n',xaxt='n',xlab="",ylab="")
		legend(x="topright",col=c(4,1),lty=c(1,2),lwd=3,legend=c("95% CI","Zero"),cex=1.5,bty='n')
dev.off()

pdf(file="divdiv_nbhd.pdf",width=12,height=8)
	plotMarginals(b_nbhd[,1,],predNames,c(3,4))
		plot(0,type='n',bty='n',yaxt='n',xaxt='n',xlab="",ylab="")
			legend(x="topright",col=c(4,1),lty=c(1,2),lwd=3,legend=c("95% CI","Zero"),cex=1.5,bty='n')
dev.off()

plot(colMeans(b_s[,1,]))