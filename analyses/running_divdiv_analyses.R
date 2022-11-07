################################################################
################################################################
# Running DivDiv biological model analyses
################################################################
################################################################


################################
# load libraries and define 
#	useful functions
################################
library(rstan)
library(ape)

source("divdiv_analysis_functions.R")

################################
# compile rstan models
################################
source("trait_mod_stan_blocks.R")
mvnPhyReg <- stan_model(model_code=mvnPhyReg)
betaPhyReg <- stan_model(model_code=betaPhyReg)
expPhyReg <- stan_model(model_code=expPhyReg)

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

if(FALSE){
	tmp <- predictors
	# fill in missing data w/ grand mean for that predictor
	md <- which(is.na(predictors),arr.ind=TRUE)
	for(i in 1:nrow(md)){
		predictors[md[i,1],md[i,2]] <- mean(tmp[md[i,1],],na.rm=TRUE)
	}
}

predNames <- c("mean species latitude","number of ecoregions",
			   "range extent","body size",
			   "egg size","generational structure",
			   "philopatry","spawning mode",
			   "larval feeding","pelagic larval duration",
			   "planktonicity")

nIter <- 5e3

################################
# analyze s with one predictor at a time
#	MVN model, unscaled s
################################

db <- lapply(1:nrow(predictors),
			function(i){
				makeDB(predictors[i,],z$s.wish,phyStr)
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
save(out,file="mvn_s_unscld.Robj")

pdf(file="mvn_s_unscld.pdf",width=12,height=10)
	par(mfrow=c(3,4)) ; for(i in 1:11){mvnPPS(out$db[[i]],out$fit[[i]],500,predName=names(out$fit)[i])}
dev.off()

################################
# analyze s with one predictor at a time
#	MVN model, scaled s
################################

s <- (z$s.wish - min(z$s.wish))/max((z$s.wish - min(z$s.wish)))
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

unscale <- list("mn" = min(z$s.wish),"mx" = max((z$s.wish - min(z$s.wish))))
names(db) <- predNames
names(fits) <- predNames
out <- list("db"=db,"fits"=fits,"unscale"=unscale)
save(out,file="mvn_s_scld.Robj")

pdf(file="mvn_s_scld.pdf",width=12,height=10)
	par(mfrow=c(3,4)) ; for(i in 1:11){mvnPPS(out$db[[i]],out$fit[[i]],500,predName=names(out$fit)[i])}
dev.off()


################################
# analyze s with one predictor at a time
#	beta model, unscaled s
################################

db <- lapply(1:nrow(predictors),
			function(i){
				makeDB(predictors[i,],z$s.wish,phyStr)
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
save(out,file="beta_s.Robj")

pdf(file="beta_s.pdf",width=12,height=10)
	par(mfrow=c(3,4)) ; for(i in 1:11){betaPPS(out$db[[i]],out$fit[[i]],500,predName=names(out$fit)[i])}
dev.off()


################################
# analyze Nbhd with one predictor at a time
#	MVN model
################################

db <- lapply(1:nrow(predictors),
			function(i){
				makeDB(predictors[i,],z$nbhd.wish,phyStr)
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
save(out,file="mvn_Nbhd.Robj")

pdf(file="mvn_Nbhd.pdf",width=12,height=10)
	par(mfrow=c(3,4)) ; for(i in 1:11){mvnPPS(out$db[[i]],out$fit[[i]],500,predName=names(out$fit)[i])}
dev.off()

################################
# analyze Nbhd with one predictor at a time
#	exp model
################################

db <- lapply(1:nrow(predictors),
			function(i){
				makeDB(predictors[i,],z$nbhd.wish,phyStr)
		})
fits <- lapply(1:nrow(predictors),
			function(i){
				sampling(object=expPhyReg,
					data=db[[i]],
					iter=nIter,
					thin=nIter/500,
					save_warmup=FALSE,
					chains=1,
					control=setNames(list(15),"max_treedepth"))
			})
names(db) <- predNames
names(fits) <- predNames
out <- list("db"=db,"fits"=fits)
save(out,file="exp_Nbhd.Robj")

pdf(file="exp_Nbhd.pdf",width=12,height=10)
	par(mfrow=c(3,4)) ; for(i in 1:11){expPPS(out$db[[i]],out$fit[[i]],500,predName=names(out$fit)[i])}
dev.off()


if(FALSE){
load("mvn_s_unscld.Robj")
pdf(file="mvn_s_unscld.pdf",width=12,height=10)
	par(mfrow=c(3,4)) ; for(i in 1:11){mvnPPS(out$db[[i]],out$fit[[i]],500,predName=names(out$fit)[i])}
dev.off()
pdf(file="betas_mvn_s_unscld.pdf",width=12,height=10)
	postBetaPlot(out=out,predNames=names(out$fits),reorder=TRUE,cols=NULL)
dev.off()


load("mvn_s_scld.Robj")
pdf(file="mvn_s_scld.pdf",width=12,height=10)
	par(mfrow=c(3,4)) ; for(i in 1:11){mvnPPS(out$db[[i]],out$fit[[i]],500,predName=names(out$fit)[i])}
dev.off()
pdf(file="betas_mvn_s_scld.pdf",width=12,height=10)
	postBetaPlot(out=out,predNames=names(out$fits),reorder=TRUE,cols=NULL)
dev.off()


load("beta_s.Robj")
pdf(file="beta_s.pdf",width=12,height=10)
	par(mfrow=c(3,4)) ; for(i in 1:11){betaPPS(out$db[[i]],out$fit[[i]],500,predName=names(out$fit)[i])}
dev.off()
pdf(file="betas_beta_s.pdf",width=12,height=10)
	postBetaPlot(out=out,predNames=names(out$fits),reorder=TRUE,cols=NULL)
dev.off()


load("mvn_Nbhd.Robj")
pdf(file="mvn_Nbhd.pdf",width=12,height=10)
	par(mfrow=c(3,4)) ; for(i in 1:11){mvnPPS(out$db[[i]],out$fit[[i]],500,predName=names(out$fit)[i])}
dev.off()
pdf(file="betas_mvn_Nbhd.pdf",width=12,height=10)
	postBetaPlot(out=out,predNames=names(out$fits),reorder=TRUE,cols=NULL)
dev.off()


load("exp_Nbhd.Robj")
pdf(file="exp_Nbhd.pdf",width=12,height=10)
	par(mfrow=c(3,4)) ; for(i in 1:11){expPPS(out$db[[i]],out$fit[[i]],500,predName=names(out$fit)[i])}
dev.off()
pdf(file="betas_exp_Nbhd.pdf",width=12,height=10)
	postBetaPlot(out=out,predNames=names(out$fits),reorder=TRUE,cols=NULL)
dev.off()
}