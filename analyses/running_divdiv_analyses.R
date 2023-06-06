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
library(doParallel)
library(foreach)

source("divdiv_analysis_functions.R")

################################
# compile rstan models
################################
source("trait_mod_stan_blocks.R")
betaPhyReg <- stan_model(model_code=betaPhyReg)
expPhyReg <- stan_model(model_code=expPhyReg)

################################
# get phylo structure,
# load master dataframe,
#	create dataBlocks
################################
z <- read.csv("../data/master_df.csv",header=TRUE,stringsAsFactors=FALSE)
z <- z[-which(z$species=="Pocillopora_damicornis"),]
if(any(is.na(z$s))){
	z <- z[!which(is.na(z$s)),]	
}

#dupes <- names(which(table(z$species) > 1))

# Gadus morhua
#	sampling of PRJNA521889 is w/in range of PRJNA528403
#	model fit looks equally reasonable
#	keep PRJNA528403

z <- z[-which(grepl("PRJNA528403",z$link)),]

# Sebastiscus marmoratus
#	more samples in PRJNA359404 than PRJNA392526
#	model fit isn't amazing for either
#	keep PRJNA359404

z <- z[-which(grepl("PRJNA392526",z$link)),]

# Lateolabrax maculatus
# more samples from more locations in PRJNA356786 than PRJNA314732
# model fit looks fine for PRJNA356786
# keep PRJNA356786

z <- z[-which(grepl("PRJNA314732",z$link)),]


load("../data/phylo/divdiv_phy_from_timetreebeta5.Robj")
sampPhy <- ape::keep.tip(phy,gsub("_"," ",z$species))
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
#	beta model, unscaled s
################################

db <- lapply(1:nrow(predictors),
			function(i){
				makeDB(predictors[i,],z$s,phyStr)
		})
names(db) <- predNames

cl <- parallel::makeCluster(11)
doParallel::registerDoParallel(cl)

fits <- foreach::foreach(i = 1:nrow(predictors)) %dopar% {
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
save(out,file="beta_s.Robj")

pdf(file="beta_s.pdf",width=12,height=10)
	par(mfrow=c(3,4)) ; for(i in 1:11){betaPPS(out$db[[i]],out$fit[[i]],500,predName=names(out$fit)[i])}
dev.off()

pdf(file="betas_beta_s.pdf",width=12,height=10)
	postBetaPlot(out=out,predNames=names(out$fits),reorder=TRUE,cols=NULL)
dev.off()

################################
# analyze Nbhd with one predictor at a time
#	exp model
################################

db <- lapply(1:nrow(predictors),
			function(i){
				makeDB(predictors[i,],z$nbhd,phyStr)
		})
names(db) <- predNames

cl <- parallel::makeCluster(11)
doParallel::registerDoParallel(cl)

fits <- foreach::foreach(i = 1:nrow(predictors)) %dopar% {
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
save(out,file="exp_Nbhd.Robj")

pdf(file="exp_Nbhd.pdf",width=12,height=10)
	par(mfrow=c(3,4)) ; for(i in 1:11){expPPS(out$db[[i]],out$fit[[i]],500,predName=names(out$fit)[i])}
dev.off()

pdf(file="betas_exp_Nbhd.pdf",width=12,height=10)
	postBetaPlot(out=out,predNames=names(out$fits),reorder=TRUE,cols=NULL)
dev.off()






