################################################################
################################################################
# Running DivDiv phylogenetic LOO analyses
################################################################
################################################################


################################
# load libraries and define 
#	useful functions
################################
library(rstan)
library(ape)

source("../divdiv_analysis_functions.R")

################################
# compile rstan models
################################
source("../trait_mod_stan_blocks.R")
mvnPhyReg <- stan_model(model_code=mvnPhyReg)
betaPhyReg <- stan_model(model_code=betaPhyReg)
expPhyReg <- stan_model(model_code=expPhyReg)

################################
# get phylo structure,
# load master dataframe,
#	create dataBlocks
################################
z <- read.csv("../../data/master_df.csv",header=TRUE,stringsAsFactors=FALSE)
z <- z[-which(z$species=="Pocillopora_damicornis"),]
z <- z[-which(is.na(z$s.wish)),]

load("../../data/phylo/divdiv_phy_from_timetreebeta5.Robj")
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

predNames <- c("mean_species_latitude","number_of_ecoregions",
			   "range_extent","body_size",
			   "egg_size","generational_structure",
			   "philopatry","spawning_mode",
			   "larval_feeding","pelagic_larval_duration",
			   "planktonicity")

################################
# set global LOO params
################################

nIter <- 3e3
nNodes <- 10
nCNsamples <- 500

################################
# analyze s with one predictor at a time
#	MVN model, unscaled s
################################

db <- lapply(1:nrow(predictors),
			function(i){
				makeDB(predictors[i,],z$s.wish,phyStr)
		})

for(i in 1:length(predNames)){
	phyloo(db=db[[i]],
		   mod=mvnPhyReg,
		   sampleNames=row.names(db[[i]]$relMat),
		   nIter=nIter,
		   nNodes=nNodes,
		   prefix=paste0("s_mvn_",predNames[i]))
	load(paste0("s_mvn_",predNames[i],"_loo.Robj"))
	looCNsamples <- lapply(1:db[[i]]$N,
						function(n){
							looCN(nCNsamples=nCNsamples,
						  		  sampleNames=row.names(db[[i]]$relMat),
						  		  masterDB=db[[i]],
						  		  looRep=loo$looRep[[n]],
						  		  link=identity,
						  		  invLink=identity)
					})
	save(looCNsamples,file=paste0("s_mvn_",predNames[i],"_loo_CNsamples.Robj"))
}


if(FALSE){

################################
# analyze s with one predictor at a time
#	beta model, unscaled s
################################

db <- lapply(1:nrow(predictors),
			function(i){
				makeDB(predictors[i,],z$s.wish,phyStr)
		})


################################
# analyze Nbhd with one predictor at a time
#	MVN model
################################

db <- lapply(1:nrow(predictors),
			function(i){
				makeDB(predictors[i,],z$nbhd.wish,phyStr)
		})


################################
# analyze Nbhd with one predictor at a time
#	exp model
################################

db <- lapply(1:nrow(predictors),
			function(i){
				makeDB(predictors[i,],z$nbhd.wish,phyStr)
		})



test <- looCN(nCNsamples=1e3,sampleNames=row.names(loo$db$relMat),masterDB=loo$db,looRep=loo$looReps[[2]],link=logit,invLink=invLogit)

test <- lapply(loo$looReps,function(z){
	looCN(nCNsamples=1e3,sampleNames=row.names(loo$db$relMat),masterDB=loo$db,looRep=z,link=logit,invLink=invLogit)
})

par(mfrow=c(6,6)) ; invisible(lapply(1:db[[4]]$N,function(i){hist(looCNsamples[[i]]) ; abline(v=db[[4]]$Y[i],col="red",lwd=2.5)}))
}