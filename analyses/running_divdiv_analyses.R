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
options(mc.cores=12)
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

load("../data/phylo/divdiv_phy_from_timetreebeta5.Robj")
sampPhy <- ape::keep.tip(phy,gsub("_"," ",z$species))
phyStr <- ape::vcv(sampPhy,corr=TRUE)

sp <- gsub("_"," ",z$species)
z <- z[match(row.names(phyStr),sp),]

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