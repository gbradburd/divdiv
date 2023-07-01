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
#mvnPhyReg <- stan_model(model_code=mvnPhyReg)
betaPhyReg <- stan_model(model_code=betaPhyReg)
expPhyReg <- stan_model(model_code=expPhyReg)

################################
# get phylo structure,
# load master dataframe,
#	create dataBlocks
################################
z <- read.csv("../../data/master_df.csv",header=TRUE,stringsAsFactors=FALSE)
z$species[which(z$species=="Seriola_lalandi_dorsalis")] <- "Seriola_dorsalis"
if(any(is.na(z$s))){
	z <- z[!which(is.na(z$s)),]	
}

load("../../data/phylo/divdiv_phy_from_timetreebeta5.Robj")
phy$tip.label[which(phy$tip.label=="Exaiptasia pallida")] <- "Exaiptasia diaphana"
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
#	beta model, unscaled s
################################

db <- lapply(1:nrow(predictors),
			function(i){
				makeDB(X=predictors[i,],Y=z$s,phyStr=phyStr)
		})
db <- setNames(db,predNames)
save(db,file="beta_s_DBs.Robj")

for(i in 1:length(predNames)){
	phyloo(db=db[[i]],
		   mod=betaPhyReg,
		   sampleNames=row.names(db[[i]]$relMat),
		   nIter=nIter,
		   nNodes=nNodes,
		   prefix=paste0("s_beta_",predNames[i]))
	load(paste0("s_beta_",predNames[i],"_loo.Robj"))
	looCNsamples <- lapply(1:db[[i]]$N,
						function(n){
							looCN(nCNsamples=nCNsamples,
						  		  sampleNames=row.names(db[[i]]$relMat),
						  		  masterDB=db[[i]],
						  		  looRep=loo$looRep[[n]],
						  		  link=logit,
						  		  invLink=invLogit)
					})
	save(looCNsamples,file=paste0("s_beta_",predNames[i],"_loo_CNsamples.Robj"))
	pdf(file=paste0("phylooFit_s_",predNames[i],".pdf"),width=12,height=9)
		phylooViz(db=db[[i]],CNsamples=looCNsamples,tree=phy,xlim=c(0.95,1.01))
	dev.off()
}


load("beta_s_DBs.Robj")
for(i in 1:length(predNames)){
	load(paste0("s_beta_",predNames[i],"_loo_CNsamples.Robj"))
	pdf(file=paste0("phylooFit_s_beta_",predNames[i],".pdf"),width=12,height=9)
		phylooViz(db=db[[i]],CNsamples=looCNsamples,tree=phy,xlim=c(0.95,1.01),valRange=c(0,0.01))
	dev.off()
}

################################
# analyze Nbhd with one predictor at a time
#	exp model
################################

db <- lapply(1:nrow(predictors),
			function(i){
				makeDB(predictors[i,],z$nbhd,phyStr)
		})
db <- setNames(db,predNames)
save(db,file="exp_nbhd_DBs.Robj")

for(i in 1:length(predNames)){
	phyloo(db=db[[i]],
		   mod=expPhyReg,
		   sampleNames=row.names(db[[i]]$relMat),
		   nIter=nIter,
		   nNodes=nNodes,
		   prefix=paste0("nbhd_exp_",predNames[i]))
	load(paste0("nbhd_exp_",predNames[i],"_loo.Robj"))
	looCNsamples <- lapply(1:db[[i]]$N,
						function(n){
							looCN(nCNsamples=nCNsamples,
						  		  sampleNames=row.names(db[[i]]$relMat),
						  		  masterDB=db[[i]],
						  		  looRep=loo$looRep[[n]],
						  		  link=log,
						  		  invLink=exp)
					})
	save(looCNsamples,file=paste0("nbhd_exp_",predNames[i],"_loo_CNsamples.Robj"))
}

load("exp_nbhd_DBs.Robj")
for(i in 1:length(predNames)){
	load(paste0("nbhd_exp_",predNames[i],"_loo_CNsamples.Robj"))
	pdf(file=paste0("phylooFit_nbhd_exp_",predNames[i],".pdf"),width=12,height=9)
		phylooViz(db=db[[i]],CNsamples=looCNsamples,tree=phy,xlim=c(-10,500)+c(min(db[[i]]$Y),max(db[[i]]$Y)),valRange=c(0,500))
	dev.off()
}

#GRAVEYARD
if(FALSE){
	# z <- z[-which(grepl("PRJNA528403",z$link)),]
	# z <- z[-which(grepl("PRJNA392526",z$link)),]
	# z <- z[-which(grepl("PRJNA314732",z$link)),]

}