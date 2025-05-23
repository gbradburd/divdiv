################################################################
################################################################
# Running DivDiv biological model analyses
################################################################
################################################################

set.seed(248)

################################
# load libraries and define 
#	useful functions
################################
library(rstan)
options(mc.cores=12)
library(ape)
library(doParallel)
library(foreach)
library(phytools)
library(phylosignal)

source("divdiv_analysis_functions.R")

#specify global code options
nNodes <- 15
nIter <- 3e4
interactive <- FALSE

################################
# compile rstan models
################################
source("trait_mod_stan_blocks.R")
if(interactive){
	betaPhyReg <- stan_model(model_code=betaPhyReg)	
}


################################
# get phylo structure,
# load master dataframe,
#	create dataBlocks
################################

z <- read.csv("../data/master_df.csv",header=TRUE,stringsAsFactors=FALSE)

load("../data/phylo/divdiv_phy_from_timetreebeta5.Robj")
sampPhy <- ape::keep.tip(phy,gsub("_"," ",z$species))
phyStr <- ape::vcv(sampPhy,corr=TRUE)


# reorder dataframe to match order of species
#	in phylogenetic variance-covariance matrix
sp <- gsub("_"," ",z$species)
z <- z[match(row.names(phyStr),sp),]

# prepare dataframe for analyses
z[["div"]] <- 1-z$s
z[["Log_BodySize"]] <- log(z[["Body_Size"]])
z[["max.95.sea.gbif.nrm"]] <- z[["max.95.sea.gbif"]]/max(z[["max.95.sea.gbif"]])
z[["eco_per_range"]] <- z[["n_ECOREGIONS.all"]]/z[["max.95.sea.gbif.nrm"]]
z[["mean_raw_read_cnt.nrm"]] <- z[["mean_raw_read_cnt"]]/1e7
z[["n.totalsnps.nrm"]] <- z[["n.totalsnps"]]/1e5

# make vector of predictor names
preds <- c("meanlat.gbif","n_ECOREGIONS.all","max.95.sea.gbif.nrm","eco_per_range",
			"Log_BodySize","Fecundity_EggSize","Generational_Structure",
			"ReturnToSpawningGround","Spawning_mode",
			"Larval_feeding","PLD_point2","isPlanktonic_atanypoint",
			"isBenthic","Large_Adult_Range",
			"ratio.sea.95","n_samples","mean_raw_read_cnt.nrm",
			"read_length","mean_locus_depth")

predNames <- c("mean species latitude","number of ecoregions","range extent","ecoregions/range_size",
			   "body size (log)","egg size","generational structure",
			   "philopatry","spawning mode",
			   "larval feeding","pelagic larval duration","planktonicity",
			   "benthic lifestyle (N/S/A)", "highly dispersive adults",
			   "ratio of range sampled : total",
			   "number of samples","mean raw read count",
			   "read length","mean locus depth")

sampCols <- z$cladecolor
names(sampCols) <- gsub("_"," ",z$species)
bioPreds <- preds[1:14]
bioPredNames <- predNames[1:14]
nuisPreds <- preds[16:19]
nuisPredNames <- predNames[16:19]


################################
# analyze diversity with one predictor at a time
################################

outs <- divdivAnalyses(z=z,X=as.list(preds),
						y="div",phyStr=phyStr,mod=betaPhyReg,
						outName="partA",nIter=nIter,parallel=TRUE,
						nNodes=nNodes,filterKeep=NULL,runNames=preds)

vizAllOuts(outs=outs,predNames=predNames,sampPhy=sampPhy,outName="partA",sampCols=sampCols)

################################
# analyze diversity with one biological predictor 
#	and all the "nuisance" parameters
################################

x <- lapply(1:length(bioPreds),function(i){c(bioPreds[i],nuisPreds)})

outs <- divdivAnalyses(z=z,X=x,
						y="div",phyStr=phyStr,mod=betaPhyReg,
						outName="partB",nIter=nIter,parallel=TRUE,
						nNodes=nNodes,filterKeep=NULL,runNames=bioPreds)

vizAllOuts(outs=outs,predNames=predNames,sampPhy=sampPhy,outName="partB",multiPred=TRUE,sampCols=sampCols)

################################
# BOUTIQUE ANALYSES
#	e.g., with specific combinations 
#	of predictors, or in 
#	specific taxonomic groups
################################

################
# diversity ~ planktotrophic + !isCoral
################

x <- list(c("Larval_feeding",nuisPreds))

outs <- divdivAnalyses(z=z,X=x,
						y="div",phyStr=phyStr,mod=betaPhyReg,
						outName="partC",nIter=nIter,parallel=TRUE,
						nNodes=nNodes,filterKeep=which(z$isCoral!=1),runNames="Larval_feeding")

vizAllOuts(outs=outs,predNames="Larval_feeding",sampPhy=sampPhy,outName="partC",multiPred=TRUE,sampCols=sampCols)


################
# diversity ~ pld + isBenthic
################

x <- list(c("PLD_point2","isBenthic"))

outs <- divdivAnalyses(z=z,X=x,
						y="div",phyStr=phyStr,mod=betaPhyReg,
						outName="partD",nIter=nIter,parallel=TRUE,
						nNodes=nNodes,filterKeep=NULL,runNames=c("PLD_point_and_isBenthic"))

vizAllOuts(outs=outs,predNames="PLD_point2",sampPhy=sampPhy,outName="partD",multiPred=TRUE,sampCols=sampCols)

################
# analyze diversity with one biological predictor 
#	and all the "nuisance" parameters
#	JUST IN ACTINOPTERYGII
################

bonyFishMRCA <- ape::getMRCA(sampPhy,c(which(sampPhy$tip.label=="Engraulis encrasicolus"),which(sampPhy$tip.label=="Sebastes diaconus")))
isBonyFish <- phangorn::Descendants(sampPhy,bonyFishMRCA,type="tips")[[1]]

x <- lapply(1:length(bioPreds),function(i){c(bioPreds[i],nuisPreds)})

outs <- divdivAnalyses(z=z,X=x,
						y="div",phyStr=phyStr,mod=betaPhyReg,
						outName="partE",nIter=nIter,parallel=TRUE,
						nNodes=nNodes,filterKeep=isBonyFish,runNames=bioPreds)

vizAllOuts(outs=outs,predNames=predNames,sampPhy=sampPhy,outName="partE",multiPred=TRUE,sampCols=sampCols)

################
# analyze diversity with one biological predictor 
#	and all the "nuisance" parameters
#	dropping all species that are secondarily marine
################

mammalMRCA <- ape::getMRCA(sampPhy,c(which(sampPhy$tip.label=="Halichoerus grypus atlantica"),which(sampPhy$tip.label=="Phocoena sinus")))
isMammal <- phangorn::Descendants(sampPhy,mammalMRCA,type="tips")[[1]]
birdMRCA <- ape::getMRCA(sampPhy,c(which(sampPhy$tip.label=="Aythya marila"),which(sampPhy$tip.label=="Caretta caretta")))
isBird <- phangorn::Descendants(sampPhy,birdMRCA,type="tips")[[1]]
isPlant <- c(which(z$species=="Rhizophora_mangle"),
			 which(z$species=="Laguncularia_racemosa"),
			 which(z$species=="Sargassum_muticum"))
primarilyMarine <- (1:nrow(z))[-c(isMammal,isBird,isPlant)]


x <- lapply(1:length(bioPreds),function(i){c(bioPreds[i],nuisPreds)})

outs <- divdivAnalyses(z=z,X=x,
						y="div",phyStr=phyStr,mod=betaPhyReg,
						outName="partF",nIter=nIter,parallel=TRUE,
						nNodes=nNodes,filterKeep=primarilyMarine,runNames=bioPreds)

vizAllOuts(outs=outs,predNames=predNames,sampPhy=sampPhy,outName="partF",multiPred=TRUE,sampCols=sampCols)

################
# analyze diversity with one biological predictor 
#	and all the "nuisance" parameters
#	dropping just dolphins
################

toKeep <- (1:nrow(z))[-grep("Lagenorhynchus",z$species)]

x <- lapply(1:length(bioPreds),function(i){c(bioPreds[i],nuisPreds)})

outs <- divdivAnalyses(z=z,X=x,
						y="div",phyStr=phyStr,mod=betaPhyReg,
						outName="partG",nIter=nIter,parallel=TRUE,
						nNodes=nNodes,filterKeep=toKeep,runNames=bioPreds)

vizAllOuts(outs=outs,predNames=predNames,sampPhy=sampPhy,outName="partG",multiPred=TRUE,sampCols=sampCols)

################
# analyze diversity with all predictors simultaneously
################
if(FALSE){
predSampCov <- cov(t(z[,bioPreds]),use="pairwise.complete.obs")
predPCs <- prcomp(predSampCov)
predPC

x <- list(c(bioPreds,nuisPreds))

outs <- divdivAnalyses(z=z,X=x,
						y="div",phyStr=phyStr,mod=betaPhyReg,
						outName="partH",nIter=nIter,parallel=TRUE,
						nNodes=nNodes,filterKeep=NULL,runNames="allPreds")

postBetaPlot_allPreds(outs[[1]],x[[1]])

tmp <- lapply(1:18,function(i){rstan::extract(outs[[1]]$fit,par=sprintf("beta[%s]",i),permute=FALSE,inc_warmup=TRUE)})

par(mfrow=c(4,5)) ; lapply(1:18,function(i){hist(tmp[[i]],main=x[[1]][i]) ; abline(v=0,col="red")})
}
