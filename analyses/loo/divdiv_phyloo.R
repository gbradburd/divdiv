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
betaLoo <- stan_model(model_code=betaLoo)
#expPhyReg <- stan_model(model_code=expPhyReg)

################################
# get phylo structure,
# load master dataframe,
#	create dataBlocks
################################
z <- read.csv("../../data/master_df.csv",header=TRUE,stringsAsFactors=FALSE)

load("../../data/phylo/divdiv_phy_from_timetreebeta5.Robj")
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

bioPreds <- predictors[1:11,]
bioPredNames <- predNames[1:11]
nuisPreds <- predictors[13:16,]
nuisPredNames <- predNames[13:16]

################################
# set global LOO params
################################

nIter <- 3e3
nNodes <- 12

################################
# analyze diversity with one biological predictor 
#	and all the "nuisance" parameters
#	beta model, unscaled diversity
################################

dbs <- lapply(1:nrow(bioPreds),
			function(i){
				makeMultiDB(list(predictors[i,],
								 nuisPreds[1,],
								 nuisPreds[2,],
								 nuisPreds[3,],
								 nuisPreds[4,]),
							 Y=1-z$s,
							 phyStr=phyStr)
		})

for(i in 1:length(dbs)){
	phyloo(db=dbs[[i]],
		   mod=betaLoo,
		   sampleNames=row.names(dbs[[i]]$relMat),
		   nIter=nIter,
		   nNodes=nNodes,
		   prefix=paste0("div_",predNames[i]))
}


pdf(file="div_phyloo_scatter.pdf",width=8,height=8)
	for(i in 1:length(dbs)){
		load(sprintf("div_%s_loo.Robj",bioPredNames[i]))
		phyloViz_scatter(loo=loo,predName=bioPredNames[i],
							valRange=NULL)
	}
dev.off()

pdf(file="div_phyloo.pdf",width=12,height=12)
	for(i in 1:nrow(bioPreds)){
		load(sprintf("div_%s_loo.Robj",bioPredNames[i]))
		phylooViz(loo=loo,predName=bioPredNames[i],tree=sampPhy,
					xlim=c(0,0.03),adj=1,qnt=0.99,
					logX=FALSE)
	}
dev.off()

pdf(file="div_phyloo_log.pdf",width=12,height=12)
	for(i in 1:nrow(bioPreds)){
		load(sprintf("div_%s_loo.Robj",bioPredNames[i]))
		phylooViz(loo=loo,predName=bioPredNames[i],tree=sampPhy,
					xlim=c(0,0.03),qnt=0.99,adj=1,
					logX=TRUE)
	}
dev.off()

## GRAVEYARD
if(FALSE){
	

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
}