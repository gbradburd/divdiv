################################################################
################################################################
# Running DivDiv biological model analyses
################################################################
################################################################

set.seed(123)

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
library(vioplot)
library(corrplot)

source("../divdiv_analysis_functions.R")

################################
# get phylo structure,
# load master dataframe,
#	create dataBlocks
################################

z <- read.csv("../../data/master_df.csv",header=TRUE,stringsAsFactors=FALSE)

load("../../data/phylo/divdiv_phy_from_timetreebeta5.Robj")
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
# test for phylogenetic signal
################################

# blomberg's k
blom_k <- phytools::phylosig(sampPhy,z$div,test=TRUE)
pagels_lambda <- phytools::phylosig(sampPhy,z$div,method="lambda",test=TRUE)
	#plot(blom_k)

load("../phy_cgram.Robj")
#plot(cgram$sampPhy_cgram)

################################
# analyze diversity with one biological predictor 
#	and all the "nuisance" parameters
################################

load("../partB_outs.Robj")
#paper figure
pps <- getBetaPPS(db=outs[[3]]$db,fit=outs[[3]]$fit,nPPS=1e3,multiPred=TRUE)

pdf(file="predictor_effect_sizes.pdf",width=8,height=7)
	par(cex.axis=1.5,cex.lab=1.5,mar=c(5,1,1,1))
	postBetaPlot(outs=outs[-c(4,14)],
				 predNames=c("Latitude","Ecoregions","Range extent",	#"ecoregions/range_size",
			 				 "Body size","Egg size","Iteroparity",
							 "Philopatry","Spawning mode","Planktotrophy",
							 "Pelagic larval duration","Planktonicity","Benthicity"),
				 reorder=TRUE,
				 cols=NULL,
				 stdize=TRUE,
				 multiPred=TRUE,
				 qnt=1)
dev.off()

pdf(file="range.pdf",width=7,height=7)
	par(cex.axis=1.5,cex.lab=1.5,mar=c(5,5,2,1))
	plotBetaPPS(db=outs[[3]]$db,ppsOut=pps,sampCols=sampCols,
				xlab="Range Extent (km)",ylab="Genetic Diversity",trp=0.7,xaxt='n')
				max(z[["max.95.sea.gbif"]])
		axis(side=1,at=seq(0,1,by=5e3/max(z[["max.95.sea.gbif"]])),
			labels=round(seq(0,1,by=5e3/max(z[["max.95.sea.gbif"]]))*max(z[["max.95.sea.gbif"]])))
dev.off()

pdf(file="planktonic.pdf",width=7,height=7)
	par(cex.axis=1.5,cex.lab=1.5,mar=c(5,5,2,1))
	discreteViolPlot(z,2,"isPlanktonic_atanypoint",xAxLabs=c("Non-planktonic","Planktonic"),logY=TRUE)
dev.off()

pdf(file="range_planktonic.pdf",width=14,height=7)
	par(cex.axis=1.5,cex.lab=1.5,mar=c(5,5,2,1),mfrow=c(1,2))
	plotBetaPPS(db=outs[[3]]$db,ppsOut=pps,sampCols=sampCols,
				xlab="Range Extent (km)",ylab="Genetic Diversity",trp=0.7,xaxt='n')
				max(z[["max.95.sea.gbif"]])
		axis(side=1,at=seq(0,1,by=5e3/max(z[["max.95.sea.gbif"]])),
			labels=round(seq(0,1,by=5e3/max(z[["max.95.sea.gbif"]]))*max(z[["max.95.sea.gbif"]])))
	discreteViolPlot(z,2,"isPlanktonic_atanypoint",xAxLabs=c("Non-planktonic","Planktonic"),logY=TRUE)
dev.off()

xx <- z[bioPreds][,-c(4,14)]
xxnames <- c("latitude","ecoregions","range extent",	#"ecoregions/range_size",
			 "body size","egg size","iteroparity",
			 "philopatry","spawning mode","planktotrophy",
			 "PLD","planktonicity","benthicity") #,"dispersive adults"
reord <- c(6:9,11:12,1:2,10,4:5,3)
xx <- xx[,reord]
xxnames <- xxnames[reord]
pdf(file="all_predictors.pdf",width=14,height=10)
	phyViz(db=outs[[3]]$db,fit=outs[[3]]$fit,
			XX=xx,
			predNames=xxnames,
		   tree=sampPhy,xlim=c(0,0.035),tipcols=z$cladecolor,
		   valRange=NULL,adj=0.5,logX=FALSE,rounding=0.05)
dev.off()


################################
# visualize predictor corrleations
################################

M <- cor(xx,use="pairwise.complete.obs",method="kendall")
row.names(M) <- xxnames
colnames(M) <- xxnames

pdf(file="predictor_corrs.pdf",width=6,height=6)
	par(cex.axis=1.5,cex.lab=1.5)
	corrplot(M,method="ellipse",diag=TRUE,type="lower",tl.col="black")
dev.off()




################################
# report statistics for paper
################################

# diversity statistics
range(z$div)
mean(z$div)
#phylo corrected mean
I <- matrix(rep(1,nrow(z)),nrow(z),1)
C_inv <- solve(phyStr)
1/(t(I) %*% C_inv %*% I) * t(I) %*% C_inv %*% z$div
mean(solve(chol(phyStr)) %*% z$div)


# effect of range extent
b <- extract(outs[[3]]$fit,"beta[1]",permute=FALSE,inc_warmup=FALSE)
mean(b)
quantile(b,c(0.025,0.975))

# effect of number of ecoregions
b <- extract(outs[[2]]$fit,"beta[1]",permute=FALSE,inc_warmup=FALSE)
mean(b)
quantile(b,c(0.025,0.975))

# effect of planktonicity
b <- extract(outs[[12]]$fit,"beta[1]",permute=FALSE,inc_warmup=FALSE)
mean(b)
quantile(b,c(0.025,0.975))

# effect of number of ecoregions PER range extent
b <- extract(outs[[4]]$fit,"beta[1]",permute=FALSE,inc_warmup=FALSE)
mean(b)
quantile(b,c(0.025,0.975))

# correlation between ecoregions & range extent
cor(z$n_ECOREGIONS.all,z$max.95.sea.gbif.nrm)