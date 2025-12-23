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
zCols <- read.csv("../../data/master_tax_color_key.csv",header=TRUE,stringsAsFactors=FALSE)
zCols <- zCols[-which(!zCols$species %in% gsub("_"," ",z$species)),]
zCols <- zCols[match(gsub("_"," ",z$species),zCols$species),]
z$newTaxCol <- zCols$cladecolor[match(gsub("_"," ",z$species),zCols$species)]

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
			   "philopatry","spawning",
			   "larval feeding","pelagic larval duration","planktonicity",
			   "benthic lifestyle (N/S/A)", "highly dispersive adults",
			   "ratio of range sampled : total",
			   "number of samples","mean raw read count",
			   "read length","mean locus depth")

sampCols <- z$newTaxCol	#z$cladecolor
names(sampCols) <- gsub("_"," ",z$species)
bioPreds <- preds[1:14]
bioPredNames <- predNames[1:14]
nuisPreds <- preds[16:19]
nuisPredNames <- predNames[16:19]

################################
# visualizing diversity
################################

pdf(file="suppmat_diversity.pdf",width=8.5,height=11)
	zsort <- order(z$div)
	spNames <- gsub("_"," ",z$species)
	s1 <- seq(1,length(spNames),by=2)
	s2 <- seq(2,length(spNames),by=2)
	plot(x=log(z$div)[zsort],y=1:nrow(z),
		yaxt='n',ylab="",xlab="",xaxt='n',
		col=z$newTaxCol[zsort],pch=19,
		xlim=log(c(1.5e-4,7e-2)))
		axis(side=1,log(c(seq(3e-4,3e-3,length.out=10),
			seq(3e-3,3e-2,length.out=10))),labels=FALSE)
		axis(side=1,at=log(3*c(1e-4,1e-3,1e-2,1e-1)),
				labels=format(3*c(1e-4,1e-3,1e-2,1e-1),scientific=FALSE))
		axis(side=3,log(c(seq(3e-4,3e-3,length.out=10),
			seq(3e-3,3e-2,length.out=10))),labels=FALSE)
		axis(side=3,at=log(3*c(1e-4,1e-3,1e-2,1e-1)),
				labels=format(3*c(1e-4,1e-3,1e-2,1e-1),scientific=FALSE))
		text(x=log(z$div)[zsort][s1],xpd=TRUE,y=s1,labels=spNames[zsort][s1],pos=2,col=z$newTaxCol[zsort][s1],cex=0.7)
		text(x=log(z$div)[zsort][s2],xpd=TRUE,y=s2-0.2,labels=spNames[zsort][s2],pos=4,col=z$newTaxCol[zsort][s2],cex=0.7)
		addPhylopics(unique(z$newTaxCol)[c(8,10,9,7,6,4,3,2,1,5)],
				ysize=3.1,tipcols=z$newTaxCol,xleft=log(1.9e-4),ytop=90,dshift=4.5,legendshift=0.3,txt.cex=0.9,pos=4)
dev.off()
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
							 "Philopatry","Spawning","Planktotrophy",
							 "Pelagic larval duration","Planktonicity","Benthicity"),
				 reorder=TRUE,
				 cols=NULL,
				 stdize=TRUE,
				 multiPred=TRUE,
				 qnt=1)
	arrows(x0=-1.4,x1=-1.9,y0=-1,y1=-1,lwd=4,xpd=TRUE,col="coral3",length=0.1)
		text(x=-1.6,y=-1.6,xpd=TRUE,labels="decreasing\ndiversity",col="coral3")
	arrows(x0=1,x1=1.5,y0=-1,y1=-1,lwd=4,xpd=TRUE,col="deepskyblue4",length=0.1)
		text(x=1.2,y=-1.6,xpd=TRUE,labels="increasing\ndiversity",col="deepskyblue4")
dev.off()


expFromMod <- simFromMod(outs[[3]],nXseq=100,nSampIters=1e5)
yrange <- range(c(log(outs[[3]]$db$Y),
				log(apply(expFromMod$y,1,quantile,probs=c(0.025,0.975)))))
pdf(file="range.pdf",width=7,height=7)
	par(cex.axis=1.5,cex.lab=1.5,mar=c(5,5,2,1))
	plot(outs[[3]]$db$X[1,],log(outs[[3]]$db$Y),type='n',
			xlab="Range Extent (km)",ylab="Genetic Diversity",xaxt='n',yaxt='n',
			ylim=yrange)
		axis(side=1,at=seq(0,1,by=5e3/max(z[["max.95.sea.gbif"]])),
			labels=round(seq(0,1,by=5e3/max(z[["max.95.sea.gbif"]]))*max(z[["max.95.sea.gbif"]])))
		axis(side=2,at=log(3*c(1e-4,1e-3,1e-2,1e-1)),
			labels=format(3*c(1e-4,1e-3,1e-2,1e-1),scientific=FALSE))
		axis(side=2,log(c(seq(3e-4,3e-3,length.out=10),
						  seq(3e-3,3e-2,length.out=10))),labels=FALSE)
		addExpPredPolygon(expFromMod,qnt=c(0.025,0.975),polyCol=adjustcolor("black",0.1),log=TRUE)
		addExpPredPolygon(expFromMod,qnt=c(0.10,0.90),polyCol=adjustcolor("black",0.2),log=TRUE)
		addExpPredPolygon(expFromMod,qnt=c(0.25,0.75),polyCol=adjustcolor("black",0.3),log=TRUE)
		points(outs[[3]]$db$X[1,],log(outs[[3]]$db$Y),col="black",bg=adjustcolor(sampCols,0.8),pch=21,cex=1.7)
		addExpPredLine(expFromMod,col="black",lwd=3,log=TRUE)
dev.off()

pdf(file="range_pps.pdf",width=7,height=7)
	par(cex.axis=1.5,cex.lab=1.5,mar=c(5,5,2,1))
	plotBetaPPS(db=outs[[3]]$db,ppsOut=pps,sampCols=sampCols,
				xlab="Range Extent (km)",ylab="Genetic Diversity",trp=0.7,xaxt='n')
				max(z[["max.95.sea.gbif"]])
		axis(side=1,at=seq(0,1,by=5e3/max(z[["max.95.sea.gbif"]])),
			labels=round(seq(0,1,by=5e3/max(z[["max.95.sea.gbif"]]))*max(z[["max.95.sea.gbif"]])))
dev.off()

pdf(file="planktonic.pdf",width=7,height=7)
	par(cex.axis=1.5,cex.lab=1.5,mar=c(5,5,2,1))
	discreteViolPlot(z,2,"isPlanktonic_atanypoint",xAxLabs=c("Non-planktonic","Planktonic"),sampCols=sampCols,
					logY=TRUE,pchMed=18,lineCol="gray98",yRange=NULL)
		axis(side=2,at=log(3*c(1e-4,1e-3,1e-2,1e-1)),
			labels=format(3*c(1e-4,1e-3,1e-2,1e-1),scientific=FALSE))
		axis(side=2,log(c(seq(3e-4,3e-3,length.out=10),
						  seq(3e-3,3e-2,length.out=10))),labels=FALSE)		
dev.off()

pdf(file="range_planktonic.pdf",width=14,height=7)
	par(cex.axis=1.5,cex.lab=1.5,mar=c(5,5,2,1),mfrow=c(1,2))
	plot(outs[[3]]$db$X[1,],log(outs[[3]]$db$Y),type='n',
			xlab="Range Extent (km)",ylab="Genetic Diversity",xaxt='n',yaxt='n',
			ylim=yrange)
		axis(side=1,at=seq(0,1,by=5e3/max(z[["max.95.sea.gbif"]])),
			labels=round(seq(0,1,by=5e3/max(z[["max.95.sea.gbif"]]))*max(z[["max.95.sea.gbif"]])))
		axis(side=2,at=log(3*c(1e-4,1e-3,1e-2,1e-1)),
			labels=format(3*c(1e-4,1e-3,1e-2,1e-1),scientific=FALSE))
		axis(side=2,log(c(seq(3e-4,3e-3,length.out=10),
						  seq(3e-3,3e-2,length.out=10))),labels=FALSE)
		addExpPredPolygon(expFromMod,qnt=c(0.025,0.975),polyCol=adjustcolor("black",0.1),log=TRUE)
		addExpPredPolygon(expFromMod,qnt=c(0.10,0.90),polyCol=adjustcolor("black",0.2),log=TRUE)
		addExpPredPolygon(expFromMod,qnt=c(0.25,0.75),polyCol=adjustcolor("black",0.3),log=TRUE)
		points(outs[[3]]$db$X[1,],log(outs[[3]]$db$Y),col="black",bg=adjustcolor(sampCols,0.8),pch=21,cex=1.7)
		addExpPredLine(expFromMod,col="black",lwd=3,log=TRUE)
	discreteViolPlot(z,2,"isPlanktonic_atanypoint",
						xAxLabs=c("Non-planktonic","Planktonic"),sampCols=sampCols,logY=TRUE,
						pchMed=18,lineCol="gray98",yRange=yrange)
		axis(side=2,at=log(3*c(1e-4,1e-3,1e-2,1e-1)),
			labels=format(3*c(1e-4,1e-3,1e-2,1e-1),scientific=FALSE))
		axis(side=2,log(c(seq(3e-4,3e-3,length.out=10),
						  seq(3e-3,3e-2,length.out=10))),labels=FALSE)		
dev.off()

pdf(file="range_planktonic_pps.pdf",width=14,height=7)
	par(cex.axis=1.5,cex.lab=1.5,mar=c(5,5,2,1),mfrow=c(1,2))
	plotBetaPPS(db=outs[[3]]$db,ppsOut=pps,sampCols=sampCols,
				xlab="Range Extent (km)",ylab="Genetic Diversity",trp=0.7,xaxt='n')
				max(z[["max.95.sea.gbif"]])
		axis(side=1,at=seq(0,1,by=5e3/max(z[["max.95.sea.gbif"]])),
			labels=round(seq(0,1,by=5e3/max(z[["max.95.sea.gbif"]]))*max(z[["max.95.sea.gbif"]])))
	discreteViolPlot(z,2,"isPlanktonic_atanypoint",xAxLabs=c("Non-planktonic","Planktonic"),sampCols=sampCols,logY=TRUE,pchMed=18,lineCol="gray98")
dev.off()

xx <- z[bioPreds][,-c(4,14)]
xxnames <- c("latitude","ecoregions","range extent",	#"ecoregions/range_size",
			 "body size","egg size","iteroparity",
			 "philopatry","spawning","planktotrophy",
			 "PLD","planktonicity","benthicity") #,"dispersive adults"
reord <- c(6:9,11:12,1:2,10,4:5,3)
xx <- xx[,reord]
xxnames <- xxnames[reord]
#mapFig <- png::readPNG("../../figures/world_map_genetic_pts-EckertIV.png")
mapFig <- png::readPNG("../../figures/world_map_genetic_pts-phy-Spilhaus-option2b.png")
#z$cladecolor[z$species=="Sargassum_muticum"] <- z$cladecolor[which(z$species=="Eukrohnia_hamata")]
pdf(file="all_predictors.pdf",width=14,height=10)
	phyViz(db=outs[[3]]$db,fit=outs[[3]]$fit,
			XX=xx,
			predNames=xxnames,
			cladeCols=unique(z$newTaxCol)[c(8,10,9,7,6,4,3,2,1,5)],
		   tree=sampPhy,xlim=c(1e-3,0.035),tipcols=z$newTaxCol,
		   valRange=NULL,adj=0.5,logX=FALSE,rounding=0.05)
#	par(mfg=c(1,1))
#		addMap2fig(mapFig=mapFig,xl=0.002,xr=0.027,yt=96,yb=57) #xr=0.027,yt=96,yb=50) #xl=0,xr=0.032,yt=96,yb=57)
dev.off()

################################
# visualize predictor correlations
################################

M <- cor(xx,use="pairwise.complete.obs",method="kendall")
row.names(M) <- xxnames
colnames(M) <- xxnames

pdf(file="predictor_corrs.pdf",width=6,height=6)
	par(cex.axis=1.5,cex.lab=1.5)
	corrplot(M,method="ellipse",diag=TRUE,type="lower",tl.col="black")
dev.off()

################################
# visualize phylogenetic correlogram
################################

load("../phy_perm_cgram.Robj")
r <- Reduce("cbind",lapply(permNullSims,"[[","r"))
lowerCI <- apply(r,1,quantile,probs=0.025)
upperCI <- apply(r,1,quantile,probs=0.975)
pdf(file="phylo_correlogram.pdf",width=7,height=7)
	plot(0,
		xlim=c(0,max(cgram$res[,1])),ylim=c(-0.15,0.2),
		type='n',
		xlab="Phylogenetic distance (my)",ylab="Correlation")
	invisible(lapply(permNullSims,function(x){
		lines(x$d,x$r,col=adjustcolor(1,0.05))})) ; 
	lines(cgram$res[,1],cgram$res[,4],col="red",lwd=3)
	lines(cgram$res[,1],lowerCI,col="black",lwd=1.5,lty=2)
	lines(cgram$res[,1],upperCI,col="black",lwd=1.5,lty=2)
	legend(x="topright",
			cex=1.25,
			legend=c("permuted trait\nnull correlogram",
		 			 "95% confidence interval",
					  "observed correlogram"),
			lty=c(1,2,1),
			col=c(adjustcolor("black",0.5),"black","red"),
			lwd=c(1.5,2,3))
dev.off()
################################
# report statistics for paper
################################
if(FALSE){

################################
# test for phylogenetic signal
################################

# blomberg's k
blom_k <- phytools::phylosig(sampPhy,z$div,test=TRUE)
pagels_lambda <- phytools::phylosig(sampPhy,z$div,method="lambda",test=TRUE)
	#plot(blom_k)
# when does the correlogram become significant relative to the null?
max(cgram$res[,1][which(cgram$res[,4]>upperCI)])

# diversity statistics
range(z$div)
mean(z$div)
#phylo corrected mean
I <- matrix(rep(1,nrow(z)),nrow(z),1)
C_inv <- solve(phyStr)
1/(t(I) %*% C_inv %*% I) * t(I) %*% C_inv %*% z$div
mean(solve(chol(phyStr)) %*% z$div)

################################
# report R^2 values
################################
load("../partB_outs.Robj")
rangeR2l <- bayesR2(outs[[3]],latent=TRUE)
plankR2l <- bayesR2(outs[[12]],latent=TRUE)
mean(rangeR2l)
quantile(rangeR2l,c(0.025,0.975))
mean(plankR2l)
quantile(plankR2l,c(0.025,0.975))



################################
# export effect sizes
################################

outTab <- reportEffectSizes(outs)
cat(capture.output(outTab), file = 'effect_size_summary.txt', sep = '\n')
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
}