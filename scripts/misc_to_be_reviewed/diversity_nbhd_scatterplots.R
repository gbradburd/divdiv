################################################################
################################################################
#	script for generating species scatterplots of diversity & nbhd size
################################################################
################################################################

################################
#	read in metadata
################################
z <- read.csv("~/Dropbox/divdiv/data/master_df.csv")

z <- z[,c("species","cladecolor","s","nbhd",
			"meanlat.gbif",
			"n_ECOREGIONS.all",
			"max.95.sea.gbif",
			"Body_Size",
			"Fecundity_EggSize",
			"Generational_Structure",
			"ReturnToSpawningGround",
			"Spawning_mode",
			"Larval_feeding",
			"PLD_point2",
			"isPlanktonic_atanypoint",
			"ratio.sea.95",
			"n_samples",
			"mean_raw_read_cnt",
			"read_length",
			"mean_locus_depth",
			"n.totalsnps")]

z[["log(deep-time diversity)"]] <- log(1-z$s)


pdf(file="../figures/variable_scatterplots.pdf",width=12,height=5,pointsize=11)
for(i in 4:ncol(z)){
	par(mar=c(8,5,3,1),bty="n")
	x <- z[,i]
	anyMD <- any(is.na(x))
	if(anyMD){
		md <- which(is.na(x))
		mdSp <- z$species[md]
		x <- x[-md]
		sp <- z$species[-md]
		spCols <- z$cladecolor[-md]
	} else {
		sp <- z$species
		spCols <- z$cladecolor
	}
	plot(x[order(x)],
			xlab="",ylab="",
			xaxt='n',pch=19,cex=2,main=names(z)[i],
			col=adjustcolor(spCols[order(x)],0.5))
	segments(x0=1:length(x),y0=min(x)-1,x1=1:length(x),y1=x[order(x)],col=adjustcolor(spCols[order(x)],0.5),lty=2)
	text(1:length(sp),par("usr")[3],labels=gsub("_"," ",sp[order(x)]),srt=50,adj=c(1.1,1.1),xpd=TRUE,cex=0.7,col=adjustcolor(spCols[order(x)]))
	if(anyMD){
		legend(x="topleft",col="black",title="species with NA values",pch=NULL,legend=gsub("_"," ",mdSp),cex=0.5)
	}
}
dev.off()

pdf(file="../figures/predictor_corrs.pdf",width=15,height=15,pointsize=11)
	par(mar=c(8,5,3,1),mfrow=c(4,5))
	for(i in 4:ncol(z)){
		for(j in c(4:ncol(z))[-(i-3)]){
			x1 <- z[[i]]
			x2 <- z[[j]]
			message(sprintf("plotting %s against %s",i,j))
			if(any(is.na(c(x1,x2)))){
				md <- unique(c(which(is.na(x1)),which(is.na(x2))))
				x1 <- x1[-md]
				x2 <- x2[-md]
				sp <- z$species[-md]
				spCols <- z$cladecolor[-md]
			} else {
				sp <- z$species
				spCols <- z$cladecolor
			}
			x1x2cor <- round(cor(x1,x2),3)
			plot(x1,x2,
					pch=19,cex=2,
					main=sprintf("%s x %s, cor = %s",
									names(z)[i],
									names(z)[j],
									x1x2cor),
					col=adjustcolor(spCols,0.5),
					xlab=names(z)[i],
					ylab=names(z)[j])
			if(abs(x1x2cor) > 0.5){
				box(lwd=2,col="red")
			}
		}
		plot.new()
		plot.new()
	}
dev.off()

pdf(file="../figures/pi_across_spp.pdf",width=10,height=5,pointsize=11)
par(mar=c(8,5,3,1))
	plot(log(1-z$s[order(log(1-z$s))]),
			xlab="",ylab=expression(paste("log(deep-time ",pi,")")),
			xaxt='n',pch=19,cex=2,main="Genetic diversity",
			ylim=range(log(1-z$s))+c(-0.001,0.001),col=adjustcolor(z$cladecolor[order(z$s)],0.5))
	text(1:length(z$s),par("usr")[3],labels=gsub("_"," ",z$species[order(log(1-z$s))]),srt=50,adj=c(1.1,1.1),xpd=TRUE,cex=0.7)
dev.off()

pdf(file="../figures/nbhd_across_spp.pdf",width=10,height=5,pointsize=11)
par(mar=c(8,5,3,1))
	plot(z$nbhd[order(z$nbhd)],
			xlab="",ylab="neighborhood size",
			xaxt='n',pch=19,cex=2,main="Neighborhood size",
			ylim=range(z$nbhd)+c(-3,5),col=adjustcolor(z$cladecolor[order(z$nbhd)],0.5))
	text(1:length(z$nbhd),par("usr")[3],labels=gsub("_"," ",z$species[order(z$nbhd)]),srt=50,adj=c(1.1,1.1),xpd=TRUE,cex=0.7)
dev.off()

pdf(file="../figures/nbhd_vs_pi.pdf",width=7,height=7,pointsize=11)
plot(z$nbhd,log(1-z$s),col=adjustcolor(z$cladecolor[order(z$nbhd)],0.5),pch=19,cex=2,
	xlab="neighborhood size",
	ylab=expression(paste("log(deep-time ",pi,")")))
dev.off()

# phylogenetic distance vs. predictor and response variable distance
load("../data/phylo/divdiv_phy_from_timetreebeta5.Robj")
sampPhy <- ape::keep.tip(phy,gsub("_"," ",z$species))
phyStr <- ape::cophenetic.phylo(sampPhy)
upTri <- upper.tri(phyStr)
colmat1 <- matrix(z$cladecolor,nrow=nrow(z),ncol=nrow(z))[upTri]
colmat2 <- matrix(z$cladecolor,nrow=nrow(z),ncol=nrow(z),byrow=TRUE)[upTri]

sp <- gsub("_"," ",z$species)
z <- z[match(row.names(phyStr),sp),]


pdf(file="../figures/phylodist.pdf",width=7,height=7,pointsize=11)
	par(cex.lab=2,mar=c(6,6,5,1))
	for(i in 4:ncol(z)){
		plot(phyStr[upTri],fields::rdist(z[,i])[upTri],
				pch=21,bg=adjustcolor(colmat1,0.4),col=adjustcolor(colmat2,0.4),
				main=names(z)[i],
				xlab="phylogenetic distance",ylab="trait distance")
	}
dev.off()

#graveyard
if(FALSE){
spnames <- gsub("_"," ",z$species)
tax <- unique(z$taxclade)
ccols <- unique(z$cladecolor)
plot(z$s,type='n',pch=19,col=z$cladecolor,ylab="genetic similarity")
	xx <- 0
	for(i in 1:length(tax)){
		inTax <- which(z$taxclade==tax[i])
		text(x=xx + 1:length(inTax),y=z$s[inTax],
			 labels=spnames[inTax],
			 col=ccols[i],cex=0.5)
		xx <- xx + length(inTax)
	}	
}

