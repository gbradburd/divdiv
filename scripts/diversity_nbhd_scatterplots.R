################################################################
################################################################
#	script for generating species scatterplots of diversity & nbhd size
################################################################
################################################################

################################
#	read in metadata
################################
z <- read.csv("~/Dropbox/divdiv/data/master_df.csv")
z$species[which(z$species=="Seriola_lalandi_dorsalis")] <- "Seriola_lalandi"

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

