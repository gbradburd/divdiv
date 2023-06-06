mvnPhyReg <- stan_model(model_code=mvnPhyReg)

################################
# analyze s with one predictor at a time
#	MVN model, unscaled s
################################

db <- lapply(1:nrow(predictors),
			function(i){
				makeDB(predictors[i,],z$s,phyStr)
		})
fits <- lapply(1:nrow(predictors),
			function(i){
				sampling(object=mvnPhyReg,
					data=db[[i]],
					iter=nIter,
					thin=nIter/500,
					chains=1,
					control=setNames(list(15),"max_treedepth"))
			})
names(db) <- predNames
names(fits) <- predNames
out <- list("db"=db,"fits"=fits)
save(out,file="mvn_s_unscld.Robj")

pdf(file="mvn_s_unscld.pdf",width=12,height=10)
	par(mfrow=c(3,4)) ; for(i in 1:11){mvnPPS(out$db[[i]],out$fit[[i]],500,predName=names(out$fit)[i])}
dev.off()

################################
# analyze s with one predictor at a time
#	MVN model, scaled s
################################

s <- (z$s - min(z$s))/max((z$s - min(z$s)))
db <- lapply(1:nrow(predictors),
			function(i){
				makeDB(predictors[i,],s,phyStr)
		})
fits <- lapply(1:nrow(predictors),
			function(i){
				sampling(object=mvnPhyReg,
					data=db[[i]],
					iter=nIter,
					thin=nIter/500,
					chains=1,
					control=setNames(list(15),"max_treedepth"))
			})

unscale <- list("mn" = min(z$s),"mx" = max((z$s - min(z$s))))
names(db) <- predNames
names(fits) <- predNames
out <- list("db"=db,"fits"=fits,"unscale"=unscale)
save(out,file="mvn_s_scld.Robj")

pdf(file="mvn_s_scld.pdf",width=12,height=10)
	par(mfrow=c(3,4)) ; for(i in 1:11){mvnPPS(out$db[[i]],out$fit[[i]],500,predName=names(out$fit)[i])}
dev.off()

################################
# analyze Nbhd with one predictor at a time
#	MVN model
################################

db <- lapply(1:nrow(predictors),
			function(i){
				makeDB(predictors[i,],z$nbhd,phyStr)
		})
fits <- lapply(1:nrow(predictors),
			function(i){
				sampling(object=mvnPhyReg,
					data=db[[i]],
					iter=nIter,
					thin=nIter/500,
					chains=1,
					control=setNames(list(15),"max_treedepth"))
			})
names(db) <- predNames
names(fits) <- predNames
out <- list("db"=db,"fits"=fits)
save(out,file="mvn_Nbhd.Robj")

pdf(file="mvn_Nbhd.pdf",width=12,height=10)
	par(mfrow=c(3,4)) ; for(i in 1:11){mvnPPS(out$db[[i]],out$fit[[i]],500,predName=names(out$fit)[i])}
dev.off()


load("mvn_s_unscld.Robj")
pdf(file="mvn_s_unscld.pdf",width=12,height=10)
	par(mfrow=c(3,4)) ; for(i in 1:11){mvnPPS(out$db[[i]],out$fit[[i]],500,predName=names(out$fit)[i])}
dev.off()
pdf(file="betas_mvn_s_unscld.pdf",width=12,height=10)
	postBetaPlot(out=out,predNames=names(out$fits),reorder=TRUE,cols=NULL)
dev.off()

load("mvn_s_scld.Robj")
pdf(file="mvn_s_scld.pdf",width=12,height=10)
	par(mfrow=c(3,4)) ; for(i in 1:11){mvnPPS(out$db[[i]],out$fit[[i]],500,predName=names(out$fit)[i])}
dev.off()
pdf(file="betas_mvn_s_scld.pdf",width=12,height=10)
	postBetaPlot(out=out,predNames=names(out$fits),reorder=TRUE,cols=NULL)
dev.off()

load("mvn_Nbhd.Robj")
pdf(file="mvn_Nbhd.pdf",width=12,height=10)
	par(mfrow=c(3,4)) ; for(i in 1:11){mvnPPS(out$db[[i]],out$fit[[i]],500,predName=names(out$fit)[i])}
dev.off()
pdf(file="betas_mvn_Nbhd.pdf",width=12,height=10)
	postBetaPlot(out=out,predNames=names(out$fits),reorder=TRUE,cols=NULL)
dev.off()


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