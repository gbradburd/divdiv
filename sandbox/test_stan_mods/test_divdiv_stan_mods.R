################################################################
################################################################
#	testing DivDiv stan models
################################################################
################################################################

################################
# load libraries and define 
#	useful functions
################################
library(rstan)
library(ape)
library(MASS)

colFunc <- function (x, cols, nCols, valRange){
    if (is.null(valRange)) {
        valRange <- c(min(x), max(x))
    }
    cols <- (grDevices::colorRampPalette(cols))(nCols)[findInterval(x, 
        seq(valRange[1], valRange[2], length.out = nCols))]
    return(cols)
}

################################
# compile rstan models
################################
source("../../analyses/trait_mod_stan_block.R")
mvn <- stan_model(model_code=mvn)
bin <- stan_model(model_code=bin)

################################
# get phylo structure
################################
load("../../data/phylo/divdiv_phy_from_timetreebeta5.Robj")
phyStr <- ape::vcv(phy,corr=TRUE)

################################
# simulate a brownian character
################################
simTrt1 <- lapply(1:100,function(i){MASS::mvrnorm(1,mu=rep(0,nrow(phyStr)),Sigma=phyStr)})
simTrt2 <- lapply(1:100,function(i){MASS::mvrnorm(1,mu=rep(0,nrow(phyStr)),Sigma=0.5*phyStr)})
simTrt3 <- lapply(1:100,function(i){MASS::mvrnorm(1,mu=rep(0,nrow(phyStr)),Sigma=0.1*phyStr)})
ut <- upper.tri(phyStr,diag=TRUE)

st1dist <- lapply(simTrt1,function(x){fields::rdist(x)[ut]})
st2dist <- lapply(simTrt2,function(x){fields::rdist(x)[ut]})
st3dist <- lapply(simTrt3,function(x){fields::rdist(x)[ut]})

png(file="phylo_sim_test.png",res=100,width=15*100,height=7*100)
par(mfrow=c(1,3))
	plot(0,type='n',xlab="phylo distance",ylab="trait distance",xlim=c(0,1),ylim=range(st1dist),main="sigma = 1")
		invisible(lapply(1:100,function(i){points(phyStr[ut],st1dist[[i]],pch=20,col=adjustcolor(1,0.01))}))
	plot(0,type='n',xlab="phylo distance",ylab="trait distance",xlim=c(0,1),ylim=range(st2dist),main="sigma = 0.5")
		invisible(lapply(1:100,function(i){points(phyStr[ut],st2dist[[i]],pch=20,col=adjustcolor(1,0.01))}))
	plot(0,type='n',xlab="phylo distance",ylab="trait distance",xlim=c(0,1),ylim=range(st3dist),main="sigma = 0.1")
		invisible(lapply(1:100,function(i){points(phyStr[ut],st3dist[[i]],pch=20,col=adjustcolor(1,0.01))}))
dev.off()

################################
# simulate continuous data with 1 predictor
################################

x <- runif(nrow(phyStr))

beta <- rexp(1,0.5)
simTrt <- MASS::mvrnorm(1,mu=beta*x,Sigma=phyStr)
plot(x,simTrt,pch=19)

db <- list("N" = nrow(phyStr),
		   "Y" = simTrt,
		   "nX"= 1,
		   "X" = matrix(x,nrow=1,ncol=nrow(phyStr)),
		   "relMat" = phyStr)

fit <- sampling(object=mvn,
				data=db,
				iter=2e3,
				chains=1)
b <- rstan::extract(fit,"beta",permute=FALSE)
plot(b,type='l') ; abline(h=beta,col="red")


################################
# simulate continuous data with 2 predictors
################################

x1 <- runif(nrow(phyStr))
x2 <- runif(nrow(phyStr))

beta1 <- rexp(1,rate=0.1)
beta2 <- rexp(1,rate=0.1)
simTrt <- MASS::mvrnorm(1,mu=beta1*x1+beta2*x2,Sigma=phyStr)
par(mfrow=c(1,2))
	plot(x1,simTrt,col=colFunc(x2,cols=c("blue","red"),nCols=nrow(phyStr),valRange=range(x2)),pch=19,main="colored by x2")
	plot(x2,simTrt,col=colFunc(x1,cols=c("blue","red"),nCols=nrow(phyStr),valRange=range(x1)),pch=19,main="colored by x1")

db <- list("N" = nrow(phyStr),
		   "Y" = simTrt,
		   "nX"= 2,
		   "X" = rbind(x1,x2),
		   "relMat" = phyStr)

fit <- sampling(object=mvn,
				data=db,
				iter=2e3,
				chains=1)
b <- rstan::extract(fit,"beta",permute=FALSE)
matplot(b[,1,],type='l',col="black") ; abline(h=beta1,col="red") ; abline(h=beta2,col="red")