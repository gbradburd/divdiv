library(rstan)
library(ape)
library(MASS)
source("../../analyses/trait_mod_stan_blocks.R")
betaPhyReg <- stan_model(model_code=betaPhyReg)

shapes2muphi <- function(shape1, shape2) {
  mu <- shape1 / (shape1 + shape2)
  phi <- shape1 + shape2
  return(list(mu = mu, phi = phi))
}

muphi2shapes <- function(mu, phi) {
  shape1 <- mu * phi
  shape2 <- (1 - mu) * phi
  return(list(shape1 = shape1, shape2 = shape2))
}

invlogit <- function(x){
	return(exp(x)/(1+exp(x)))
}

load("../data/phylo/divdiv_phy_from_timetreebeta5.Robj")
phyStr <- ape::vcv(phy,corr=TRUE)

bad <- TRUE
counter <- 0
while(bad){
	n <- nrow(phyStr)
	x <- runif(n,0,100)
	alpha <- runif(1)
	gamma <- rnorm(n=1,mean=2)
	phi <- rexp(1)
	bt <- rnorm(1,sd=0.1)
	mvn_mean <- gamma + bt*x
	theta <- MASS::mvrnorm(1,mu=mvn_mean,Sigma=alpha*phyStr)
	lnk_beta_mean <- invlogit(theta)
	ab <- muphi2shapes(lnk_beta_mean,phi)
	sim <- rbeta(n,ab[[1]],ab[[2]])
	if((min(sim) > 0.85)){ #!any(sim==1) & !any(sim==0) & 
		bad <- FALSE
	}
	counter <- counter + 1
	print(counter)
}

# par(mfrow=c(1,3))
# plot(x,sim)
# plot(x,lnk_beta_mean,ylim=c(0.99,1))
# plot(x,sim,ylim=c(0.99,1))

sim[which(sim==1)] <- 1-1e-6
sim[which(sim==0)] <- 1+1e-6
db <- list("N" = nrow(phyStr),
		   "Y" = sim,
		   "nX"= 1,
		   "X" = matrix(x,nrow=1,ncol=nrow(phyStr)),
		   "relMat" = phyStr)

fit <- sampling(object=betaPhyReg,
				data=db,
				iter=2e3,
				chains=1)

alpha_e <- extract(fit,pars="alpha",inc_warmup=FALSE,permute=FALSE)
gamma_e <- extract(fit,pars="gamma",inc_warmup=FALSE,permute=FALSE)
phi_e <- extract(fit,pars="phi",inc_warmup=FALSE,permute=FALSE)
beta_e <- extract(fit,pars="beta",inc_warmup=FALSE,permute=FALSE)

par(mfrow=c(1,4))
plot(alpha_e) ; abline(h=alpha,col="red")
plot(gamma_e) ; abline(h=gamma,col="red")
plot(phi_e) ; abline(h=phi,col="red")
plot(beta_e) ; abline(h=bt,col="red")




invLogit <- function(x){
	return(exp(x)/(1+exp(x)))
}

betaPPS <- function(db,fit,nPPS){
	nIter <- length(get_logposterior(fit)[[1]])
	sampledIter <- sample(1:(nIter/2),nPPS)
	beta <- extract(fit,"beta",permute=FALSE)[,1,1]
	sig <- 1-2*(abs(0.5-ecdf(beta)(0)))
	siglwd <- ifelse(sig<0.05,4,0.5)
	beta <- beta[sampledIter]
	alpha <- extract(fit,"alpha",permute=FALSE)[sampledIter,1,1]
	gamma <- extract(fit,"gamma",permute=FALSE)[sampledIter,1,1]
	phi <- extract(fit,"phi",permute=FALSE)[sampledIter,1,1]
	mvnMean <- lapply(1:nPPS,function(i){c(gamma[i] + beta[i]*db$X)})
	theta <- lapply(1:nPPS,function(i){MASS::mvrnorm(1,mu=mvnMean[[i]],Sigma=alpha[i]*db$relMat)})
	mu <- lapply(theta,function(x){invLogit(x)})
	shape1 <- lapply(1:nPPS,function(i){phi[[i]]*mu[[i]]})
	shape2 <- lapply(1:nPPS,function(i){phi[[i]]*(1-mu[[i]])})
#	shape1 <- extract(fit,"shape1",permute=FALSE)[sampledIter,1,]
#	shape2 <- extract(fit,"shape2",permute=FALSE)[sampledIter,1,]
	pps <- lapply(1:nPPS,function(i){rbeta(db$N,shape1[[i]],shape2[[i]])})
	ppsCIs <- getPPSci(db,pps)
	plot(db$X,db$Y,
		ylim=range(c(db$Y,unlist(ppsCIs))),
		xlab="predictor",ylab="response",
		type='n',main=sprintf("posterior predictive sampling (p=%s)",round(sig,3)))
		invisible(
			lapply(1:db$N,
				function(n){
					segments(x0=db$X[n],
							 y0=ppsCIs[[n]][1],
							 x1=db$X[n],
							 y1=ppsCIs[[n]][2],
							 lwd=0.75)
				})
		)
	points(db$X,db$Y,col="red",pch=19,cex=1)
	lines(0:max(db$X),invLogit(mean(beta)*(0:max(db$X))+mean(gamma)),lwd=siglwd)
	box(lwd=2,col=ifelse(sig<0.05,"red","black"))
}

getPPSci <- function(db,pps){
	ppsByN <- lapply(1:db$N,
				function(n){
					quantile(unlist(lapply(pps,"[[",n)),c(0.025,0.975))
				})
	return(ppsByN)
}

betaPPS(db,fit,nPPS=100)