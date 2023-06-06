library(rstan)
library(ape)
library(MASS)

load("../../data/phylo/divdiv_phy_from_timetreebeta5.Robj")
phyStr <- ape::vcv(phy,corr=TRUE)

source("../../analyses/trait_mod_stan_blocks.R")
#betaPhyReg <- stan_model(model_code=betaPhyReg)
expPhyReg <- stan_model(model_code=expPhyReg)

n <- nrow(phyStr)
x <- rnorm(n,sd=10)
alpha <- runif(1)
gamma <- rnorm(1)
bt <- rnorm(1,sd=0.1)
mn <- gamma + bt*x
theta <- MASS::mvrnorm(1,mu=mn,Sigma=alpha*phyStr)
lnk_mn <- exp(theta)
sim <- rexp(n,lnk_mn)

# par(mfrow=c(1,3))
# plot(x,sim)
# plot(x,lnk_beta_mean,ylim=c(0.99,1))
# plot(x,sim,ylim=c(0.99,1))



db <- list("N" = nrow(phyStr),
		   "Y" = sim,
		   "nX"= 1,
		   "X" = matrix(x,nrow=1,ncol=nrow(phyStr)),
		   "relMat" = phyStr)

fit <- sampling(object=expPhyReg,
				data=db,
				iter=5e3,
				chains=1)

alpha_e <- extract(fit,pars="alpha",inc_warmup=FALSE,permute=FALSE)
gamma_e <- extract(fit,pars="gamma",inc_warmup=FALSE,permute=FALSE)
beta_e <- extract(fit,pars="beta",inc_warmup=FALSE,permute=FALSE)
lp <- get_logposterior(fit,inc_warmup=FALSE)[[1]]

par(mfrow=c(1,4))
plot(lp)
plot(alpha_e) ; abline(h=alpha,col="red")
plot(gamma_e) ; abline(h=gamma,col="red")
plot(beta_e) ; abline(h=bt,col="red")



