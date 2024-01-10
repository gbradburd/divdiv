library(rstan)

betaPhyReg2 <- 
"functions {
	vector getMean(int N, int nX, vector[] X, real[] beta) {
		vector[N] myMean;
		myMean = rep_vector(0,N);
		for(i in 1:nX){
			myMean += beta[i] * X[i];
		}
		return(myMean);
	}
}
data {
	int<lower=2> N; 	  			// number of samples
	int<lower=2> Nstar;			// when inferring a missing value, Nstar=N+1
	vector[N] Y;						// vector of response data (ranging between 0 and 1)
	int<lower=0> nX;				// number of predictor variables
	vector[Nstar] X[nX];				// array of predictors with N rows and nX columns with focal sample as last row/column
	matrix[Nstar,Nstar] relMat; 			// relatedness matrix with focal sample as the last row/column
}
transformed data {
	real mnTrY=logit(mean(Y));			// mean transformed response variable
	real varTrY=variance(Y);					// mean transformed response variable
}
parameters {
	real<lower=0> alpha;			// effect of relatedness on similarity in response variable
	real gamma;						// mean response variable
	real<lower=0> phi;			// dispersion parameter of beta distribution
	real beta[nX];					// effects of predictors on response variable
	vector[Nstar] theta;				// latent variable
	vector<lower=0,upper=1>[Nstar] y;	// expected value of species' response variable
	real<lower=0> sigma;				// variance of Y around y
}
transformed parameters {
	vector[Nstar] meanVec = gamma + getMean(Nstar,nX,X,beta);	// logit mean response variable
	vector[Nstar] mu = inv_logit(theta);						// mean response variable
	vector[Nstar] shape1 = mu * phi;						// beta shape1 parameter
	vector[Nstar] shape2 = (1-mu) * phi;					// beta shape2 parameter
}
model {
	alpha ~ normal(0,10);									// prior on effect of relatedness matrix
	phi ~ exponential(varTrY);								// prior on effect of dispersion
	gamma ~ normal(mnTrY,1);									// prior on global mean
	beta ~ normal(0,10);									// prior on effects of predictors
	theta ~ multi_normal(meanVec,alpha * relMat);			// MVN likelihood of transformed p(response variable)
	y ~ beta(shape1,shape2);									// beta likelihood of response data
	sigma ~ normal(0,1);									// prior on sigma
	for(i in 1:N) Y[i] ~ normal(y[i],sigma);				// observed data are drawn from a normal distribution around their predicted value
}
"

mod <- stan_model(model_code=betaPhyReg2)


simData <- function(N,nX){
	X <- matrix(rnorm(N*nX),nrow=nX,ncol=N)
	X[,N] <- 0
	phy <- ape::rtree(N)
	phyCov <- ape::vcv(phy)
	alpha <- abs(rnorm(1))
	gamma <- rnorm(1,mean=-6,sd=0.01)
	beta <- matrix(rnorm(nX,mean=0,sd=1),nrow=nX,ncol=1)
	phi <- rexp(1,rate=1/100)
	meanVec <- gamma + t(beta) %*% X
	theta <- MASS::mvrnorm(n=1,mu=meanVec,Sigma=alpha*phyCov)
	mu <- 1/(1+exp(-theta))
	shape1 <- mu*phi
	shape2 <- (1-mu)*phi
	y <- rbeta(n=N,shape1=shape1,shape2=shape2)
	sigma <- abs(rnorm(1,mean=0,0.001))
	Y <- rnorm(n=N,mean=y,sd=sigma)
	dataBlock <- list("N"=N-1,
					  "Nstar"=N,
					  "Y"=Y[1:(N-1)],
					  "nX"=nX,
					  "X"=X,
					  "relMat"=phyCov)
	simPars <- list("phy"=phy,
					"alpha"=alpha,
					"gamma"=gamma,
					"beta"=beta,
					"phi"=phi,
					"meanVec"=meanVec,
					"theta"=theta,
					"mu"=mu,
					"shape1"=shape1,
					"shape2"=shape2,
					"y"=y,
					"sigma"=sigma,
					"Y"=Y)
	out <- list("db"=dataBlock,
				"pars"=simPars)
	return(out)
}

sim <- simData(50,2)
fit <- rstan::sampling(object=mod,
						data=sim$db,
						iter=2e3,
						thin=2e3/500,
						chains=1,
						control=setNames(list(15),"max_treedepth"))
plot(get_logposterior(fit)[[1]])
alpha <- rstan::extract(fit,"alpha",inc_warmup=FALSE,permute=FALSE)
gamma <- rstan::extract(fit,"gamma",inc_warmup=FALSE,permute=FALSE)
beta <- rstan::extract(fit,"beta",inc_warmup=FALSE,permute=FALSE)
phi <- rstan::extract(fit,"phi",inc_warmup=FALSE,permute=FALSE)
sigma <- rstan::extract(fit,"sigma",inc_warmup=FALSE,permute=FALSE)
y <- rstan::extract(fit,"y",inc_warmup=FALSE,permute=FALSE)[,1,]


layout(matrix(c(1:5,rep(6,4),7),nrow=2,ncol=5,byrow=TRUE))
	hist(alpha)
		abline(v=sim$pars$alpha,col="red")
	hist(gamma)
		abline(v=sim$pars$gamma,col="red")
	hist(beta)
		abline(v=sim$pars$beta,col="red")
	hist(phi)
		abline(v=sim$pars$phi,col="red")
	hist(sigma)
		abline(v=sim$pars$sigma,col="red")
	plot(sim$db$X[1:sim$db$N],sim$db$Y,type='n',
			xlab="predictor",ylab="response")
		invisible(
			lapply(1:sim$db$N,
				function(n){
					qnt <- quantile(y[,n],c(0.025,0.975))
					segments(x0=sim$db$X[n],
							 x1=sim$db$X[n],
							 y0=qnt[1],
							 y1=qnt[2],
							 col=adjustcolor("red",0.9))
					}))
	points(sim$db$X[1:sim$db$N],sim$db$Y,pch=20)
	hist(y[,sim$db$Nstar],xlab="values",main="predicting missing\ndatum")
		abline(v=sim$pars$Y[sim$db$Nstar],col="red")








