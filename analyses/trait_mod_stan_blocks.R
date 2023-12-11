################################################################
################################################################
#	modeling trait data while controlling for relatedness
################################################################
################################################################

mvnPhyReg <- 
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
	vector[N] Y;				// vector of response data
	int<lower=0> nX;					// number of predictor variables
	vector[N] X[nX];				// array of predictors with N rows and nX columns
	matrix[N,N] relMat; 			// relatedness matrix
}
parameters {
	real<lower=0> alpha;			// effect of relatedness on similarity in response variable
	real gamma;						// mean response variable
	real beta[nX];					// effects of predictors on response variable
}
transformed parameters {
	vector[N] meanVec;						// mean response variable
	meanVec = gamma + getMean(N,nX,X,beta);
}
model {
	alpha ~ normal(0,10);									// prior on effect of relatedness matrix
	gamma ~ normal(0,100);										// prior on global mean
	beta ~ normal(0,10);									// prior on effects of predictors
	Y ~ multi_normal(meanVec,alpha * relMat);		// MVN likelihood of transformed p(response variable)
}
"

betaPhyReg <- 
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
	vector[N] Y;						// vector of response data (ranging between 0 and 1)
	int<lower=0> nX;				// number of predictor variables
	vector[N] X[nX];				// array of predictors with N rows and nX columns
	matrix[N,N] relMat; 			// relatedness matrix
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
	vector[N] theta;				// latent variable
}
transformed parameters {
	vector[N] meanVec = gamma + getMean(N,nX,X,beta);	// logit mean response variable
	vector[N] mu = inv_logit(theta);						// mean response variable
	vector[N] shape1 = mu * phi;						// beta shape1 parameter
	vector[N] shape2 = (1-mu) * phi;					// beta shape2 parameter
}
model {
	alpha ~ normal(0,10);									// prior on effect of relatedness matrix
	phi ~ exponential(varTrY);								// prior on effect of dispersion
	gamma ~ normal(mnTrY,1);									// prior on global mean
	beta ~ normal(0,10);									// prior on effects of predictors
	theta ~ multi_normal(meanVec,alpha * relMat);			// MVN likelihood of transformed p(response variable)
	Y ~ beta(shape1,shape2);									// beta likelihood of response data
}
"


expPhyReg <- 
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
	vector[N] Y;				// vector of response data
	int<lower=0> nX;					// number of predictor variables
	vector[N] X[nX];				// array of predictors with N rows and nX columns
	matrix[N,N] relMat; 			// relatedness matrix
}
parameters {
	real<lower=0> alpha;			// effect of relatedness on similarity in response variable
	real gamma;						// mean response variable
	real beta[nX];					// effects of predictors on response variable
	vector[N] theta;
}
transformed parameters {
	vector[N] meanVec = gamma + getMean(N,nX,X,beta);
	vector[N] lambda = exp(theta);						// mean response variable
}
model {
	alpha ~ normal(0,10);									// prior on effect of relatedness matrix
	gamma ~ normal(0,10);									// prior on global mean
	beta ~ normal(0,10);									// prior on effects of predictors
	theta ~ multi_normal(meanVec,alpha * relMat);			// MVN likelihood of transformed p(response variable)
	Y ~ exponential(lambda);
}
"

betaLoo <- 
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