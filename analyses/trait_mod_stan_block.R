################################################################
################################################################
#	modeling trait data while controlling for relatedness
################################################################
################################################################

library(rstan)

mvn <- 
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
	real mu;						// mean response variable
	real beta[nX];					// effects of predictors on response variable
}
transformed parameters {
	vector[N] meanVec;						// mean response variable
	meanVec = mu + getMean(N,nX,X,beta);
}
model {
	alpha ~ normal(0,10);									// prior on effect of relatedness matrix
	mu ~ normal(0,10);										// prior on global mean
	beta ~ normal(0,10);									// prior on effects of predictors
	Y ~ multi_normal(meanVec,alpha * relMat);		// MVN likelihood of transformed p(response variable)
}
"