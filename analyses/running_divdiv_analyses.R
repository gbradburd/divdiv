################################################################
################################################################
# Running DivDiv analyses
################################################################
################################################################


################################
# load libraries and define 
#	useful functions
################################
library(rstan)
library(ape)

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
source("trait_mod_stan_block.R")
mvn <- stan_model(model_code=mvn)

################################
# get phylo structure,
# load master dataframe,
#	create dataBlocks
################################
z <- read.csv("../data/master_df.csv",header=TRUE,stringsAsFactors=FALSE)
z <- z[-which(z$species=="Pocillopora_damicornis"),]

load("../data/phylo/divdiv_phy_from_timetreebeta5.Robj")
sampPhy <- ape::keep.tip(sampPhy,gsub("_"," ",z$species))
phyStr <- ape::vcv(sampPhy,corr=TRUE)

predictors <- rbind(z[["meanlat.gbif"]],
					z[["n_ECOREGIONS.all"]],
					z[["maxgbif.sea_km"]],
					z[["Body_Size"]],
					z[["Fecundity_EggSize"]],
					z[["Generational_Structure"]],
					z[["ReturnToSpawningGround"]],
					z[["Spawning_mode"]],
					z[["Larval_feeding"]],
					z[["PLD_point2"]],
					z[["isPlanktonic_atanypoint"]])

tmp <- predictors
# fill in missing data w/ grand mean for that predictor
md <- which(is.na(predictors),arr.ind=TRUE)
for(i in 1:nrow(md)){
	predictors[md[i,1],md[i,2]] <- mean(tmp[md[i,1],],na.rm=TRUE)
}

# modeling collecting phase pi
s <- (z$s - min(z$s))/max((z$s - min(z$s)))
db_s <- list("N" = nrow(phyStr),
		   	 "Y" = s,
		   	 "nX"= nrow(predictors),
		   	 "X" = predictors,
		     "relMat" = phyStr)

# modeling neighborhood sizes
db_nbhd <- list("N" = nrow(phyStr),
		   	    "Y" = z$nbhd,
		   	    "nX"= nrow(predictors),
		   	    "X" = predictors,
		        "relMat" = phyStr)

fit_s <- sampling(object=mvn,
					data=db_s,
					iter=2e3,
					chains=2)

fit_nbhd <- sampling(object=mvn,
						data=db_nbhd,
						iter=2e3,
						chains=2)

b_s <- rstan::extract(fit_s,"beta",permute=FALSE)
b_nbhd <- rstan::extract(fit_nbhd,"beta",permute=FALSE)

par(mfrow=c(1,2))
	matplot(b_s[,1,],type='l')
	matplot(b_s[,2,],type='l')

par(mfrow=c(1,2))
	matplot(b_nbhd[,1,],type='l')
	matplot(b_nbhd[,2,],type='l')