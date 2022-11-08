#idea: trying to make a density ridge plot to show leave-one-out fit	

library(ggplot)
library(ggridges)

rm(list = ls())
gc()

setwd("/Users/rachel/Desktop/DivDiv/divdiv/sandbox/troubleshooting")


	#load list of datablocks for leave one out analyses
	load("db_MVN_unscaled_s.Robj")
	#focus on a single analysis using body size as a predictor
	db <- db[[4]]
	#load the associated LOO conditional normal (CN) samples
	load("s_mvn_body_size_loo_CNsamples.Robj")
	CNsamples <- looCNsamples
	# species names
	spNames <- row.names(db$relMat)
	#describe conditional normal model fit for each species
	cnFits <- unlist(
				lapply(1:db$N,
					function(n){
						2*(abs(0.5-ecdf(CNsamples[[n]])(db$Y[n])))
				}))
	# turn it all into a disgusting dataframe so that it "works" with ggplot
	df <- Reduce("rbind",
			lapply(1:length(CNsamples[[1]]),
				function(i){data.frame("species"=spNames,
									  "cns"=unlist(lapply(CNsamples,"[[",i)),
									  "cnf"=cnFits,
									  "s"=db$Y)}))
	# order this most glorious beauty of a way to store data so lines match ridges
	df <- df %>% arrange(species) %>% group_by(species) %>% mutate(o = cur_group_id())
	
	# try to make a density ridge figure, w/ density plots colored by how far off the CN fit is
	df %>% ggplot() + 
	  geom_density_ridges(aes(x=cns,y=species,fill=cnf), scale=2,rel_min_height=0.01) +
		geom_segment(aes(x=s, xend=s, y=o, yend=o+2), colour="red") + 
	  xlim(c(0.97,1.02)) +
	  theme_bw() +
	  theme(panel.grid = element_blank()) +
	  labs(y = "Species")
	
	
	
	
#graveyard -------
	df %>% ggplot() +
	  geom_density(aes(x=cns, group=species)) +
	  facet_wrap(~species, ncol = 1) +
	  #geom_segment(aes(x=s, xend=s, y=o, yend=o+2), colour="red") + 
	  xlim(c(0.97,1.02)) +
	  theme_bw() +
	  theme(panel.grid = element_blank())
	
	