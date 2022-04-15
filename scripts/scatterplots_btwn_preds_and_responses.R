#idea: look at scatterplots between our model predictors and responses to get a sense for expected results

#load libraries
library(dplyr)
library(tidyr)
library(ggplot2)

rm(list = ls())
gc()

#read in data
df <- read.csv("data/master_df.csv")

#for now, drop one outlier species
df <- df %>% filter(species != "Pocillopora_damicornis")

#get preds we're using
preds <- c("meanlat.gbif","n_ECOREGIONS.all","maxgbif.sea_km","Body_Size",
           "Fecundity_EggSize","Generational_Structure",
           "ReturnToSpawningGround","Spawning_mode","Larval_feeding",
           "PLD_point2","isPlanktonic_atanypoint")

df <- df %>% dplyr::select(species, all_of(preds), s, nbhd, inDeme)
df <- df %>% pivot_longer(., names_to = "predictor", values_to = "value", cols = 2:12)

#plot
#response = s
df %>% ggplot(aes(x = value, y = s)) + 
  geom_point() +
  geom_smooth(formula = y~x, method = "lm") +
  geom_text(aes(label = species), size = 1, vjust = 0, nudge_y = 0.00025) +
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  labs(x = "Predictor variable") +
  facet_wrap(~ predictor, scales = "free")

ggsave(paste0("figures/scatterplot-s_vs_preds.pdf") , width = 280, height = 216, dpi = 600, units = c("mm"))


#response = nbhd
df %>% ggplot(aes(x = value, y = nbhd)) + 
  geom_point() +
  geom_smooth(formula = y~x, method = "lm") +
  geom_text(aes(label = species), size = 1, vjust = 0, nudge_y = 3) +
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  labs(x = "Predictor variable") +
  facet_wrap(~ predictor, scales = "free")

ggsave(paste0("figures/scatterplot-nbhd_vs_preds.pdf"), width = 280, height = 216, dpi = 600, units = c("mm"))


#response = inDeme
df %>% ggplot(aes(x = value, y = inDeme)) + 
  geom_point() +
  geom_smooth(formula = y~x, method = "lm") +
  geom_text(aes(label = species), size = 1, vjust = 0, nudge_y = 0.02) +
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  labs(x = "Predictor variable") +
  facet_wrap(~ predictor, scales = "free")

ggsave(paste0("figures/scatterplot-inDeme_vs_preds.pdf"), width = 280, height = 216, dpi = 600, units = c("mm"))




#graveyard ------------

#for making one pdf/plot per pred/repsonse combo
for ( p in preds ){
  
  #response = s
  print(paste0("p is ", p))
  df %>% ggplot(aes(x = .data[[p]], y = s)) + geom_point() +
    geom_smooth(formula = y~x, method = "lm") +
    theme_bw() + 
    theme(panel.grid = element_blank())

  ggsave(paste0("figures/scatterplot-s_vs_",p,".pdf") , width = 80, height = 80, dpi = 600, units = c("mm"))
  
  #response = nbhd
  print(paste0("p is ", p))
  df %>% ggplot(aes(x = .data[[p]], y = nbhd)) + geom_point() +
    geom_smooth(formula = y~x, method = "lm") +
    theme_bw() + 
    theme(panel.grid = element_blank())
  
  ggsave(paste0("figures/scatterplot-nbhd_vs_",p,".pdf") , width = 80, height = 80, dpi = 600, units = c("mm"))
  
  #response = inDeme
  print(paste0("p is ", p))
  df %>% ggplot(aes(x = .data[[p]], y = inDeme)) + geom_point() +
    geom_smooth(formula = y~x, method = "lm") +
    theme_bw() + 
    theme(panel.grid = element_blank())
  
  ggsave(paste0("figures/scatterplot-inDeme_vs_",p,".pdf") , width = 80, height = 80, dpi = 600, units = c("mm"))
  
}
