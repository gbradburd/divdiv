#idea: combine all predictors and response dataframes into one master df for input to models

#load libraries
library(dplyr)
library(tidyr)

rm(list = ls())
gc()

#read in data
lat <- read.csv("data/abiotic/mean_min_max_latitude_stats-wide.csv") %>% dplyr::select(-X) %>% mutate(link = paste0("bioprj_",link))
ecor <- read.csv("data/abiotic/number_of_ecoregions_stats-wide.csv") %>% mutate(species = gsub(" ","_",species))
rangesize <- read.csv("data/abiotic/range_lengths_and_ratios.csv") %>% dplyr::rename("link"="run_name")

biotic <- read.csv("data/biotic/cleaned_numeric_biotic_traits.csv") %>% dplyr::rename("species"="organism_biosamp") %>% 
  mutate(species = gsub(" ","_",species))

popg <- read.csv("data/popgen/r80_popgen_WM_stats-wide.csv") %>% dplyr::rename("link"="run_name") %>% dplyr::select(-species)


#keep just traits in hypoth bingo
biotic <- biotic %>% dplyr::select(species, Body_Size, Fecundity_EggSize, Generational_Structure, 
                                   ReturnToSpawningGround, Spawning_mode, Larval_feeding,
                                   PLD_point2, isPlanktonic_atanypoint)


#name issues to correct/sync
#Exaiptasia diaphana in gbif ecor
#Exaiptasia diaphana in gbif mean/min/max latitudes
#Exaiptasia pallida in full list

#Seriola dorsalis in full list
#Seriola lalandi dorsalis in gbif ecor
#Seriola lalandi dorsalis in gbif mean/min/max latitudes


#merge
df <- merge(lat, ecor, by = "species", all = F)
df <- merge(df, rangesize, by = "link", all = F)
df <- merge(df, biotic, by = "species", all = F)
df <- merge(df, popg, by = "link", all = F)

names(df)


#save
#write.csv(df, paste0("data/master_df-",Sys.Date(),".csv"), row.names = FALSE)
write.csv(df, "data/master_df.csv", row.names = FALSE)


