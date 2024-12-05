#idea: combine all predictors and response dataframes into one master df for input to models

#load libraries
library(dplyr)
library(tidyr)
library(ggplot2)

rm(list = ls())
gc()

setwd("/Users/rachel/divdiv")

#-----------------------------------------------------------------------------------------------------------------
#read in data
lat <- read.csv("data/abiotic/mean_min_max_latitude_stats-wide.csv") %>% 
  dplyr::select(-X) %>% mutate(run_name = paste0("bioprj_",link)) %>% 
  dplyr::select(-link)
ecor <- read.csv("data/abiotic/number_of_ecoregions_stats-wide.csv") %>% mutate(species = gsub(" ","_",species))
rangesize <- read.csv("data/abiotic/range_lengths_and_ratios.csv")
methods <- read.csv("data/methodological/methodological_predictors-wide.csv") %>% 
  dplyr::select(run_name, n_samples, mean_raw_read_cnt, read_length, 
                mean_locus_depth, littlem, bigm, n, mean_depth_snp_filtered, 
                n.totalsnps)

biotic <- read.csv("data/biotic/cleaned_numeric_biotic_traits.csv") %>% dplyr::rename("species"="organism_biosamp") %>% 
  mutate(species = gsub(" ","_",species))

popgWM <- read.csv("data/popgen/r80_WM_stats-wide.csv") %>% dplyr::select(-species)
popgsummary <- read.csv("data/popgen/r80_popgen_stats-wide.csv") %>% 
  dplyr::select(-species) %>% dplyr::select(-stacksparams)

taxcolorkey <- read.csv("data/master_tax_color_key.csv") %>% mutate(species = gsub(" ","_",species))

#-----------------------------------------------------------------------------------------------------------------

#drop a few biotic preds we don't want
biotic <- biotic %>% dplyr::select(-PLD_point, -PLD_Min, -PLD_Max)

#make all latitude values positive
lat %>% ggplot() + geom_histogram(aes(x = meanlat.gbif), bins = 30) + scale_x_continuous(limits = c(-90,90))
lat <- lat %>% mutate(meanlat.gbif = abs(meanlat.gbif))
lat %>% ggplot() + geom_histogram(aes(x = meanlat.gbif), bins = 30) + scale_x_continuous(limits = c(-90,90))

#sync some species names
biotic$species[biotic$species == "Exaiptasia_pallida"] = "Exaiptasia_diaphana" 
biotic$species[biotic$species == "Seriola_dorsalis"] = "Seriola_lalandi_dorsalis" 
taxcolorkey$species[taxcolorkey$species == "Seriola_dorsalis"] = "Seriola_lalandi_dorsalis"

#name issues to correct/sync
#Exaiptasia diaphana in gbif ecor
#Exaiptasia diaphana in gbif mean/min/max latitudes
#Exaiptasia pallida in full list

#Seriola dorsalis in full list
#Seriola lalandi dorsalis in gbif ecor
#Seriola lalandi dorsalis in gbif mean/min/max latitudes


#do fake merge first to check that name mismatches are corrected
df <- merge(lat, ecor, by = "species", all = T)
df <- merge(df, rangesize, by = "run_name", all = T)
df <- merge(df, biotic, by = "species", all = T)
df <- merge(df, popgWM, by = "run_name", all = T)
df %>% filter(grepl("Exaiptasia",run_name)) %>% dplyr::select(run_name,species)
df %>% filter(grepl("Exaiptasia",species)) %>% dplyr::select(run_name,species)
df %>% filter(grepl("Seriola",run_name)) %>% dplyr::select(run_name,species)
df %>% filter(grepl("Seriola",species)) %>% dplyr::select(run_name,species)


#do real merge
df <- merge(lat, ecor, by = "species", all = F)
df <- merge(df, rangesize, by = "run_name", all = F)
df <- merge(df, biotic, by = "species", all = F)
df <- merge(df, methods, by = "run_name", all = F)
df.preds <- df
df <- merge(df, popgWM, by = "run_name", all = F)
df <- merge(df, popgsummary, by = "run_name", all = F)
#add colors for tax. groups (that match phy)
df <- merge(df, taxcolorkey, by = "species", all = F)
names(df)

#check if there are any duplicate rows/datasets (should return no rows if no dups)
df %>% group_by(run_name) %>% summarise(n=n()) %>% filter(n > 1)

#change a few species names (AGAIN...), for our sanity
df$species[df$species == "Seriola_lalandi_dorsalis"] = "Seriola_dorsalis"
df$species[df$species == "Halichoerus-grypus-atlantica"] = "Halichoerus_grypus"



# save --------
#write.csv(df, paste0("data/master_df-",Sys.Date(),".csv"), row.names = FALSE)
write.csv(df, "data/master_df.csv", row.names = FALSE)



