#idea: save out final clean .csv of biotic traits used in end to include in suppmat

#load libraries
library(dplyr)

rm(list=ls())
gc()

#define dir paths
indir = "/Users/rachel/divdiv"
outdir = "/Users/rachel/divdiv/data/biotic"

#get cleaned final trait data
df <- read.csv(paste0(workdir,"/data/biotic/inputs_and_working/marinerds_traits_01-30-2024.csv"))

#get final master df used in analyses
master <- read.csv(paste0(workdir,"/data/master_df.csv"))

#keep just species and traits 

# START HERE - need to get match up between Fig 4 names and master df (waiting for answer from Gideon)
df <- df %>% filter(species %in% master$species) %>% 
  dplyr::select(species, X, X, Spawning_mode, isPlanktonic_atanypoint)

#save
write.csv(df, paste0(outdir, "/final_biotic_traits.csv"), row.names = FALSE)


