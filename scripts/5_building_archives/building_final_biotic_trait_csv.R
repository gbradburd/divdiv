#idea: save out final clean .csv of biotic traits used in end to include in suppmat

#load libraries
library(dplyr)

rm(list=ls())
gc()

#define dir paths
indir = "/Users/rachel/divdiv"
outdir = "/Users/rachel/divdiv/data/biotic"

#get cleaned final trait data
trts.og <- read.csv(paste0(indir,"/data/biotic/inputs_and_working/marinerds_traits_01-30-2024.csv"))

#get final master df used in analyses
master <- read.csv(paste0(indir,"/data/master_df.csv")) %>% 
  mutate(organism_biosamp = gsub("_"," ",species))
master$organism_biosamp[master$organism_biosamp == "Exaiptasia diaphana"] = "Exaiptasia pallida"

#keep just species and traits 

df <- trts.og %>% filter(organism_biosamp %in% master$organism_biosamp) %>% 
  dplyr::select(organism_biosamp,
                Generational_Structure,ReturnToSpawningGround,
                Spawning_mode,Larval_feeding,isPlanktonic_atanypoint,
                isBenthic,
                PLD_point2,
                Body_Size, Fecundity_EggSize)

#save
write.csv(df, paste0(outdir, "/final_biotic_traits.csv"), row.names = FALSE)

#used combo of an Rmarkdown script and quite a bit of manual processing to add 
#numerical citation value columns to final_biotic_traits.csv and 
#generate a full citation text doc

#notes:

#name matchup between trait work (and master df) and pred names in final paper
"Generational_Structure" -> "iteroparity"
"ReturnToSpawningGround" -> "philopatry"
"Spawning_mode" -> "spawning mode"
"Larval_feeding" -> "planktotrophy"
"isPlanktonic_atanypoint" -> "planktonicity"
"isBenthic" -> "benthicity"
"PLD_point2" -> "PLD"
"BodySize" -> "body size"
"Fecundity_EggSize" -> "egg size"

#how we recoded some variables for modeling
#for dispersal related traits, larger values mean higher dispersal
isPlanktonic_atanypoint
"TRUE" = 1
"FALSE" = 0

ReturnToSpawningGround
"TRUE" = 1
"FALSE" = 0
NA

Spawning_mode
"I" = 0
"N" = 1
"F" = 2

Larval_feeding
"N" = 0
"L" = 1
"P" = 2

Generational_Structure
"I" = 0
"S" = 1

isBenthic
"N" = 0
"S" = 1
"A" = 2

str(df)

for ( trait in names(df) ) {
  print(df %>% group_by(!!as.symbol(trait)) %>% summarise(n=n()))
}








