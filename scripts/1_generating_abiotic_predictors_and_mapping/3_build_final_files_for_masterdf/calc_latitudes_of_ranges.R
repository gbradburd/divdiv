#idea: calc min, max, and mean latitude for both GBIF occurence records and sampled genetic points

#load libraries
library(dplyr)
library(tidyr)
library(ggplot2)


rm(list = ls())
gc()

indir.gbif = "/Users/rachel/divdiv/data/abiotic/input_and_working/speciesLocDat_Sep08_2021_withFilters-cleaned"
indir.genetic = "/Users/rachel/divdiv/data/abiotic/input_and_working/lat_long_tables_per_dataset"
outdir = "/Users/rachel/divdiv/data/abiotic"



#get lats for each GBIF dataset - range wide -------------
fileslist <- list.files(path = indir.gbif, pattern = ".csv", full.names = TRUE)


df <- data.frame("species"=NA, "meanlat.gbif"=NA, "minlat.gbif"=NA, "maxlat.gbif"=NA)
for (loop.iter in 1:length(fileslist)) {
  
  print(paste0("iteration ", loop.iter, " of ", length(fileslist)))
  file = fileslist[loop.iter]
  out <- read.csv(file)
  
  #get species name
  species = file %>% strsplit(., split = "/") %>% as.data.frame() %>% .[nrow(.),] %>% 
    gsub("GBIFLocations_","",.) %>% gsub(".csv","",.) %>% gsub(" ", "_", .)
  
  #calc lats
  lats <- out %>% mutate(meanlat.gbif = mean(decimalLatitude),
                         minlat.gbif = min(decimalLatitude),
                         maxlat.gbif = max(decimalLatitude),
                         species = species) %>% 
    dplyr::select(species, contains("gbif")) %>% distinct()

  df <- rbind(df, lats)
  rm(file,out,lats)
  
}
df <- df %>% filter(is.na(species)==F)

df %>% ggplot() + geom_histogram(aes(x = meanlat.gbif), binwidth = 5) + coord_flip() + geom_vline(xintercept = 0, colour = "blue", linetype = "dashed") + 
  labs(x = "Mean range latitude",
       y = "Number of datasets")

df.gbif <- df



#get lats for each genetic dataset - sampled points -------------
fileslist <- list.files(path = indir.genetic, pattern = "lat_long_table", full.names = TRUE)


df <- data.frame("link"=NA, "species"=NA, "meanlat.genetic"=NA, "minlat.genetic"=NA, "maxlat.genetic"=NA)
for (loop.iter in 1:length(fileslist)) {
  
  print(paste0("iteration ", loop.iter, " of ", length(fileslist)))
  file = fileslist[loop.iter]
  out <- read.delim(file)
  
  #get species name
  species = file %>% strsplit(., split = "/") %>% as.data.frame() %>% .[nrow(.),] %>% 
    gsub("lat_long_table-bioprj_","",.) %>% gsub(".txt","",.) %>% strsplit(., split = "_") %>% 
    as.data.frame() %>% .[nrow(.),] %>% gsub("-","_",.)
  
  #calc lats
  lats <- out %>% mutate(meanlat.genetic = mean(lat),
                         minlat.genetic = min(lat),
                         maxlat.genetic = max(lat),
                         species = species) %>% 
    dplyr::select(link, species, contains("genetic")) %>% distinct()
  
  df <- rbind(df, lats)
  rm(file,out,lats)
  
}
df <- df %>% filter(is.na(species)==F)

df %>% ggplot() + geom_histogram(aes(x = meanlat.genetic), binwidth = 5) + coord_flip() + geom_vline(xintercept = 0, colour = "blue", linetype = "dashed") + 
  labs(x = "Mean range latitude",
       y = "Number of datasets")

df.genetic <- df


#join GBIF and genetic
head(df.gbif)
head(df.genetic)

df <- merge(df.genetic, df.gbif, by = "species", all.x = T) %>% dplyr::select(link, everything())

#save
write.csv(df, paste0(outdir, "/mean_min_max_latitude_stats-wide.csv"))






