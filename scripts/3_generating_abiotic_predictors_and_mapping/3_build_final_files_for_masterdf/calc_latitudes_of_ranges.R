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


# END FORMAL CODE


# side bar - effect of using original GBIF locations rather than relocated locations for "land points" ---------------

#get lats for each GBIF dataset AFTER MOVING POINTS - range wide -------------
fileslist <- list.files(path = indir.gbif.post, pattern = "gbif", full.names = TRUE)


df.gbif.post <- data.frame("species"=NA, "meanlat.gbif.post"=NA, "minlat.gbif.post"=NA, "maxlat.gbif.post"=NA)
for (loop.iter in 1:length(fileslist)) {
  
  print(paste0("iteration ", loop.iter, " of ", length(fileslist)))
  file = fileslist[loop.iter]
  out <- read.csv(file)
  
  #get species name
  species = file %>% strsplit(., split = "/") %>% as.data.frame() %>% .[nrow(.),] %>%
    as.data.frame() %>% separate(., ".", into = c("garbage","garbage2","sp"), sep = "_") %>%
    dplyr::select(sp) %>% 
    gsub(".csv","",.) %>% gsub("-", "_", .)
  
  #calc lats
  lats <- out %>% 
    mutate(end.lat = ifelse("end.lat" %in% names(.), end.lat, NA)) %>% 
    mutate(end.lat = ifelse(is.na(end.lat)==T, y.orig, end.lat)) %>% 
    mutate(meanlat.gbif.post = mean(end.lat),
           minlat.gbif.post = min(end.lat),
           maxlat.gbif.post = max(end.lat),
           species = species) %>% 
    dplyr::select(species, contains("gbif")) %>% distinct()
  
  df.gbif.post <- rbind(df.gbif.post, lats)
  rm(file,out,lats)
  
}
df.gbif.post <- df.gbif.post %>% filter(is.na(species)==F)

df.gbif.post %>% ggplot() + geom_histogram(aes(x = meanlat.gbif.post), binwidth = 5) + coord_flip() + geom_vline(xintercept = 0, colour = "blue", linetype = "dashed") + 
  labs(x = "Mean range latitude",
       y = "Number of datasets")

#join onto existing df
head(df.gbif)
head(df.genetic)
head(df)
head(df.gbif.post)
df <- merge(df, df.gbif.post, by = "species", all.x = T) %>% dplyr::select(link, everything())
head(df)

#correlation
df %>% ggplot() + geom_point(aes(x=meanlat.gbif, y=meanlat.gbif.post)) + 
  geom_abline(slope = 1, intercept = 0, colour = "blue")

cor(df$meanlat.gbif, df$meanlat.gbif.post, method = "pearson", use = "complete.obs")


# side bar - effect of not combining genetic and GBIF to calc range lat -----------

fileslist <- list.files(path = indir.genetic, pattern = "lat_long_table", full.names = TRUE)


df <- data.frame("link"=NA, "species"=NA, 
                 "meanlat.both"=NA, "minlat.both"=NA, "maxlat.both"=NA,
                 "meanlat.gbif"=NA, "minlat.gbif"=NA, "maxlat.gbif"=NA)
for (loop.iter in 1:length(fileslist)) {
  
  print(paste0("iteration ", loop.iter, " of ", length(fileslist)))
  file = fileslist[loop.iter]
  
  #get species name
  species = file %>% strsplit(., split = "/") %>% as.data.frame() %>% .[nrow(.),] %>% 
    gsub("lat_long_table-bioprj_","",.) %>% gsub(".txt","",.) %>% strsplit(., split = "_") %>% 
    as.data.frame() %>% .[nrow(.),] %>% gsub("-","_",.)
  
  #get genetic data
  genetic <- read.delim(file)
  link = genetic$link[1]
  genetic <- genetic %>% dplyr::select(lat, long) %>% mutate(species = species) %>% 
    mutate(type = "genetic")
  
  #get GBIF data
  species.gbif <- species %>% gsub("_", " ", .)
  gbif <- list.files(path = indir.gbif, pattern = species.gbif, full.names = TRUE)
  gbif <- read.csv(gbif) %>% rename("long"="decimalLongitude", "lat"="decimalLatitude") %>% 
    dplyr::select(lat, long) %>% mutate(species = species) %>% mutate(type = "gbif")
  gbiflats <- gbif %>% mutate(meanlat.gbif = mean(lat),
                              minlat.gbif = min(lat),
                              maxlat.gbif = max(lat)) %>% 
    dplyr::select(meanlat.gbif, minlat.gbif, maxlat.gbif) %>% distinct()
  
  #join
  out <- rbind(genetic, gbif)
  
  out <- out %>% mutate(link = link,
                        species = species,
                        meanlat.both = mean(lat),
                        minlat.both = min(lat),
                        maxlat.both = max(lat)) %>%
    dplyr::select(link,species,meanlat.both,minlat.both,maxlat.both) %>% 
    distinct()
  out <- cbind(out, gbiflats)
  
  df <- rbind(df, out)
  rm(file,species,genetic,species.gbif,gbif,gbiflats,out)
  
}
df <- df %>% filter(is.na(species)==F)

df %>% ggplot() + geom_histogram(aes(x = meanlat.gbif), binwidth = 5) + coord_flip() + geom_vline(xintercept = 0, colour = "blue", linetype = "dashed") + 
  labs(x = "Mean range latitude",
       y = "Number of datasets")


#plot

df %>% ggplot() + geom_point(aes(x=meanlat.both, y=meanlat.gbif)) + 
  geom_abline(slope = 1, intercept = 0, colour = "blue")

cor(df$meanlat.both, df$meanlat.gbif, method = "pearson", use = "complete.obs")

df %>% mutate(mid = abs(maxlat.gbif-minlat.gbif)/2) %>% 
  ggplot() + geom_point(aes(x = abs(meanlat.gbif), y=mid))


