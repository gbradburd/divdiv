#idea: drop spatial outliers and points in invasives ranges from downloaded GBIF points
#built the list of points to drop by looking visually at a map of the GBIF points for each species
#and comparing this to range maps/natural history knowledge and circling on the map which points to keep or drop
#then found IDs for points to drop and compiled them all into one master "to drop" list


#load libraries
library(dplyr)
library(tidyr)

rm(list = ls())
gc()

#define where original files of lat/long GBIF points live
indir = "/Users/rachel/divdiv/data/abiotic/input_and_working/speciesLocDat_Sep08_2021_withFilters"

#define where to saved new lat/long GBIF points with invasive and spatial outliers dropped
outdir = "/Users/rachel/divdiv/data/abiotic/input_and_working/speciesLocDat_Sep08_2021_withFilters-cleaned"

#optional - define where to save figures that visualize which points were dropped for each species
figdir = "/Users/rachel/divdiv/figures/maps_to_check_cleaning"

#optional - define where country outlines Robj lives, for making visualization maps
countries = "/Users/rachel/divdiv/data/abiotic/input_and_working/countries.outlines.Robj"

#get master list of points to drop
todrop <- read.csv("master_list_GBIF_points_to_drop.csv")



#read in original GBIF points that we downloaded and saved for each species
#filter out any points in our to drop list
#save new, cleaned df's

#get list of all GBIF files to process
file_list <- list.files(indir, pattern="GBIFLocations", full.names = TRUE)

#drop points from each df and save
for (loop.iter in 1:length(file_list)) {
  
  print(paste0("loop ", loop.iter, " of ", length(file_list)))
  file <- file_list[loop.iter]
  speciesX = file %>% stringr::str_split(., "/") %>% as.data.frame() %>% .[nrow(.),] %>% gsub("GBIFLocations_","",.) %>% gsub(".csv","",.)
  file <- read.csv(file) %>% mutate(sampid = paste0("gbif",X)) %>% dplyr::select(-X) %>% dplyr::select(sampid, everything())
  
  todrop.speciesX <- todrop %>% filter(species == speciesX)
  
  if (length(todrop.speciesX > 0)) {
  out <- file %>% filter(!sampid %in% todrop.speciesX$sampid)
  print(paste0("dropped ", nrow(file)-nrow(out), " points"))
  } else {
  print("no points to drop")
  out <- file
  }
  
  write.csv(out, paste0(outdir,"/GBIFLocations_",speciesX,".csv"), 
            row.names = FALSE)

}



#check that dropping points went correctly by mapping ---------
file_list <- list.files(outdir, pattern="GBIFLocations", full.names = TRUE)
df <- data.frame("sampid"=NA, "decimalLongitude"=NA, "decimalLatitude"=NA, "establishmentMeans"=NA, "occurrenceStatus"=NA, "status_x"=NA, "species"=NA)
for (loop.iter in 1:length(file_list)) {
  print(paste0("loop ", loop.iter, " of ", length(file_list)))
  file <- file_list[loop.iter]
  species = file %>% stringr::str_split(., "/") %>% as.data.frame() %>% .[nrow(.),] %>% gsub("GBIFLocations_","",.) %>% gsub(".csv","",.)
  file <- read.csv(file) %>% mutate(status_x = "a cleaned dataset", species = species)
  df <- rbind(df, file)
}
df.cleaned <- df %>% filter(is.na(status_x)==FALSE)


file_list <- list.files(indir, pattern="GBIFLocations", full.names = TRUE)
df <- data.frame("sampid"=NA, "decimalLongitude"=NA, "decimalLatitude"=NA, "establishmentMeans"=NA, "occurrenceStatus"=NA, "status_y"=NA, "species"=NA)
for (loop.iter in 1:length(file_list)) {
  print(paste0("loop ", loop.iter, " of ", length(file_list)))
  file <- file_list[loop.iter]
  species = file %>% stringr::str_split(., "/") %>% as.data.frame() %>% .[nrow(.),] %>% gsub("GBIFLocations_","",.) %>% gsub(".csv","",.)
  file <- read.csv(file) %>% mutate(status_y = "an og dataset") %>% mutate(sampid = paste0("gbif",X)) %>% mutate(species = species) %>%
    dplyr::select(-X) %>% dplyr::select(sampid, everything())
  df <- rbind(df, file)
}
df.og <- df %>% filter(is.na(status_y)==FALSE)

df.cleaned <- df.cleaned %>% dplyr::select(-establishmentMeans, -occurrenceStatus, -decimalLongitude, -decimalLatitude) %>% filter(species %in% todrop$species)
df.og <- df.og %>% dplyr::select(-establishmentMeans, -occurrenceStatus) %>% filter(species %in% todrop$species)

nrow(df.og) - nrow(df.cleaned) 
df <- merge(df.cleaned, df.og, by = c("sampid","species"), all = T)
df %>% group_by(status_x, status_y) %>% summarise(n=n())
df <- df %>% mutate(status = ifelse(is.na(status_x)==T,"dropped","kept"))
df %>% group_by(status_x, status_y, status) %>% summarise(n=n())


load(countries)
for ( sp in unique(df$species) ) {
  
  df.spX <- df %>% filter(species == sp)
  df.spX %>% ggplot() + 
    geom_polygon(data = countries, aes(long, lat, group = group), colour = "gray50", fill = "white", size = 0.25) +
    geom_point(aes(x=decimalLongitude, y=decimalLatitude, colour=status, shape=status)) +
    scale_x_continuous(breaks = seq(-180,180,20)) +
    scale_y_continuous(breaks = seq(-90,90,20)) +
    scale_shape_manual(values = c(21, 22)) + 
    labs(title = sp) +
    coord_fixed()
  ggsave(paste0(figdir,"/",sp,".png"), width = 12, height = 8, units = c("in"))
  
}
