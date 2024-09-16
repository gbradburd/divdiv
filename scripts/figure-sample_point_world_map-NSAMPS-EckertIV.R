#idea: make a map of the world and plot heat map of N samples on it

#https://spatialreference.org - to look up projection codes, like PROJ.4 codes

#load libraries -------
library(sf)      
library(ggplot2)
library(dplyr)
library(tidyr)

rm(list=ls())
gc()


setwd("/Users/rachel/divdiv/")
sampkeypath = "/Users/rachel/divdiv/data/all_samplenamekeys"
wmpath = "/Users/rachel/divdiv/data/popgen/input_and_working/ALL_r80_gendiv_data/"


#load files for world map --------

load("data/abiotic/input_and_working/earth_map_objects.RData")

## This will load 5 R objects:
##   xbl.X & lbl.Y are two data.frames that contain labels for graticule lines
##   NE_box is a SpatialPolygonsDataFrame object and represents a bounding box for Earth 
##   NE_countries is a SpatialPolygonsDataFrame object representing countries 
##   NE_graticules is a SpatialLinesDataFrame object that represents 10 dg latitude lines and 20 dg longitude lines
##   NOTE: data downloaded from http://www.naturalearthdata.com/

#also get coastline outline
NE_coastlines <- rnaturalearth::ne_coastline(scale = "medium", returnclass = 'sf')
NE_countries <- rnaturalearth::ne_countries(scale = "medium", returnclass = 'sf')

#also get Great Lakes outlines
lakes <- rnaturalearth::ne_download(scale="medium", category = 'physical', type = "lakes", returnclass = "sf")
#subset world lake outlines to just GLs
gls <- c("Lake Superior", "Lake Huron",
         "Lake Michigan", "Lake Erie", "Lake Ontario")
gls <- subset(lakes, name %in% gls)



# get the genetic sample points we want to plot on map ----------

#get all lat/longs
datasets <- list.files(path = "data/abiotic/input_and_working/lat_long_tables_per_dataset",
                       pattern = "lat_long_table", full.names = TRUE)
sites <- data.frame("link"=NA, "run_acc_sra"=NA, "lat"=NA, "long"=NA, "coordinateUncertaintyInMeters"=NA, "origin"=NA)
for ( d in datasets ) {
  out <- read.delim(d)
  sites <- rbind(sites,out)
}
sites <- sites %>% filter(is.na(link)==F) %>% dplyr::rename("lon" = "long")


#filter to just those that made it thru popgen
wmlist <- list.files(path = wmpath, 
                     pattern="initPars", 
                     full.names = TRUE)

sites.kept <- data.frame(sampid_assigned_for_bioinf = NA, run_name = NA, run_acc_sra = NA)
for (i in 1:length(wmlist)) {
  
  wmFile <- wmlist[i]
  dataset = wmFile %>% strsplit(., "/") %>% as.data.frame() %>% .[nrow(.),] %>% 
    gsub("WMfitwishart-","",.) %>% gsub("_initPars.Robj","",.)
  run_name = dataset %>% strsplit(., "_") %>% unlist() %>% .[1:3] %>% 
    paste(., sep="", collapse="_")
  
  load(wmFile)
  if (dim(initPars$parHom)[1] != dim(initPars$parHom)[2]) {"error - matrix dims are not the same"}
  df.i = rownames(initPars$parHom) %>% as.data.frame() %>% dplyr::rename("sampid_assigned_for_bioinf"=".") %>% 
    mutate(run_name = run_name)
  
  sampkeyFile <- paste0(sampkeypath,"/samplenamekey-",run_name,".txt")
  sampkey <- read.delim(sampkeyFile)
  df.i <- merge(df.i, 
                sampkey %>% dplyr::select(sampid_assigned_for_bioinf, run_acc_sra), 
                by = "sampid_assigned_for_bioinf", all.x = TRUE)
  sites.kept <- rbind(sites.kept, df.i)
  
}
sites.kept <- sites.kept %>% filter(is.na(sampid_assigned_for_bioinf)==F)

#are run_acc_sra's all unique in two dfs about to use/filter by?
#should return 0 if so
nrow(sites) - sites %>% dplyr::select(run_acc_sra) %>% distinct() %>% nrow()
nrow(sites.kept) - sites.kept %>% dplyr::select(run_acc_sra) %>% distinct() %>% nrow()

sites <- sites %>% filter(run_acc_sra %in% sites.kept$run_acc_sra)

#and filter by datasets in master list of datasets we're using as final check
using <- read.csv("data/master_df.csv") %>% mutate(link = gsub("bioprj_","",run_name))
using <- using %>% dplyr::select(species,link)
sites <- sites %>% filter(link %in% using$link)



#add clade groups on for coloring
taxcolorkey <- read.csv("data/master_tax_color_key.csv") %>% mutate(species = gsub(" ","_",species))
#recode to match phy fig
taxcolorkey <- read.csv("/Users/rachel/divdiv/data/master_tax_color_key.csv") %>% mutate(species = gsub(" ","_",species))
taxcolorkey$simplerclades <- taxcolorkey$taxclade
taxcolorkey$simplerclades[taxcolorkey$taxclade == "sauropsida"] = "Birds"
taxcolorkey$simplerclades[taxcolorkey$taxclade == "chondrichthyes"] = "Cartilaginous fishes"
taxcolorkey$simplerclades[taxcolorkey$taxclade == "not assigned"] = "Other"
taxcolorkey$simplerclades[taxcolorkey$taxclade == "ochrophyta"] = "Other"
taxcolorkey$simplerclades[taxcolorkey$species == "Caretta_caretta"] = "Other"
#merge
sites <- sites %>%
  separate(., link, into = c("garbage","species"), sep = "_", remove = F) %>% 
  dplyr::select(-garbage) %>% mutate(species = gsub("-","_",species))
sites <- merge(sites, 
               taxcolorkey, 
               by = "species", all.x = T)



# set up projections for map and our points --------------
## give the PORJ.4 string for Eckert IV projection
PROJ <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" 
## or use the short form "+proj=eck4"

## convert shapefiles to sf objects and ensure all in same projection
NE_countries <- sf::st_as_sf(NE_countries)
NE_countries.prj <- sf::st_transform(NE_countries, crs = PROJ)
NE_graticules <- sf::st_as_sf(NE_graticules)
NE_graticules.prj <- sf::st_transform(NE_graticules, crs = PROJ)
NE_box <- sf::st_as_sf(NE_box)
NE_box.prj <- sf::st_transform(NE_box, crs = PROJ)
gls.prj <- sf::st_transform(gls, crs = PROJ)
NE_coastlines.prj <- sf::st_transform(NE_coastlines, crs = PROJ)

#project sites
PROJ.pts <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" #EPSG:4326, WGS 84
#first make dataframe into spatial object, using projection points are in
sites.prj <- sf::st_as_sf(sites, coords = c("lon", "lat"), crs = PROJ.pts)
#then transform projection to match other map layers
sites.prj <- sf::st_transform(sites.prj, crs = PROJ)



# do projecting and heat/raster layer making -----------
#make heatmap raster - N number of points in each box (large res = fewer/bigger squares)
rast <- raster::raster(res = 5)
raster::extent(rast) <- raster::extent(sites.prj)
rast
raster.Nindivs <- raster::rasterize(sites.prj, rast, fun = "count")
raster.Nindivs <- raster::rasterToPoints(raster.Nindivs) %>% as.data.frame()



# make heatmap - cells colored by N individuals - white ocean --------------
ggplot() +
  
  ## add projected bounding box, (a blue for ocean = #abd9e9)
  geom_sf(data = NE_box.prj,
          colour = "black", fill = "white", size = .25) +
  
  ## add graticules
  geom_sf(data = NE_graticules.prj,
          linetype = "dotted", colour = "grey70", size = .25) +
  
  ## add projected countries
  geom_sf(data = NE_countries.prj, 
          colour = "gray90", fill = "gray90", size = .25) +
  
  # add land outlines
  geom_sf(data = NE_coastlines.prj, 
          colour = "gray50", linewidth = .25) +
  
  ## add great lakes
  geom_sf(data = gls.prj, 
          colour = "gray50", fill = "white", size = .25) +
  
  ## add N indivs raster
  geom_tile(data = raster.Nindivs, aes(x = x, y = y, fill = ID)) +
  
  viridis::scale_fill_viridis(name = "N. samples", trans = "log",
                              breaks = c(0,1,10,100,1000), labels = c(0,1,10,100,1000),
                              option="mako") +
  
  ## Set empty theme
  theme_void() + # remove the default background, gridlines & default gray color around legend's symbols
  theme(legend.title = element_text(colour="black", size=10, face="bold"), # adjust legend title
        legend.position = c(1.035, 0.2), # relative position of legend
        plot.margin = unit(c(t=0, r=3, b=0, l=0), unit="cm"), # adjust margins
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA))

ggsave(paste("figures/world_map_genetic_pts-NSAMPS-EckertIV.pdf",sep="") , width = 30, height = 15, units = c("cm"))
ggsave(paste("figures/world_map_genetic_pts-NSAMPS-EckertIV.png",sep="") , width = 30, height = 15, units = c("cm"))


