#idea: make a map of the world and plot NCBI current winners marine data on it

#https://spatialreference.org - to look up projection codes, like PROJ.4 codes

#load libraries -------
library(sf)      
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)

rm(list=ls())
gc()



# load files for world map --------

load("/Users/rachel/divdiv/data/abiotic/input_and_working/earth_map_objects.RData")

## This will load 5 R objects:
##   lbl.X & lbl.Y are two data.frames that contain labels for graticule lines
##   NE_box is a SpatialPolygonsDataFrame object and represents a bounding box for Earth 
##   NE_countries is a SpatialPolygonsDataFrame object representing countries 
##   NE_graticules is a SpatialLinesDataFrame object that represents 10 dg latitude lines and 20 dg longitude lines
##   NOTE: data downloaded from http://www.naturalearthdata.com/

# get Great Lakes outlines
lakes <- rnaturalearth::ne_download(scale="medium", category = 'physical', type = "lakes", returnclass = "sf")
#subset world lake outlines to just GLs
gls <- c("Lake Superior", "Lake Huron",
         "Lake Michigan", "Lake Erie", "Lake Ontario")
gls <- subset(lakes, name %in% gls)
#how does it look?
ggplot() +
  geom_sf(data = gls)

# get the genetic sample points we want to plot on map ----------

#get master list of datasets we're using
using <- read.csv("/Users/rachel/divdiv/data/master_df.csv")
using <- using %>% dplyr::select(species,link) %>% mutate(link = gsub("bioprj_","",link))

#get lat/longs
datasets <- list.files(path = "/Users/rachel/divdiv/data/abiotic/input_and_working/lat_long_tables_per_dataset",
                       pattern = "lat_long_table", full.names = TRUE)
sites <- data.frame("link"=NA, "run_acc_sra"=NA, "lat"=NA, "long"=NA, "coordinateUncertaintyInMeters"=NA, "origin"=NA)
for ( d in datasets ) {
  out <- read.delim(d)
  sites <- rbind(sites,out)
}
sites <- sites %>% filter(is.na(link)==F)

#keep ones we want
sites <- sites %>% filter(link %in% using$link) %>% rename("lon" = "long")

#add clade groups on for coloring
taxcolorkey <- read.csv("/Users/rachel/divdiv/data/master_tax_color_key.csv")
sites <- merge(sites %>% separate(., link, into = c("garbage","species"), sep = "_", remove = F) %>% dplyr::select(-garbage) %>% mutate(species = gsub("-","_",species)), 
               taxcolorkey %>% mutate(species = gsub(" ","_",species)), 
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

## project long-lat coordinates columns for data frames
lbl.X.prj <- st_as_sf(lbl.X, coords = c("lon", "lat"), crs = PROJ)
lbl.Y.prj <- st_as_sf(lbl.Y, coords = c("lon", "lat"), crs = PROJ)

#project sites
PROJ.pts <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" #EPSG:4326, WGS 84
#first make dataframe into spatial object, using projection points are in
sites.prj <- st_as_sf(sites, coords = c("lon", "lat"), crs = PROJ.pts)
#then transform projection to match other map layers
sites.prj <- st_transform(sites.prj, crs = PROJ)



# plot ----------------------
ggplot() +
  
  ## add projected bounding box, (a blue for ocean = #abd9e9)
  geom_sf(data = NE_box.prj,
          colour = "black", fill = "white", size = .25) +
  
  ## add graticules
  geom_sf(data = NE_graticules.prj,
          linetype = "dotted", colour = "grey70", size = .25) +
  
  ## add projected countries
  geom_sf(data = NE_countries.prj, 
          colour = "gray50", fill = "gray90", size = .25) +

  ## add great lakes
  geom_sf(data = gls.prj, 
          colour = "gray50", fill = "white", size = .25) +
  
  ## add genetic sample points
  geom_sf(data = st_jitter(sites.prj, 5000), aes(fill = taxclade), shape=21, colour="black", size = 4, alpha = 0.5) +
  scale_fill_manual(values = c("#A6CEE3","#33A02C","#FF7F00","#FDBF6F","#FB9A99","#B2DF8A","#E31A1C","#C2C0C0","#6A3D9A","#1F78B4","#CAB2D6")) +

  ## Set empty theme
  theme_void() + # remove the default background, gridlines & default gray color around legend's symbols
  theme(
    legend.position = "none",
    #legend.text = element_text(size = 20),
    #legend.title = element_blank(),
    #plot.margin = unit(c(t=0, r=3, b=0, l=0), unit="cm"), # adjust margins
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent")
  )

ggsave(filename = "/Users/rachel/divdiv/figures/world_map_genetic_pts.pdf", width = 30, height = 15, units = c("cm"))



