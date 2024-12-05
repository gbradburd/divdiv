#idea: make a map of the world and plot genetic samples on it, colored by phy aka taxon

#load libraries -------
library(ggplot2)
library(ncdf4)
library(dplyr)
library(tidyr)

rm(list=ls())
gc()



setwd("/Users/rachel/divdiv/")

#load Spilhaus projection functions
source("scripts/5_building_figures/spilhaus_functions.R") 



# get spatial data -----------------
#going to convert everything into spatial vects so we can use terra funcs

#define a simple lat/long projection so can ensure all in same proj to start 
PROJ.latlon <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

#get land
NE_countries <- rnaturalearth::ne_countries(scale = "small", returnclass = 'sf') %>% 
  sf::st_transform(., crs = PROJ.latlon) %>% terra::vect()

#get world coastlines
NE_coastlines <- rnaturalearth::ne_coastline(scale = "small", returnclass = 'sf') %>% 
  sf::st_transform(., crs = PROJ.latlon) %>% terra::vect()

#get great lakes
lakes <- rnaturalearth::ne_download(scale="medium", category = 'physical', type = "lakes", returnclass = "sf")
gls <- c("Lake Superior", "Lake Huron",
         "Lake Michigan", "Lake Erie", "Lake Ontario")
gls <- subset(lakes, name %in% gls) %>% 
  sf::st_transform(., crs = PROJ.latlon) %>% terra::vect()

#erase Great Lakes from land
land <- terra::erase(NE_countries, gls)

#build some gridlines (custom func)
grid <- gen.gridlines(c(-180,-90),c(180,90),c(36,18),c(3000,1800)) %>% 
  terra::vect(., "lines")



# project spatial data into spilhaus ----------------

#preprocess (split regions/lines along seam, etc.) and pull out geom info
#width controls the width of the seam (ensures split polygons end up on correct side)
#might need to be tweaked in some cases
#found that as low as 100 m worked for buffered V of data (that extends beyond coasts)!
#width of 2.5km aka width=2500 needed for V of ecor's that stop at coastline
coast.coords <- preproc.split(NE_coastlines, NA.radius=1066000) %>% terra::geom() %>% as.data.frame()
grid.coords <- preproc.split(grid, NA.radius=1066000) %>% terra::geom() %>% as.data.frame()
land.coords <- preproc.split(land, NA.radius=1066000) %>% terra::geom() %>% as.data.frame()

#build new cols of coords in spilhaus proj
coast.coords <- from_lonlat_to_spilhaus_xy(coast.coords$x, coast.coords$y) %>% cbind(., coast.coords)
grid.coords <- from_lonlat_to_spilhaus_xy(grid.coords$x, grid.coords$y) %>% cbind(., grid.coords)
land.coords <- from_lonlat_to_spilhaus_xy(land.coords$x, land.coords$y) %>% cbind(., land.coords)

# convert layers back to SpatVector (using spilhaus coords)
coast.prj <- as.matrix(coast.coords[,c("geom","part","spilhaus_x","spilhaus_y","hole")]) %>% 
  terra::vect(., "lines", atts=terra::values(NE_coastlines))
grid.prj <- as.matrix(grid.coords[,c("geom","part","spilhaus_x","spilhaus_y","hole")]) %>% 
  terra::vect(., "lines", atts=terra::values(grid))
land.prj <- as.matrix(land.coords[,c("geom","part","spilhaus_x","spilhaus_y","hole")]) %>% 
  terra::vect(., "polygons", atts=terra::values(land))



# fix the map corners for polygon layers ----------------

#For now, patching corners erases any visible Pacific seams...
#(Had to buffer lands.prj with 0 width because it was forming invalid geoms after reprojection)
#Note that patch.corners() also takes NA.radius, SA.radius, etc. arguments now
# these radii should be slightly bigger than whatever was used for preproc.split()!
# (you can also set a radius to 0 to only patch particular corners)

land.prj <- patch.corners(terra::buffer(land.prj,0), 
                          NA.radius=1066000+10,
                          SA.radius=1e6+10,
                          Asia.radius=1e6+10) %>% 
  terra::aggregate()



# erase some stuff --------------------
#create gridlines just for ocean (aka erase over land)
grid.ocean <- terra::erase(grid.prj, land.prj)



# expand extent of map ---------------
#save og not expanded Vs first
coast.prj.og <- coast.prj
land.prj.og <- land.prj
grid.ocean.og <- grid.ocean

#amount varies from 0 to 1 and controls how far the edges are expanded out
# (in terms of fraction of the base map's dimensions)
#exp.fudge.fac buffers any polygons to ensure shapes aggregate across map borders
coast.prj <- expand.map(coast.prj, amount=0.05, exp.fudge.fac=10)
land.prj <- expand.map(land.prj, amount=0.05, exp.fudge.fac=10)
grid.ocean <- expand.map(grid.ocean, amount=0.05, exp.fudge.fac=10)

ggplot() + 
  tidyterra::geom_spatvector(data=grid.ocean.og, linewidth = 0.1) +
  tidyterra::geom_spatvector(data=land.prj.og, fill = "forestgreen") +
  tidyterra::geom_spatvector(data=coast.prj.og, colour = "blue")

ggplot() + 
  tidyterra::geom_spatvector(data=grid.ocean, linewidth = 0.1) +
  tidyterra::geom_spatvector(data=land.prj, fill = "forestgreen") +
  tidyterra::geom_spatvector(data=coast.prj, colour = "blue")



# get the genetic sample points ----------
sites <- read.csv("data/final_genetic_latlongs.csv", header = T)

#add clade groups on for coloring
taxcolorkey <- read.csv("/Users/rachel/divdiv/data/master_tax_color_key.csv") %>% mutate(species = gsub(" ","_",species))
sites <- merge(sites, taxcolorkey, 
               by = "species", all.x = T)

sites <- sites %>% mutate(pt_ID = paste0("pt_",1:nrow(.)))

#project genetic lat longs to spilhaus
sites.prj <- as.data.frame(from_lonlat_to_spilhaus_xy(sites$lon, sites$lat))
sites.prj <- cbind(sites.prj, sites)



# make maps ! -------------------------

# option 1 -------------

ggplot() +
  
  # color oceans
  geom_rect(aes(xmin = -11825474-1182547, xmax = 11825474+1182547, ymin = -11825474-1182547, ymax = 11825474+1182547), fill = "skyblue", alpha = 0.3) +
  
  # add graticules over oceans
  tidyterra::geom_spatvector(data=grid.ocean, colour = "gray80", alpha = 0.5, linewidth = 0.25) +
  
  # add land
  tidyterra::geom_spatvector(data = land.prj, fill = "gray95", colour = "gray45", alpha = 1) +
  
  ## add genetic points
  geom_point(data = sites.prj, aes(x=spilhaus_x, y=spilhaus_y, fill=taxclade), 
             #position = position_jitter(width = 1000000, height = 1000000),
             shape = 21, colour="black", size = 4, alpha = 0.5) +
  
  scale_fill_manual(values = c("#1F78B4","#A6CEE3","#33A02C","#FF7F00","#FDBF6F","#FB9A99","#B2DF8A","#E31A1C","#a4afb0","#772E25")) +
  
  ## Set empty theme
  theme_void() +
  theme(legend.title = element_text(colour="black", size=10, face="bold", margin = margin(b = 10)), # margin is space between title and rest of legend
        plot.margin = unit(c(t=1, r=1, b=1, l=1), unit="cm"), # adjust margins
        panel.background = element_rect(fill = 'white', color = 'black'),
        plot.background = element_blank(),
        legend.position = "none")

ggsave(paste("figures/world_map_genetic_pts-phy-Spilhaus-option1.png",sep="") , width = 15, height = 14.5, units = c("cm"))



# option 2 -------------

ggplot() +
  
  # add graticules over oceans
  #tidyterra::geom_spatvector(data=grid.ocean, colour = "gray80", alpha = 0.5, linewidth = 0.25) +
  tidyterra::geom_spatvector(data=grid.ocean, colour = "gray80", alpha = 0.8, linewidth = 0.17, linetype = "dashed") +
  
  # add land
  tidyterra::geom_spatvector(data = land.prj, fill = "gray90", colour = "gray45", alpha = 1) +
  
  ## add genetic points
  geom_point(data = sites.prj, aes(x=spilhaus_x, y=spilhaus_y, fill=taxclade), 
             #position = position_jitter(width = 1000000, height = 1000000),
             shape = 21, colour="black", size = 4, alpha = 0.5) +
  
  scale_fill_manual(values = c("#1F78B4","#A6CEE3","#33A02C","#FF7F00","#FDBF6F","#FB9A99","#B2DF8A","#E31A1C","#a4afb0","#772E25")) +
  
  ## Set empty theme
  theme_void() +
  theme(legend.title = element_text(colour="black", size=10, face="bold", margin = margin(b = 10)), # margin is space between title and rest of legend
        plot.margin = unit(c(t=1, r=1, b=1, l=1), unit="cm"), # adjust margins
        panel.background = element_rect(fill = 'white', color = 'black'),
        plot.background = element_blank(),
        legend.position = "none")

ggsave(paste("figures/world_map_genetic_pts-phy-Spilhaus-option2b.png",sep="") , width = 15, height = 14.5, units = c("cm"))



# option 3 -------------

ggplot() +
  
  # color oceans
  geom_rect(aes(xmin = -11825474-1182547, xmax = 11825474+1182547, ymin = -11825474-1182547, ymax = 11825474+1182547), fill = "#55e2e9", alpha = 0.25) +
  
  # add graticules over oceans
  tidyterra::geom_spatvector(data=grid.ocean, colour = "gray80", alpha = 0.5, linewidth = 0.25) +
  
  # add land
  tidyterra::geom_spatvector(data = land.prj, fill = "white", colour = "gray45", alpha = 1) +
  
  ## add genetic points
  geom_point(data = sites.prj, aes(x=spilhaus_x, y=spilhaus_y, fill=taxclade), 
             #position = position_jitter(width = 1000000, height = 1000000),
             shape = 21, colour="black", size = 4, alpha = 0.5) +
  
  scale_fill_manual(values = c("#1F78B4","#A6CEE3","#33A02C","#FF7F00","#FDBF6F","#FB9A99","#B2DF8A","#E31A1C","#a4afb0","#772E25")) +
  
  ## Set empty theme
  theme_void() +
  theme(legend.title = element_text(colour="black", size=10, face="bold", margin = margin(b = 10)), # margin is space between title and rest of legend
        plot.margin = unit(c(t=1, r=1, b=1, l=1), unit="cm"), # adjust margins
        panel.background = element_rect(fill = 'white', color = 'black'),
        plot.background = element_blank(),
        legend.position = "none")

ggsave(paste("figures/world_map_genetic_pts-phy-Spilhaus-option3.png",sep="") , width = 15, height = 14.5, units = c("cm"))



# option 4 -------------

ggplot() +
  
  # color oceans
  geom_rect(aes(xmin = -11825474-1182547, xmax = 11825474+1182547, ymin = -11825474-1182547, ymax = 11825474+1182547), fill = "#b7ecea", alpha = 0.4) +
  
  # add graticules over oceans
  tidyterra::geom_spatvector(data=grid.ocean, colour = "gray80", alpha = 0.8, linewidth = 0.17, linetype = "dashed") +
  
  # add land
  tidyterra::geom_spatvector(data = land.prj, fill = "white", colour = "gray45", alpha = 1) +
  
  ## add genetic points
  geom_point(data = sites.prj, aes(x=spilhaus_x, y=spilhaus_y, fill=taxclade), 
             #position = position_jitter(width = 1000000, height = 1000000),
             shape = 21, colour="black", size = 4, alpha = 0.5) +
  
  scale_fill_manual(values = c("#1F78B4","#A6CEE3","#33A02C","#FF7F00","#FDBF6F","#FB9A99","#B2DF8A","#E31A1C","#a4afb0","#772E25")) +
  
  ## Set empty theme
  theme_void() +
  theme(legend.title = element_text(colour="black", size=10, face="bold", margin = margin(b = 10)), # margin is space between title and rest of legend
        plot.margin = unit(c(t=1, r=1, b=1, l=1), unit="cm"), # adjust margins
        panel.background = element_rect(fill = 'white', color = 'black'),
        plot.background = element_blank(),
        legend.position = "none")

ggsave(paste("figures/world_map_genetic_pts-phy-Spilhaus-option4.png",sep="") , width = 15, height = 14.5, units = c("cm"))



