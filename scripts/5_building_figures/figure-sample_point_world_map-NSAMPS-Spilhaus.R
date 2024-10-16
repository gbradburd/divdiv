#idea: make a map of the world and plot heat map of N samples on it

#https://spatialreference.org - to look up projection codes, like PROJ.4 codes

#load libraries -------
library(ggplot2)
library(ncdf4)
library(dplyr)
library(tidyr)

rm(list=ls())
gc()



setwd("/Users/rachel/divdiv/")
ecoregionspath = "data/abiotic/input_and_working/Marine_Ecoregions_Of_the_World_(MEOW)-shp/"

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

#get marine ecoregions
#V we actually used in analyses - regions have buffer that extends inland (to avoid weirdness/dropping points right on coasts)
ecor <- terra::vect(paste0(ecoregionspath,"/Marine_Ecoregions_Of_the_World__MEOW_.shp")) %>% 
  terra::project(., PROJ.latlon)
#V that stops at coastline already, that Bruce found online
#aggregate some indiv parts of ecoregions that are split out
#ecor <- terra::aggregate(ecor, by="ECOREGION",fun="modal", count=FALSE)
#ecor<-terra::aggregate(terra::vect("/Users/rachel/data0/meow_ecos_expl_clipped_expl.shp"),
#                       by="ECOREGION",fun="modal",count=FALSE)



# project spatial data into spilhaus ----------------

#preprocess (split regions/lines along seam, etc.) and pull out geom info
#width controls the width of the seam (ensures split polygons end up on correct side)
#might need to be tweaked in some cases
#found that as low as 100 m worked for buffered V of data (that extends beyond coasts)!
#width of 2.5km aka width=2500 needed for V of ecor's that stop at coastline
ecor.coords <- preproc.split(ecor, NA.radius=1066000) %>% terra::geom() %>% as.data.frame()
coast.coords <- preproc.split(NE_coastlines, NA.radius=1066000) %>% terra::geom() %>% as.data.frame()
grid.coords <- preproc.split(grid, NA.radius=1066000) %>% terra::geom() %>% as.data.frame()
land.coords <- preproc.split(land, NA.radius=1066000) %>% terra::geom() %>% as.data.frame()

#build new cols of coords in spilhaus proj
ecor.coords <- from_lonlat_to_spilhaus_xy(ecor.coords$x, ecor.coords$y) %>% cbind(., ecor.coords)
coast.coords <- from_lonlat_to_spilhaus_xy(coast.coords$x, coast.coords$y) %>% cbind(., coast.coords)
grid.coords <- from_lonlat_to_spilhaus_xy(grid.coords$x, grid.coords$y) %>% cbind(., grid.coords)
land.coords <- from_lonlat_to_spilhaus_xy(land.coords$x, land.coords$y) %>% cbind(., land.coords)

# convert layers back to SpatVector (using spilhaus coords)
ecor.prj <- as.matrix(ecor.coords[,c("geom","part","spilhaus_x","spilhaus_y","hole")]) %>% 
  terra::vect(., "polygons", atts=terra::values(ecor))
coast.prj <- as.matrix(coast.coords[,c("geom","part","spilhaus_x","spilhaus_y","hole")]) %>% 
  terra::vect(., "lines", atts=terra::values(NE_coastlines))
grid.prj <- as.matrix(grid.coords[,c("geom","part","spilhaus_x","spilhaus_y","hole")]) %>% 
  terra::vect(., "lines", atts=terra::values(grid))
land.prj <- as.matrix(land.coords[,c("geom","part","spilhaus_x","spilhaus_y","hole")]) %>% 
  terra::vect(., "polygons", atts=terra::values(land))


#get rid of Pacific seams (only really necessary for polygons, not lines)
#three parameters to potentially tweak:
# - intersection.check.tol --> Controls the thickness of line used to check
#     for overlap with "seam" of lonlat projection (180 longitude line).
#     Haven't played around with too much, but 1000 seems fine.
# - snap.tol --> Controls the radius at which borders are "snapped" together
#     (because imprecision of transform makes it so that borders along the
#     seam aren't perfectly touching and/or overlapping).
#     Lower might be ideal, but 5000 was the lowest I could get it to work for
#     all cases!
# - exp.fudge.fac --> Controls how much regions are buffered (i.e., expanded)
#     following snapping--this was necessary in some cases to get things
#     working correctly.
#     Obviously not ideal since it slightly distorts the geometry, but otherwise
#     you need to use some crazy high snap.tol values that ends up making things
#     much worse!
#     100 was the lowest value I could get away with
#all units are in the Spilhaus projection coordinate units, which obviously
# isn't ideal for interpretation...
#the values below are the default values for the function
ecor.prj <- fix.pacific.seams(ecor.prj,
                              intersection.check.tol=1000,
                              snap.tol=5000,
                              exp.fudge.fac=100)



# fix the map corners for polygon layers ----------------

#I don't recommend the above approach for large, complex geometries with many isolated pieces
# becomes unreasonably slow (perhaps add condition to check for number of pieces in function?)
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

ecor.prj <- patch.corners(ecor.prj, 
                          NA.radius=1066000+10,
                          SA.radius=1e6+10,
                          Asia.radius=1e6+10)



# erase some stuff --------------------
#create gridlines just for ocean (aka erase over land)
grid.ocean <- terra::erase(grid.prj, land.prj)

# use landmask to erase inland buffer parts of marine ecoregions -----------
ecor.prj.1 <- terra::erase(ecor.prj, land.prj)
ggplot() + tidyterra::geom_spatvector(data = ecor.prj, fill = "blue", alpha = 0.2) + theme_bw()
ggplot() + tidyterra::geom_spatvector(data = land.prj, fill = "gray80", colour = "black") + 
  tidyterra::geom_spatvector(data = ecor.prj.1, fill = "blue", alpha = 0.2) + theme_bw()



# expand extent of map ---------------
#save og not expanded Vs first
ecor.prj.1.og <- ecor.prj.1
coast.prj.og <- coast.prj
land.prj.og <- land.prj
grid.ocean.og <- grid.ocean

#amount varies from 0 to 1 and controls how far the edges are expanded out
# (in terms of fraction of the base map's dimensions)
#exp.fudge.fac buffers any polygons to ensure shapes aggregate across map borders
ecor.prj.1 <- expand.map(ecor.prj.1, amount=0.05, exp.fudge.fac=10)
coast.prj <- expand.map(coast.prj, amount=0.05, exp.fudge.fac=10)
land.prj <- expand.map(land.prj, amount=0.05, exp.fudge.fac=10)
grid.ocean <- expand.map(grid.ocean, amount=0.05, exp.fudge.fac=10)

ggplot() + 
  tidyterra::geom_spatvector(data=grid.ocean.og, linewidth = 0.1) +
  tidyterra::geom_spatvector(data=ecor.prj.1.og, fill = "blue") +
  tidyterra::geom_spatvector(data=land.prj.og, fill = "forestgreen") +
  tidyterra::geom_spatvector(data=coast.prj.og, colour = "red")

ggplot() + 
  tidyterra::geom_spatvector(data=grid.ocean, linewidth = 0.1) +
  tidyterra::geom_spatvector(data=ecor.prj.1, fill = "blue") +
  tidyterra::geom_spatvector(data=land.prj, fill = "forestgreen") +
  tidyterra::geom_spatvector(data=coast.prj, colour = "red")


# get the genetic sample points and build heatmap of sampling intensity ----------
sites <- read.csv("data/final_genetic_latlongs.csv", header = T)

#project into spilhaus
sites.prj = from_lonlat_to_spilhaus_xy(sites$lon, sites$lat)

#manually build grid and count how many points in each cell
#(can't use raster::rasterize(sites.prj, rast, fun = "count") here)

#set grid size (number of grid cells in Spilhaus x and y coordinates)
n_grid = 50 #30 largest sq size to go probably, 50 smallest
#set limits for grid aka max extent of spilhaus proj
lim <- 11825474
x_lo = -lim
x_hi = lim
y_lo = -lim
y_hi = lim
#build grid
gridded_x_pos = ceiling(n_grid * (sites.prj[,1] - x_lo) / (x_hi - x_lo))
gridded_y_pos = ceiling(n_grid * (sites.prj[,2] - y_lo) / (y_hi - y_lo))
sites_grid = expand.grid(x=seq(x_lo, x_hi, len=n_grid), y=seq(y_lo, y_hi, len=n_grid))
sites_grid$z = 0
#build heat map aka count how many pts in each grid cell
for (i in seq(1, nrow(sites.prj))) {
  pos = gridded_x_pos[i] + (gridded_y_pos[i] - 1) * n_grid
  sites_grid$z[pos] = sites_grid$z[pos] + 1
}
sites_grid$l = (sites_grid$z == 0)
# pretty up
nsamps.res50 = pretify_spilhaus_df(sites_grid)

#build another size of grid so we have options
#set grid size (number of grid cells in Spilhaus x and y coordinates)
n_grid = 36 #30 largest sq size to go probably, 50 smallest
#build grid
gridded_x_pos = ceiling(n_grid * (sites.prj[,1] - x_lo) / (x_hi - x_lo))
gridded_y_pos = ceiling(n_grid * (sites.prj[,2] - y_lo) / (y_hi - y_lo))
sites_grid = expand.grid(x=seq(x_lo, x_hi, len=n_grid), y=seq(y_lo, y_hi, len=n_grid))
sites_grid$z = 0
#build heat map aka count how many pts in each grid cell
for (i in seq(1, nrow(sites.prj))) {
  pos = gridded_x_pos[i] + (gridded_y_pos[i] - 1) * n_grid
  sites_grid$z[pos] = sites_grid$z[pos] + 1
}
sites_grid$l = (sites_grid$z == 0)
# pretty up
nsamps.res36 = pretify_spilhaus_df(sites_grid)



# make maps ! -------------------------

# option 1 -------------
ggplot() +
  
  # add graticules over oceans
  tidyterra::geom_spatvector(data=grid.ocean, colour = "gray70", alpha = 1, linewidth = 0.19, linetype = "dotted") +

  # add land
  tidyterra::geom_spatvector(data = land.prj, fill = "gray80", colour = "gray55", alpha = 0.3) +
  
  # add ecoregions
  tidyterra::geom_spatvector(data = ecor.prj.1, fill = "skyblue", colour = "skyblue", alpha = 0.3, linewidth = 0.25) +
  
  ## add N indivs heatmap/raster
  geom_tile(data = nsamps.res36, aes(x = x, y = y, fill = z)) +
  
  viridis::scale_fill_viridis(name = "N. sequenced\nindividuals", trans = "log",
                              breaks = c(0,1,10,100,450), labels = c(0,1,10,100,450),
                              option="magma") +
  
  ## Set empty theme
  theme_void() +
  theme(legend.title = element_text(colour="black", size=10, face="bold", margin = margin(b = 10)), # margin is space between title and rest of legend
        plot.margin = unit(c(t=1, r=1, b=1, l=1), unit="cm"), # adjust margins
        panel.background = element_rect(fill = 'white', color = 'black'),
        plot.background = element_blank(),
        legend.position = c(0.16, 0.23), legend.direction = "vertical")

ggsave(paste("figures/world_map_genetic_pts-NSAMPS-Spilhaus-option1.png",sep="") , width = 15, height = 14.5, units = c("cm"))



# option 2 -------------
ggplot() +
  
  # add graticules over oceans
  tidyterra::geom_spatvector(data=grid.ocean, colour = "gray80", alpha = 0.5, linewidth = 0.25) +
  
  # add land
  tidyterra::geom_spatvector(data = land.prj, fill = "gray80", colour = "gray45", alpha = 0.3) +
  
  ## add N indivs heatmap/raster
  geom_tile(data = nsamps.res50, aes(x = x, y = y, fill = z)) +
  
  viridis::scale_fill_viridis(name = "N. sequenced\nindividuals", trans = "log",
                              breaks = c(0,1,10,100,450), labels = c(0,1,10,100,450),
                              option="magma") +
  
  ## Set empty theme
  theme_void() +
  theme(legend.title = element_text(colour="black", size=10, face="bold", margin = margin(b = 10)), # margin is space between title and rest of legend
        plot.margin = unit(c(t=1, r=1, b=1, l=1), unit="cm"), # adjust margins
        panel.background = element_rect(fill = 'white', color = 'black'),
        plot.background = element_blank(),
        legend.position = c(0.16, 0.23), legend.direction = "vertical")

ggsave(paste("figures/world_map_genetic_pts-NSAMPS-Spilhaus-option2.png",sep="") , width = 15, height = 14.5, units = c("cm"))



# option 3 -------------
ggplot() +
  
  # add land
  tidyterra::geom_spatvector(data = land.prj, fill = "gray80", colour = "gray55", alpha = 0.3) +
  
  # add ecoregions
  tidyterra::geom_spatvector(data = ecor.prj.1, fill = "skyblue", colour = "skyblue", alpha = 0.3, linewidth = 0.25) +
  
  ## add N indivs heatmap/raster
  geom_tile(data = nsamps.res36, aes(x = x, y = y, fill = z)) +
  
  viridis::scale_fill_viridis(name = "N. sequenced\nindividuals", trans = "log",
                              breaks = c(0,1,10,100,450), labels = c(0,1,10,100,450),
                              option="magma") +
  
  ## Set empty theme
  theme_void() +
  theme(legend.title = element_text(colour="black", size=10, face="bold", margin = margin(b = 10)), # margin is space between title and rest of legend
        plot.margin = unit(c(t=1, r=1, b=1, l=1), unit="cm"), # adjust margins
        panel.background = element_rect(fill = 'white', color = 'black'),
        plot.background = element_blank(),
        legend.position = c(0.16, 0.23), legend.direction = "vertical")

ggsave(paste("figures/world_map_genetic_pts-NSAMPS-Spilhaus-option3.png",sep="") , width = 15, height = 14.5, units = c("cm"))



# option 4 -------------
ggplot() +
  
  # add graticules over oceans
  tidyterra::geom_spatvector(data=grid.ocean, colour = "gray80", alpha = 0.8, linewidth = 0.17, linetype = "dashed") +

  # add land
  tidyterra::geom_spatvector(data = land.prj, fill = "white", colour = "gray55", alpha = 0.3) +
  
  # add ecoregions
  tidyterra::geom_spatvector(data = ecor.prj.1, fill = "skyblue", colour = "skyblue", alpha = 0.3, linewidth = 0.25) +
  
  ## add N indivs heatmap/raster
  geom_tile(data = nsamps.res36, aes(x = x, y = y, fill = z)) +
  
  viridis::scale_fill_viridis(name = "N. sequenced\nindividuals", trans = "log",
                              breaks = c(0,1,10,100,450), labels = c(0,1,10,100,450),
                              option="magma") +
  
  # add coast outlines
  geom_path(data=coast.coords %>% mutate(id = paste0(geom,"_",part)), aes(x=spilhaus_x, y=spilhaus_y, group=id), 
            colour = "black", alpha = 1, linewidth = 0.25) +
  
  ## Set empty theme
  theme_void() +
  theme(legend.title = element_text(colour="black", size=10, face="bold", margin = margin(b = 10)), # margin is space between title and rest of legend
        plot.margin = unit(c(t=1, r=1, b=1, l=1), unit="cm"), # adjust margins
        panel.background = element_rect(fill = 'white', color = 'black'),
        plot.background = element_blank(),
        legend.position = c(0.16, 0.23), legend.direction = "vertical")

ggsave(paste("figures/world_map_genetic_pts-NSAMPS-Spilhaus-option4.png",sep="") , width = 15, height = 14.5, units = c("cm"))



# option 5 -------------
ggplot() +
  
  # add graticules over oceans
  tidyterra::geom_spatvector(data=grid.ocean, colour = "gray80", alpha = 0.5, linewidth = 0.25) +
  
  # add land
  tidyterra::geom_spatvector(data = land.prj, fill = "transparent", colour = "black", alpha = 1, linewidth = 0.25) +
  
  ## add N indivs heatmap/raster
  geom_tile(data = nsamps.res36, aes(x = x, y = y, fill = z)) +
  
  viridis::scale_fill_viridis(name = "N. sequenced\nindividuals", trans = "log",
                              breaks = c(0,1,10,100,450), labels = c(0,1,10,100,450),
                              option="magma") +
  
  # add coast outlines
  geom_path(data=coast.coords %>% mutate(id = paste0(geom,"_",part)), aes(x=spilhaus_x, y=spilhaus_y, group=id), 
            colour = "black", alpha = 1, linewidth = 0.2) +
  
  ## Set empty theme
  theme_void() +
  theme(legend.title = element_text(colour="black", size=10, face="bold", margin = margin(b = 10)), # margin is space between title and rest of legend
        plot.margin = unit(c(t=1, r=1, b=1, l=1), unit="cm"), # adjust margins
        panel.background = element_rect(fill = 'white', color = 'black'),
        plot.background = element_blank(),
        legend.position = c(0.16, 0.23), legend.direction = "vertical")

ggsave(paste("figures/world_map_genetic_pts-NSAMPS-Spilhaus-option5.png",sep="") , width = 15, height = 14.5, units = c("cm"))



# option 6 -------------
ggplot() +
  
  # add ecoregions
  tidyterra::geom_spatvector(data = ecor.prj.1, fill = "skyblue", colour = "skyblue", alpha = 0.3, linewidth = 0.25) +
  
  # add land
  tidyterra::geom_spatvector(data = land.prj, fill = "gray80", colour = "black", alpha = 0.3) +
  
  ## add N indivs heatmap/raster
  geom_tile(data = nsamps.res36, aes(x = x, y = y, fill = z)) +
  
  viridis::scale_fill_viridis(name = "N. sequenced\nindividuals", trans = "log",
                              breaks = c(0,1,10,100,450), labels = c(0,1,10,100,450),
                              option="magma") +
  
  # add coast outlines
  geom_path(data=coast.coords %>% mutate(id = paste0(geom,"_",part)), aes(x=spilhaus_x, y=spilhaus_y, group=id), 
            colour = "black", alpha = 1, linewidth = 0.3) +
  
  ## Set empty theme
  theme_void() +
  theme(legend.title = element_text(colour="black", size=10, face="bold", margin = margin(b = 10)), # margin is space between title and rest of legend
        plot.margin = unit(c(t=1, r=1, b=1, l=1), unit="cm"), # adjust margins
        panel.background = element_rect(fill = 'white', color = 'black'),
        plot.background = element_blank(),
        legend.position = c(0.16, 0.23), legend.direction = "vertical")

ggsave(paste("figures/world_map_genetic_pts-NSAMPS-Spilhaus-option6.png",sep="") , width = 15, height = 14.5, units = c("cm"))
ggsave(paste("figures/world_map_genetic_pts-NSAMPS-Spilhaus-option6.pdf",sep="") , width = 15, height = 14.5, units = c("cm"))



# graveyard -------------

#old way of getting land
#download marine data to locate ocean and
#create mask and use to pull land shape out
landmask <- download_land_mask()
land <- terra::rast(t(ncvar_get(landmask)))
terra::ext(land) <- c(-180,180,-90,90)
land <- land==1
land <- terra::as.polygons(land)[2]



ggplot() +
  
  # add graticules
  
  #full globe
  #geom_path(data=grid.coords %>% mutate(id = paste0(geom,"_",part)), aes(x=spilhaus_x, y=spilhaus_y, group=id), 
  #             colour = "gray80", alpha = 0.5, linewidth = 0.25) +
  
  #just over oceans
  tidyterra::geom_spatvector(data=grid.ocean, colour = "gray80", alpha = 0.5, linewidth = 0.25) +
  
  # add land
  tidyterra::geom_spatvector(data = land.prj, fill = "gray80", colour = "gray55", alpha = 0.3) +
  
  # add great lakes
  tidyterra::geom_spatvector(data = gls.prj, fill = "white", colour = "gray55", linewidth = 0.15) +
  
  # add ecoregions
  tidyterra::geom_spatvector(data = ecor.prj.1, fill = "skyblue", colour = "skyblue", alpha = 0.3, linewidth = 0.25) +
  #tidyterra::geom_spatvector(data = ecor.prj.2, fill = "skyblue", colour = "skyblue", alpha = 0.3, linewidth = 0.3) +
  
  # add coast outlines
  #geom_path(data=coast.coords %>% mutate(id = paste0(geom,"_",part)), aes(x=spilhaus_x, y=spilhaus_y, group=id), 
  #          colour = "gray80", alpha = 0.8, linewidth = 0.15) +
  
  ## add N indivs heatmap/raster
  geom_tile(data = pretty_sites, aes(x = x, y = y, fill = z)) +
  
  viridis::scale_fill_viridis(name = "N. sequenced\nindividuals", trans = "log",
                              breaks = c(0,1,10,100,450), labels = c(0,1,10,100,450),
                              option="magma") +
  
  ## Set empty theme
  theme_void() +
  theme(legend.title = element_text(colour="black", size=10, face="bold", margin = margin(b = 10)), # margin is space between title and rest of legend
        plot.margin = unit(c(t=1, r=1, b=1, l=1), unit="cm"), # adjust margins
        panel.background = element_rect(fill = 'white', color = 'black'),
        plot.background = element_blank(),
        legend.position = c(0.13, 0.2), legend.direction = "vertical")


