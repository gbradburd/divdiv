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
sampkeypath = "data/all_samplenamekeys"
wmpath = "data/popgen/input_and_working/ALL_r80_gendiv_data/"
ecoregionspath = "data/abiotic/input_and_working/Marine_Ecoregions_Of_the_World_(MEOW)-shp/"


#load Spilhaus projection functions
source("scripts/spilhaus_functions.R") 



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

#set corners so we can stretch/extend layers to edges of map later
lim <- 11825474
wd1 <- 0.575*lim
wd2 <- 0.15*lim
corners <- list(cbind(c(-lim,-lim+wd1,-lim+wd1,-lim),
                      c(-lim,-lim,-lim+wd1,-lim+wd1)),
                cbind(c(lim-wd1,lim,lim,lim-wd1),
                      c(lim-wd1,lim-wd1,lim,lim)),
                cbind(c(-lim,-lim+wd2,-lim+wd2,-lim),
                      c(lim-wd2,lim-wd2,lim,lim)),
                cbind(c(lim-wd2,lim,lim,lim-wd2),
                      c(-lim,-lim,-lim+wd2,-lim+wd2)))



# project spatial data into spilhaus ----------------

#preprocess (split regions/lines along seam, etc.) and pull out geom info
#width controls the width of the seam (ensures split polygons end up on correct side)
#might need to be tweaked in some cases
#found that as low as 100 m worked for buffered V of data (that extends beyond coasts)!
#width of 2.5km aka width=2500 needed for V of ecor's that stop at coastline
ecor.coords <- preproc.split(ecor, width=100) %>% terra::geom() %>% as.data.frame()
coast.coords <- preproc.split(NE_coastlines, width=100) %>% terra::geom() %>% as.data.frame()
grid.coords <- preproc.split(grid, width=100) %>% terra::geom() %>% as.data.frame()
land.coords <- preproc.split(NE_countries, width=100) %>% terra::geom() %>% as.data.frame()
gls.coords <- preproc.split(gls, width=100) %>% terra::geom() %>% as.data.frame()

#build new cols of coords in spilhaus proj
ecor.coords <- from_lonlat_to_spilhaus_xy(ecor.coords$x, ecor.coords$y) %>% cbind(., ecor.coords)
coast.coords <- from_lonlat_to_spilhaus_xy(coast.coords$x, coast.coords$y) %>% cbind(., coast.coords)
grid.coords <- from_lonlat_to_spilhaus_xy(grid.coords$x, grid.coords$y) %>% cbind(., grid.coords)
land.coords <- from_lonlat_to_spilhaus_xy(land.coords$x, land.coords$y) %>% cbind(., land.coords)
gls.coords <- from_lonlat_to_spilhaus_xy(gls.coords$x, gls.coords$y) %>% cbind(., gls.coords)



# convert great lakes back to SpatVector (using spilhaus coords)
gls.prj <- as.matrix(gls.coords[,c("geom","part","spilhaus_x","spilhaus_y","hole")]) %>% 
  terra::vect(., "polygons", atts=terra::values(gls))



# tidy up land -----------------------
#polygons won't plot correctly using coords and geom_poly (I think issue with holes)
#so convert back to SpatVector (using spilhaus coords)
#ggplot() + geom_polygon(data=land.coords %>% mutate(id = paste0(geom,"_",part)), aes(x=x, y=y, group=id), fill = "gray60", alpha = 1)
#ggplot() + geom_polygon(data=land.coords %>% mutate(id = paste0(geom,"_",part)), aes(x=spilhaus_x, y=spilhaus_y, group=id), fill = "gray60", alpha = 1)
land.prj <- as.matrix(land.coords[,c("geom","part","spilhaus_x","spilhaus_y","hole")]) %>% 
  terra::vect(., "polygons", atts=terra::values(NE_countries))
ggplot() + tidyterra::geom_spatvector(data = land.prj) + theme_bw()
#and fix corners of land
land.prj <- terra::buffer(land.prj,0) %>% terra::union(., terra::vect(corners,"polygons")) %>% terra::aggregate()
ggplot() + tidyterra::geom_spatvector(data = land.prj) + theme_bw()



# tidy up ecoregions -----------
#convert back to SpatVector (using spilhaus coords)
ecor.prj <- as.matrix(ecor.coords[,c("geom","part","spilhaus_x","spilhaus_y","hole")]) %>% 
  terra::vect(., "polygons", atts=terra::values(ecor))

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

#option1: use landmask to erase inland buffer parts of marine ecoregions
ecor.prj.1 <- terra::erase(ecor.prj, land.prj)
ggplot() + tidyterra::geom_spatvector(data = ecor.prj, fill = "blue", alpha = 0.2) + theme_bw()
ggplot() + tidyterra::geom_spatvector(data = land.prj, fill = "gray80", colour = "black") + 
  tidyterra::geom_spatvector(data = ecor.prj.1, fill = "blue", alpha = 0.2) + theme_bw()

#option2: use above corner-padding trick to fix corner of offending ecoregion
corners<-cbind(c(lim-wd1,lim,lim,lim-wd2),
               c(lim,lim,lim-wd1,lim-wd2))
polygon(corners)
prob <- which(terra::values(ecor.prj)$ECOREGION=="East China Sea")
tmp <- terra::aggregate(terra::union(ecor.prj[prob],
                                     terra::vect(corners,"polygon")))
terra::values(tmp) <- terra::values(ecor.prj)[prob,]
ecor.prj.2 <- rbind(ecor.prj,tmp)
ecor.prj.2 <- ecor.prj.2[-prob]
ggplot() + tidyterra::geom_spatvector(data = land.prj, fill = "gray80", colour = "black") + 
  tidyterra::geom_spatvector(data = ecor.prj.2, fill = "blue", alpha = 0.2) + theme_bw()



# get the genetic sample points ----------

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



# make raster heatmap to show sampling intensity and project into sphilhaus -----------
sites.prj = from_lonlat_to_spilhaus_xy(sites$lon, sites$lat)

#manually build grid and count how many points in each cell
#(can't use raster::rasterize(sites.prj, rast, fun = "count") here)

#set grid size (number of grid cells in Spilhaus x and y coordinates)
n_grid = 50 #30 largest sq size to go probably, 50 smallest
#build grid
x_lo = -lim
x_hi = lim
y_lo = -lim
y_hi = lim
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

#set grid size (number of grid cells in Spilhaus x and y coordinates)
n_grid = 36 #30 largest sq size to go probably, 50 smallest
#build grid
x_lo = -lim
x_hi = lim
y_lo = -lim
y_hi = lim
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


# create gridlines just for ocean (aka erase over land) ------------
grid.prj <- as.matrix(grid.coords[,c("geom","part","spilhaus_x","spilhaus_y","hole")]) %>% 
  terra::vect(., "lines", atts=terra::values(grid))
#erase Great Lakes from land
land.sans.gls.prj <- terra::erase(land.prj, gls.prj)
#erase land sans GLs from gridlines
grid.ocean <- terra::erase(grid.prj, land.sans.gls.prj)



# make maps ! -------------------------

# option 1 -------------
ggplot() +
  
  # add graticules over oceans
  tidyterra::geom_spatvector(data=grid.ocean, colour = "gray70", alpha = 1, linewidth = 0.19, linetype = "dotted") +

  # add land
  tidyterra::geom_spatvector(data = land.sans.gls.prj, fill = "gray80", colour = "gray55", alpha = 0.3) +
  
  # add great lakes
  tidyterra::geom_spatvector(data = gls.prj, fill = "transparent", colour = "gray55", linewidth = 0.15) +
  
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
        legend.position = c(0.13, 0.2), legend.direction = "vertical")

ggsave(paste("figures/world_map_genetic_pts-NSAMPS-Spilhaus-option1.png",sep="") , width = 15, height = 14.5, units = c("cm"))



# option 2 -------------
ggplot() +
  
  # add graticules over oceans
  tidyterra::geom_spatvector(data=grid.ocean, colour = "gray80", alpha = 0.5, linewidth = 0.25) +
  
  # add land
  tidyterra::geom_spatvector(data = land.sans.gls.prj, fill = "gray80", colour = "gray45", alpha = 0.3) +
  
  # add great lakes
  tidyterra::geom_spatvector(data = gls.prj, fill = "transparent", colour = "gray45", linewidth = 0.15) +
  
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
        legend.position = c(0.13, 0.2), legend.direction = "vertical")

ggsave(paste("figures/world_map_genetic_pts-NSAMPS-Spilhaus-option2.png",sep="") , width = 15, height = 14.5, units = c("cm"))



# option 3 -------------
ggplot() +
  
  # add land
  tidyterra::geom_spatvector(data = land.sans.gls.prj, fill = "gray80", colour = "gray55", alpha = 0.3) +
  
  # add great lakes
  tidyterra::geom_spatvector(data = gls.prj, fill = "transparent", colour = "gray55", linewidth = 0.15) +
  
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
        legend.position = c(0.13, 0.2), legend.direction = "vertical")

ggsave(paste("figures/world_map_genetic_pts-NSAMPS-Spilhaus-option3.png",sep="") , width = 15, height = 14.5, units = c("cm"))



# option 4 -------------
#dashed, linewidth=0.17, alpha=0.8
#dotted, linewidth=0.19, alpha=1
ggplot() +
  
  # add graticules over oceans
  tidyterra::geom_spatvector(data=grid.ocean, colour = "gray80", alpha = 0.8, linewidth = 0.17, linetype = "dashed") +

  # add great lakes
  tidyterra::geom_spatvector(data = gls.prj, fill = "transparent", colour = "gray55", linewidth = 0.15) +
  
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
        legend.position = c(0.13, 0.2), legend.direction = "vertical")

ggsave(paste("figures/world_map_genetic_pts-NSAMPS-Spilhaus-option4.png",sep="") , width = 15, height = 14.5, units = c("cm"))



# option 5 -------------
ggplot() +
  
  # add graticules over oceans
  tidyterra::geom_spatvector(data=grid.ocean, colour = "gray80", alpha = 0.5, linewidth = 0.25) +
  
  # add great lakes
  tidyterra::geom_spatvector(data = gls.prj, fill = "transparent", colour = "gray55", linewidth = 0.15) +
  
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
        legend.position = c(0.13, 0.2), legend.direction = "vertical")

ggsave(paste("figures/world_map_genetic_pts-NSAMPS-Spilhaus-option5.png",sep="") , width = 15, height = 14.5, units = c("cm"))



# option 6 -------------
ggplot() +
  
  # add land
  tidyterra::geom_spatvector(data = land.sans.gls.prj, fill = "gray80", colour = "black", alpha = 0.3) +
  
  # add great lakes
  tidyterra::geom_spatvector(data = gls.prj, fill = "transparent", colour = "black", linewidth = 0.25) +
  
  # add ecoregions
  tidyterra::geom_spatvector(data = ecor.prj.1, fill = "skyblue", colour = "skyblue", alpha = 0.3, linewidth = 0.25) +
  
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
        legend.position = c(0.13, 0.2), legend.direction = "vertical")

ggsave(paste("figures/world_map_genetic_pts-NSAMPS-Spilhaus-option6.png",sep="") , width = 15, height = 14.5, units = c("cm"))



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


