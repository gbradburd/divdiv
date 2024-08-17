#idea: make a map of the world and plot heat map of N samples on it

#https://spatialreference.org - to look up projection codes, like PROJ.4 codes

#load libraries -------

library(rerddap)
library(ggplot2)
library(ncdf4)
library(dplyr)
library(tidyr)

rm(list=ls())
gc()


setwd("/Users/rachel/divdiv/")
sampkeypath = "/Users/rachel/divdiv/data/all_samplenamekeys"
wmpath = "/Users/rachel/divdiv/data/popgen/input_and_working/ALL_r80_gendiv_data/"


#load Spilhaus projection functions
source("scripts/spilhaus_functions.R") 



# make spilhaus map -----------------
#create a data frame with NxN pixels
#higher spilhaus_res = more defined coastline 
spilhaus_df = make_spilhaus_xy_gridpoints(spilhaus_res=1000)
#convert spilhaus coors to mercator
lonlat = from_spilhaus_xy_to_lonlat(spilhaus_df$x, spilhaus_df$y)
# download netCDF data from Coastwatch so we can denote where ocean vs land is (will only get data for ocean pts)
da = download_sst_data("2021-09-16", "2021-09-16")
# extract data for our lats and lons
spilhaus_df$z = extract_sst_data(da, lonlat)
# mask
spilhaus_df$l = is.na(spilhaus_df$z)
# prettify
pretty_spilhaus_df = pretify_spilhaus_df(spilhaus_df)



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



#code expects lat and long to be first two cols and to be called lat and lon
sites <- sites %>% dplyr::select(lon,lat,everything())
sites.df <- sites

#project sites
PROJ.pts <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" #EPSG:4326, WGS 84
#first make dataframe into spatial object, using projection points are in
sites.prj <- sf::st_as_sf(sites, coords = c("lon", "lat"), crs = PROJ.pts)



# trying something -------------

spil.latlon <- as.data.frame(from_spilhaus_xy_to_lonlat(pretty_spilhaus_df$x, pretty_spilhaus_df$y)) %>% 
  filter(is.na(longitude)==F)
spil.latlon.prj <- sf::st_as_sf(spil.latlon, coords = c("longitude", "latitude"), crs = PROJ.pts)

rast <- raster::raster(res = 5)
raster::extent(rast) <- raster::extent(spil.latlon.prj)
rast
raster.Nindivs <- raster::rasterize(sites.prj, rast, fun = "count")
raster.Nindivs <- raster::rasterToPoints(raster.Nindivs) %>% as.data.frame()

raster.Nindivs.prj <- as.data.frame(from_lonlat_to_spilhaus_xy(raster.Nindivs$x, raster.Nindivs$y))
raster.Nindivs.prj <- cbind(raster.Nindivs.prj, raster.Nindivs)

ggplot() + geom_sf(data = spil.latlon.prj)
ggplot() + geom_point(data = spil.latlon, aes(x=longitude,y=latitude))
ggplot() + geom_tile(data = raster.Nindivs, aes(x = x, y = y, fill = ID)) + coord_fixed()
ggplot() + geom_tile(data = raster.Nindivs.prj, aes(x = spilhaus_x, y = spilhaus_y, fill = ID)) + coord_fixed()


#make heatmap raster - N number of points in each box (large res = fewer/bigger squares)
rast <- raster::raster(res = 5)
raster::extent(rast) <- raster::extent(sites.prj)
rast
raster.Nindivs <- raster::rasterize(sites.prj, rast, fun = "count")
raster.Nindivs <- raster::rasterToPoints(raster.Nindivs) %>% as.data.frame()

#project genetic lat longs to spilhaus
raster.Nindivs.prj <- as.data.frame(from_lonlat_to_spilhaus_xy(raster.Nindivs$x, raster.Nindivs$y))
test <- cbind(raster.Nindivs, raster.Nindivs.prj)

ggplot() +
  geom_tile(data=pretty_spilhaus_df, aes(x=x, y=y), fill = "white") +
  geom_point(data = test, aes(x = spilhaus_x, y = spilhaus_y,  colour = ID), shape = 15, size = 3) +
  viridis::scale_colour_viridis(name = "N. samples", trans = "log",
                              breaks = c(0,1,10,100,500), labels = c(0,1,10,100,500),
                              option="mako") +
  theme(panel.background = element_rect(fill = 'gray80', color = 'gray80'),
        panel.grid = element_blank(),
        #legend.position = "none",
        plot.background = element_rect(fill = "transparent")) +
  coord_equal()

ggsave(filename = "figures/world_map_genetic_pts-NSAMPS-spilhaus.pdf", width = 30, height = 15, units = c("cm"))



