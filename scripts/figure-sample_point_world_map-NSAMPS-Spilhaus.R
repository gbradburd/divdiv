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



# make spilhaus map -----------------
#create Splihaus map grid (table of x, y points)
spilhaus_df = make_spilhaus_xy_gridpoints(spilhaus_res=1500)
x_lo = min(spilhaus_df$x); x_hi = max(spilhaus_df$x)
y_lo = min(spilhaus_df$y); y_hi = max(spilhaus_df$y)

#convert grid to lat/lon
lonlat = from_spilhaus_xy_to_lonlat(spilhaus_df$x, spilhaus_df$y)

#download marine data to locate ocean and 
#create mask aka pull ocean shape out of square grid
da = download_land_mask()
spilhaus_df$z = extract_mask(da, lonlat)
spilhaus_df$l = is.na(spilhaus_df$z)



# get world coastlines and convert to spilhaus -------------
NE_coastlines <- rnaturalearth::ne_coastline(scale = "medium", returnclass = 'sf')
#convert to spatvect so can use terra funcs
PROJ.latlon <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
NE_coastlines.prj <- sf::st_transform(NE_coastlines, crs = PROJ.latlon)
NE_coastlines.prj <- as(NE_coastlines.prj, "Spatial")
NE_coastlines.prj <- terra::vect(NE_coastlines.prj)
#pull out lat/longs
coasts <- terra::geom(NE_coastlines.prj) %>% as.data.frame() %>% mutate(grp = paste0(geom,"_",part))
#do they look right in current proj?
ggplot() + geom_polygon(data = coasts, aes(x=x, y=y, group=grp), fill = "transparent", colour = "orange", alpha = 0.2)
#convert lat/longs to spilhaus and add grouping variable back on
coasts.prj = from_lonlat_to_spilhaus_xy(coasts$x, coasts$y)
coasts.prj <- cbind(coasts.prj,coasts)
head(coasts.prj)
#deal with polygons that get torn apart due to discontinuities near Panama and Bering Straights
coasts.prj$grp2 <- c(0, cumsum(mapply(2:nrow(coasts.prj), FUN = function(i) {
  sqrt((coasts.prj$spilhaus_x[i] - coasts.prj$spilhaus_x[i-1])^2 + 
         (coasts.prj$spilhaus_y[i] - coasts.prj$spilhaus_y[i-1])^2) > 2e5
})))



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



# do site projecting and heat/raster layer making -----------
sites.prj = from_lonlat_to_spilhaus_xy(sites$lon, sites$lat)

#manually build grid and count how many points in each cell
#(can't use raster::rasterize(sites.prj, rast, fun = "count") here)
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



# pretty up projected layers ----------------
pretty_spilhaus_df = pretify_spilhaus_df(spilhaus_df)
pretty_sites = pretify_spilhaus_df(sites_grid)



# make heatmap - cells colored by N individuals - white ocean --------------
ggplot() +
  
  # add ocean
  geom_tile(data=pretty_spilhaus_df, aes(x=x, y=y), fill = "white") +
  
  # add land outlines
  #geom_sf(data = NE_coastlines.prj, 
  #        colour = "gray50", linewidth = .25) +
  geom_path(data = coasts.prj, aes(x=spilhaus_x, y=spilhaus_y, group=grp2), colour = "gray70", alpha = 1, linewidth = 0.25) +
  
  ## add N indivs raster
  geom_tile(data = pretty_sites, aes(x = x, y = y, fill = z)) +
  
  viridis::scale_fill_viridis(name = "N. sequenced\nindividuals", trans = "log",
                              breaks = c(0,1,10,100,450), labels = c(0,1,10,100,450),
                              option="mako") +
  
  ## Set empty theme
  theme_void() + # remove the default background, gridlines & default gray color around legend's symbols
  theme(legend.title = element_text(colour="black", size=10, face="bold", margin = margin(b = 10)), # margin is space between title and rest of legend
        plot.margin = unit(c(t=1, r=1, b=1, l=1), unit="cm"), # adjust margins
        panel.background = element_rect(fill = 'gray85', color = 'black'),
        plot.background = element_blank(),
        legend.position = c(0.13, 0.2), legend.direction = "vertical") +
  coord_equal()

ggsave(paste("figures/world_map_genetic_pts-NSAMPS-Spilhaus-36grid-outline.pdf",sep="") , width = 15, height = 14.5, units = c("cm"))
ggsave(paste("figures/world_map_genetic_pts-NSAMPS-Spilhaus-36grid-outline.png",sep="") , width = 15, height = 14.5, units = c("cm"))



# get marine ecoregions ------------------
ecor <- terra::vect("/Users/rachel/divdiv/data/abiotic/input_and_working/Marine_Ecoregions_Of_the_World_(MEOW)-shp/Marine_Ecoregions_Of_the_World__MEOW_.shp")
#reproject into norm lat/long and pull out
PROJ.latlon <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
ecor <- terra::project(ecor, PROJ.latlon)
#pull out lat/longs
ecor <- terra::geom(ecor) %>% as.data.frame() %>% mutate(grp = paste0(geom,"_",part))
#do they look right in current proj?
ggplot() + geom_polygon(data = ecor, aes(x=x, y=y, group=grp), fill = "transparent", colour = "orange", alpha = 0.2)
#convert lat/longs to spilhaus and add grouping variable back on
ecor.prj = from_lonlat_to_spilhaus_xy(ecor$x, ecor$y)
ecor.prj <- cbind(ecor.prj,ecor)
head(ecor.prj)
#deal with polygons that get torn apart due to discontinuities near Panama and Bering Straights
ecor.prj$grp2 <- c(0, cumsum(mapply(2:nrow(ecor.prj), FUN = function(i) {
  sqrt((ecor.prj$spilhaus_x[i] - ecor.prj$spilhaus_x[i-1])^2 + 
         (ecor.prj$spilhaus_y[i] - ecor.prj$spilhaus_y[i-1])^2) > 2e5
})))

#hmmm, not great looking
ggplot() + geom_path(data = ecor.prj, aes(x=spilhaus_x, y=spilhaus_y, group=grp2), colour = "orange", alpha = 1)
ggplot() + geom_path(data = ecor.prj, aes(x=spilhaus_x, y=spilhaus_y, group=grp), colour = "orange", alpha = 1)

#find the polygons that are problems using some math
ecor.prj %>% group_by(grp) %>% mutate(range.x = max(spilhaus_x)-min(spilhaus_x)) %>% 
  mutate(range.y = max(spilhaus_y)-min(spilhaus_y)) %>% 
  dplyr::select(geom,grp,range.x,range.y) %>% 
  distinct() %>%
  pivot_longer(., names_to = "dim", values_to = "range", cols = range.x:range.y) %>% 
  ggplot() + geom_histogram(aes(x=range)) + facet_wrap(.~dim, ncol = 1)

#pull out the prob polygons
probs <- ecor.prj %>% group_by(grp) %>% mutate(range.x = max(spilhaus_x)-min(spilhaus_x)) %>% 
  mutate(range.y = max(spilhaus_y)-min(spilhaus_y)) %>% 
  dplyr::select(geom,grp,range.x,range.y) %>% 
  distinct() %>%
  pivot_longer(., names_to = "dim", values_to = "range", cols = range.x:range.y) %>% 
  filter(range > 5e6)

#viz where prob polygons are
probs.plot <- ecor.prj %>% filter(grp %in% probs$grp)

ggplot() + 
  geom_path(data = ecor.prj, aes(x=spilhaus_x, y=spilhaus_y, group=grp), colour = "orange", alpha = 1) +
  geom_path(data = probs.plot, aes(x=spilhaus_x, y=spilhaus_y, group=grp), colour = "red", alpha = 1)

ggplot() + 
  geom_polygon(data = ecor, aes(x=x, y=y, group=grp), fill = "transparent", colour = "orange", alpha = 0.2) +
  geom_polygon(data = ecor %>% filter(grp %in% probs$grp), aes(x=x, y=y, group=grp), fill = "red", colour = "red", alpha = 0.2)



#look at just one of probs
head(probs)
prob = "34_1"
prob = "37_1"
prob = "42_1"

#what it looks like
ggplot() + 
  geom_path(data = ecor.prj %>% filter(!grp %in% probs$grp), aes(x=spilhaus_x, y=spilhaus_y, group=grp), colour = "orange", alpha = 1) +
  geom_path(data = ecor.prj %>% filter(grp == prob), aes(x=spilhaus_x, y=spilhaus_y, group=grp), colour = "red", alpha = 1)
#what it should look like
ggplot() + 
  geom_polygon(data = ecor, aes(x=x, y=y, group=grp), fill = "transparent", colour = "orange", alpha = 0.2) +
  geom_polygon(data = ecor %>% filter(grp == prob), aes(x=x, y=y, group=grp), fill = "red", colour = "red", alpha = 0.2)
#single out the specific prob points in polygon
ggplot() + 
  geom_path(data = ecor.prj %>% filter(!grp %in% probs$grp), aes(x=spilhaus_x, y=spilhaus_y, group=grp), colour = "orange", alpha = 1) +
  geom_path(data = ecor.prj %>% filter(grp == prob), aes(x=spilhaus_x, y=spilhaus_y, group=grp), colour = "red", alpha = 1) +
  geom_point(data = ecor.prj %>% filter(grp == prob) %>% filter(spilhaus_x < -1e7), aes(x=spilhaus_x, y=spilhaus_y, group=grp), colour = "purple")
#try to prettify
test <- ecor.prj %>% filter(grp == prob) %>% filter(spilhaus_x < -1e7) %>%
  rename("lon"="x", "lat" = "y", "x" = "spilhaus_x", "y" = "spilhaus_y") %>% 
  mutate(l = FALSE) %>%
  mutate(z = grp) %>% 
  dplyr::select(x, y, z, l)
test.pretty = pretify_spilhaus_df(test)
#points do move, but not to the right place
ggplot() + 
  geom_path(data = ecor.prj %>% filter(!grp %in% probs$grp), aes(x=spilhaus_x, y=spilhaus_y, group=grp), colour = "orange", alpha = 1) +
  geom_path(data = ecor.prj %>% filter(grp == prob), aes(x=spilhaus_x, y=spilhaus_y, group=grp), colour = "red", alpha = 1) +
  geom_point(data = ecor.prj %>% filter(grp == prob) %>% filter(spilhaus_x < -1e7), aes(x=spilhaus_x, y=spilhaus_y, group=grp), colour = "purple") +
  geom_point(data = test.pretty, aes(x=x, y=y), colour = "blue")
#





test <- ecor.prj %>% filter(grp == "34_1") %>% filter(spilhaus_x < -1e7)

test <- test %>% mutate(new = )



for(i in 1:length(unique(ecor.prj$grp))) {
  
  list <- unique(ecor.prj$grp)
  ecor.i <- list[i]
  print(paste0("ecor is ", ecor.i))
  tmp <- ecor.prj %>% filter(grp == ecor.i)
  out <- ggplot() +
    geom_path(data = tmp, aes(x=spilhaus_x, y=spilhaus_y, group=grp), colour = "orange", alpha = 1) +
    labs(title = ecor.i)
  print(out)

}



