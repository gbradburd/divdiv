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
sampkeypath = "/Users/rachel/divdiv/data/all_samplenamekeys"
wmpath = "/Users/rachel/divdiv/data/popgen/input_and_working/ALL_r80_gendiv_data/"


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

ggsave(paste("figures/world_map_genetic_pts-NSAMPS-Spilhaus-36grid.pdf",sep="") , width = 15, height = 14.5, units = c("cm"))
ggsave(paste("figures/world_map_genetic_pts-NSAMPS-Spilhaus-36grid.png",sep="") , width = 15, height = 14.5, units = c("cm"))

max(pretty_sites$z)


