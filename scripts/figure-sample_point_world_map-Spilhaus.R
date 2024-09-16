# EXAMPLE: 
# Extract sea surface temperature data from a remote repo and plot it using the Spilhaus projection

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
#download netCDF data from Coastwatch so we can denote where ocean vs land is (will only get data for ocean pts)
da = download_sst_data("2021-09-16", "2021-09-16")
#extract data for our lats and lons
spilhaus_df$z = extract_sst_data(da, lonlat)
#mask
spilhaus_df$l = is.na(spilhaus_df$z)
#prettify
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



#add clade groups on for coloring
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

sites <- sites %>% mutate(pt_ID = paste0("pt_",1:nrow(.)))

#project genetic lat longs to spilhaus
sites.prj <- as.data.frame(from_lonlat_to_spilhaus_xy(sites$lon, sites$lat))
sites.prj <- cbind(sites.prj, sites)



# plot ----------------------
ggplot() +
  geom_tile(data=pretty_spilhaus_df, aes(x=x, y=y), fill = "white") +
  geom_point(data = sites.prj, aes(x=spilhaus_x, y=spilhaus_y, fill=simplerclades), 
             #position = position_jitter(width = 1000000, height = 1000000),
             shape = 21, colour="black", size = 4, alpha = 0.5) +
  scale_fill_manual(values = c("#1F78B4","#A6CEE3","#33A02C","#FF7F00","#FDBF6F","#FB9A99","#B2DF8A","#E31A1C","transparent","#CAB2D6")) +
  theme(panel.background = element_rect(fill = 'gray80', color = 'gray80'),
        panel.grid = element_blank(),
        legend.position = "none",
        plot.background = element_rect(fill = "transparent"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) +
  coord_equal()

ggsave(filename = "figures/world_map_genetic_pts-spilhaus.pdf", width = 30, height = 15, units = c("cm"))




# problems -------------

probs <- sites.prj %>% filter(spilhaus_x > 0) %>% filter(spilhaus_y > 6000000)
sites.probs <- sites %>% filter(pt_ID %in% probs$pt_ID)


