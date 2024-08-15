# EXAMPLE: 
# Extract sea surface temperature data from a remote repo and plot it using the Spilhaus projection

library(rerddap)
library(ggplot2)
library(ncdf4)
library(dplyr)
library(tidyr)

rm(list=ls())
gc()

setwd("/Users/rachel/divdiv")



#load Spilhaus projection functions
source("scripts/spilhaus_functions.R") 

# auxiliary function to extract netCDF file from Coastwatch
download_sst_data = function(start_date, end_date) {
  myInfo = rerddap::info('NOAA_DHW_monthly', url='https://coastwatch.pfeg.noaa.gov/erddap/')
  myData = rerddap::griddap(datasetx = myInfo,
                            fields = "sea_surface_temperature",
                            latitude = c(-89.975, 89.975),
                            longitude = c(-179.975, 179.975),
                            time = c(start_date, end_date))
  da = ncdf4::nc_open(myData$summary$filename)
  return(da)
}

# auxiliary function to extract a list (length = number of instants) 
#of snapshots (7200x3600 pixels) from the netCDF file
extract_sst_data = function(da, lonlat) {
  lons = sort(ncdf4::ncvar_get(da, "longitude"))
  lats = sort(ncdf4::ncvar_get(da, "latitude"))
  times = ncdf4::ncvar_get(da, "time")

  dlon = lons[2] - lons[1]
  dlat = lats[2] - lats[1]
  ln = as.integer(round((lonlat[,1] - lons[1]) / dlon) + 1)
  la = as.integer(round((lonlat[,2] - lats[1]) / dlat) + 1)
  get_chunk = function(timepoint) {
    sst = ncdf4::ncvar_get(da, "sea_surface_temperature",
                           start = c(1,1,timepoint), count = c(7200, 3600, 1))
    sst = sst[,ncol(sst):1]
    chunk = sst[ln + (la - 1) * dim(sst)[1]]
    return(chunk)
  }
  sst_data = mapply(1:length(times), FUN=function(timepoint) {get_chunk(timepoint)})
  return(sst_data)
}
                    
#first crate a data frame with NxN pixels
#higher spilhaus_res = more defined coastline 
spilhaus_df = make_spilhaus_xy_gridpoints(spilhaus_res=1000)
# convert the spilhaus coordinates into Mercator coordinates
lonlat = from_spilhaus_xy_to_lonlat(spilhaus_df$x, spilhaus_df$y)

# download the netCDF data from Coastwatch
da = download_sst_data("2021-09-16", "2021-09-16")

# extract the SST data for the required lats and lons
spilhaus_df$z = extract_sst_data(da, lonlat)
# mask
spilhaus_df$l = is.na(spilhaus_df$z)
# prettify
pretty_spilhaus_df = pretify_spilhaus_df(spilhaus_df)



# get the genetic sample points we want to plot on map ----------

#get master list of datasets we're using
using <- read.csv("/Users/rachel/divdiv/data/master_df.csv") %>% 
  mutate(link = gsub("bioprj_","",run_name)) %>% dplyr::select(species,link)

#get lat/longs
datasets <- list.files(path = "/Users/rachel/divdiv/data/abiotic/input_and_working/lat_long_tables_per_dataset",
                       pattern = "lat_long_table", full.names = TRUE)
sites <- data.frame("link"=NA, "run_acc_sra"=NA, "lat"=NA, "long"=NA, "coordinateUncertaintyInMeters"=NA, "origin"=NA)
for ( d in datasets ) {
  out <- read.delim(d)
  sites <- rbind(sites,out)
}
sites <- sites %>% filter(is.na(link)==F) %>% dplyr::rename("lon" = "long")

#keep ones we want and pull out species name
sites <- sites %>% filter(link %in% using$link) %>% 
  separate(., link, into = c("garbage","species"), sep = "_", remove = F) %>% 
  dplyr::select(-garbage) %>% mutate(species = gsub("-","_",species))

#add clade groups on for coloring
taxcolorkey <- read.csv("/Users/rachel/divdiv/data/master_tax_color_key.csv") %>% mutate(species = gsub(" ","_",species))
sites <- merge(sites, 
               taxcolorkey, 
               by = "species", all.x = T)

sites <- sites %>% mutate(pt_ID = paste0("pt_",1:nrow(.)))

#project genetic lat longs to spilhaus
sites.prj <- as.data.frame(from_lonlat_to_spilhaus_xy(sites$lon, sites$lat))
sites.prj <- cbind(sites.prj, sites)

# plot ----------------------
ggplot() +
  geom_raster(data=pretty_spilhaus_df, aes(x=x, y=y), fill = "white") +
  geom_point(data = sites.prj, aes(x=spilhaus_x, y=spilhaus_y, fill=taxclade), shape = 21, colour="black", size = 4, alpha = 0.5) +
  scale_fill_manual(values = c("#A6CEE3","#33A02C","#FF7F00","#FDBF6F","#FB9A99","#B2DF8A","#E31A1C","#C2C0C0","#6A3D9A","#1F78B4","#CAB2D6")) +
  theme(panel.background = element_rect(fill = 'gray80', color = 'gray80'),
        panel.grid = element_blank(),
        legend.position = "none",
        plot.background = element_rect(fill = "transparent")) +
  coord_equal()

ggsave(filename = "figures/world_map_genetic_pts-spilhaus.pdf", width = 30, height = 15, units = c("cm"))




# problems -------------

probs <- sites.prj %>% filter(spilhaus_x > 0) %>% filter(spilhaus_y > 6000000)
sites.probs <- sites %>% filter(pt_ID %in% probs$pt_ID)

load("data/abiotic/input_and_working/earth_map_objects.RData")

# get Great Lakes outlines
lakes <- rnaturalearth::ne_download(scale="medium", category = 'physical', type = "lakes", returnclass = "sf")
#subset world lake outlines to just GLs
gls <- c("Lake Superior", "Lake Huron",
         "Lake Michigan", "Lake Erie", "Lake Ontario")
gls <- subset(lakes, name %in% gls)

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
lbl.X.prj <- sf::st_as_sf(lbl.X, coords = c("lon", "lat"), crs = PROJ)
lbl.Y.prj <- sf::st_as_sf(lbl.Y, coords = c("lon", "lat"), crs = PROJ)

#project sites
PROJ.pts <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" #EPSG:4326, WGS 84
#first make dataframe into spatial object, using projection points are in
sites.prj <- sf::st_as_sf(sites.probs, coords = c("lon", "lat"), crs = PROJ.pts)
#then transform projection to match other map layers
sites.prj <- sf::st_transform(sites.prj, crs = PROJ)



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
  geom_sf(data = sf::st_jitter(sites.prj, 5000), aes(fill = taxclade), shape=21, colour="black", size = 4, alpha = 0.5) +
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

ggsave(filename = "figures/world_map_genetic_pts.pdf", width = 30, height = 15, units = c("cm"))






