#idea: make map of species range and sampled genetic points - SPECIFIC FOR BEE

#load libraries
library(ggplot2)
library(dplyr)
library(tidyr)
#library(rnaturalearth)

rm(list=ls())
gc()



#define path to .csv's of lat/long after running marmap (and moving landpoints to water)
postcoordsdir = "/Users/rachel/divdiv/data/abiotic/input_and_working/postmarmap_coords"

#define path to save figs to
outdir = "/Users/rachel/divdiv/suppmat_figures/maps"

#get world map outlines
coasts <- rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf")
land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

#get Great Lakes outlines
lakes <- rnaturalearth::ne_download(scale="medium", category = 'physical', type = "lakes", returnclass = "sf")
#subset world lake outlines to just GLs
gls <- c("Lake Superior", "Lake Huron",
         "Lake Michigan", "Lake Erie", "Lake Ontario")
gls <- subset(lakes, name %in% gls)
#how does it look?
ggplot() +
  geom_sf(data = gls)

#get df of range mean lats
meanlat <- read.csv("/Users/rachel/divdiv/data/abiotic/mean_min_max_latitude_stats-wide.csv") %>% dplyr::select(-X)

#get final list of "in" datasets from master df
using <- read.csv("/Users/rachel/divdiv/data/master_df.csv")
dataset_list <- using$link


for ( loop.iter in 1:length(dataset_list) ) {
  
  #get datasets loaded
  dataset <- dataset_list[loop.iter]
  
  run_name = dataset
  link.i = gsub("bioprj_","",run_name)
  print(paste0("starting to plot: ", run_name, "; ", loop.iter, " of ", length(dataset_list)))
  species = run_name %>% strsplit(., "_") %>% as.data.frame() %>% .[nrow(.),] %>% gsub("-"," ",.)
  
  gbif.post <- read.csv(paste(postcoordsdir,"/gbif-coordspostmarmap-",run_name,".csv",sep=""))
  genetic.post <- read.csv(paste(postcoordsdir,"/genetic-coordspostmarmap-",run_name,".csv",sep=""))
  
  #get range mean lat of GBIF points
  meanlat.sp <- meanlat %>% filter(link == link.i)
  meanlat.sp <- meanlat.sp$meanlat.gbif
  
  #set bounds for inset map aka zoomed in map of sample locations
  maxlat <- max(gbif.post$y, genetic.post$y) + 10
  minlat <- min(gbif.post$y, genetic.post$y) - 30
  maxlon <- max(gbif.post$x, genetic.post$x) + 10
  minlon <- min(gbif.post$x, genetic.post$x) - 30
  if (maxlat > 90) {maxlat = 90}
  if (minlat < -90) {minlat = -90}
  if (maxlon > 180) {maxlon = 180}
  if (minlon < -180) {minlon = -180}
  
  #world with box for inset -----
  inset.plot <- ggplot() + 
    #build map
    geom_sf(data = land, fill = "gray90", colour = "gray90", linewidth = 0.2) +
    geom_sf(data = coasts, fill = "transparent", colour = "gray70", linewidth = 0.2) +
    geom_sf(data = gls, fill = "white", colour = "gray70", linewidth = 0.3) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    #add zoom box
    geom_rect(aes(xmin = minlon, xmax = maxlon, ymin = minlat, ymax = maxlat), colour = "black", fill = "transparent") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.background = element_rect(fill = "transparent", colour = NA),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    guides(fill = "none") +
    labs(x = "Longitude",
         y = "Latitude") +
    coord_sf()
  
  # points zoomed --------
  main.plot <- ggplot() + 
    #build map
    geom_sf(data = land, fill = "gray90", colour = "gray90", linewidth = 0.2) +
    geom_sf(data = coasts, fill = "transparent", colour = "gray70", linewidth = 0.2) +
    geom_sf(data = gls, fill = "white", colour = "gray70", linewidth = 0.3) +
    #add mean lat of GBIF points
    geom_hline(yintercept = meanlat.sp, colour = "black", linewidth = 0.4, linetype = 2) +
    #add points
    geom_point(data = gbif.post, aes(x=x, y=y, color = "blue"), alpha = 0.35, size = 1) +
    geom_point(data = genetic.post, aes(x=x, y=y, color = "magenta"), alpha = 1, size = 1) +
    #manually make a legend for points and lines
    scale_color_identity(guide = "legend",
                         name = "",
                         breaks = c("blue", "magenta"),
                         labels = c("Occurrence reports (GBIF)", "Genetic samples (INSDC)")) +
    labs(title = paste0("Dataset: ", run_name),
         subtitle = paste0("N = ",nrow(genetic.post)," genetic samples"),
         x = "Longitude",
         y = "Latitude") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = c(0.88, 0.08),
          legend.background = element_blank()) +
    guides(fill = "none") +
    coord_sf(xlim = c(minlon,maxlon),
             ylim = c(minlat,maxlat))
  
  # put together and save
  pdf(file = paste(outdir,"/suppmat-map_of_occurrence_and_genetic_samples-",run_name,".pdf", sep = ""), width = 11, height = 8.5)
  plot <- cowplot::ggdraw() +
    cowplot::draw_plot(main.plot) +
    #draw_plot(inset.plot, x = 0.04, y = 0.734, width = 0.25, height = 0.25)
    cowplot::draw_plot(inset.plot, x = 0.02, y = 0.01, width = 0.25, height = 0.25)
  print(plot)
  dev.off()
  rm(maxlat,maxlon,minlat,minlon,run_name,species,
     bathymapZOOMED,bathymapZOOMED.simple,
     gbif.post,genetic.post,
     main.plot,inset.plot)
  
}


#end



