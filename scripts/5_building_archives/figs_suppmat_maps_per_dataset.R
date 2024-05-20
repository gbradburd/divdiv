#idea: make map of species range and sampled genetic points - SPECIFIC FOR BEE

#load libraries
library(ggplot2)
library(dplyr)
library(tidyr)
#library(rnaturalearth)

rm(list=ls())
gc()

setwd("/Users/rachel/divdiv/")



#define variables

#define path to .csv's of lat/long after running marmap (and moving landpoints to water)
postcoordsdir = "data/abiotic/input_and_working/postmarmap_coords"

#define path to .csv's of lat/long for GBIF 
#(cleaned but no points moved or dropped for marmap aka the files that went into calcing n.ecoregions and mean/min/max lat)
gbifdir = "data/abiotic/input_and_working/speciesLocDat_Sep08_2021_withFilters-cleaned/"

#define path to save figs to
outdir = "suppmat_figures/maps"

#define path to WM files (so we can get genetic samples we kept in end)
wmpath = "data/popgen/input_and_working/ALL_r80_gendiv_data"

#define path to sample name keys
sampkeypath = "data/all_samplenamekeys"



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
meanlat <- read.csv("data/abiotic/mean_min_max_latitude_stats-wide.csv") %>% dplyr::select(-X)

#get final list of "in" datasets from master df
using <- read.csv("data/master_df.csv")
dataset_list <- using$run_name



#build PDF map for each dataset
for ( loop.iter in 1:length(dataset_list) ) {
  
  #get datasets loaded
  dataset <- dataset_list[loop.iter]
  
  run_name = dataset
  link.i = gsub("bioprj_","",run_name)
  print(paste0("starting to plot: ", run_name, "; ", loop.iter, " of ", length(dataset_list)))
  species = run_name %>% strsplit(., "_") %>% as.data.frame() %>% .[nrow(.),] %>% gsub("-"," ",.)
  
  #get GBIF and genetic points to map
  gbif.clean <- read.csv(paste(gbifdir,"/GBIFLocations_",species,".csv",sep=""))
  genetic.post <- read.csv(paste(postcoordsdir,"/genetic-coordspostmarmap-",run_name,".csv",sep=""))
  
  #get list of samples retained for popgen/gendiv calcs (so we only plot just genetic samples we used/kept in end)
  wmfile <- list.files(path = wmpath, pattern = run_name, full.names = TRUE)
  wmfile <- wmfile[which(grepl("initPars", wmfile)==T)]
  load(wmfile)
  if (dim(initPars$parHom)[1] != dim(initPars$parHom)[2]) {"error - matrix dims are not the same"}
  n_samples = dim(initPars$parHom)[1] #print this in title of final map PDF later in this script
  geneticsampskept <- rownames(initPars$parHom) %>% as.data.frame() %>% dplyr::rename("sampid_assigned_for_bioinf" = ".")
  #get sample name match key (to match SRR and sample0 name formats)
  sampkeyFile <- paste0(sampkeypath,"/samplenamekey-",run_name,".txt")
  sampkeyFile <- read.delim(sampkeyFile)
  geneticsampskept <- merge(geneticsampskept, 
                            sampkeyFile %>% dplyr::select(sampid_assigned_for_bioinf, run_acc_sra),
                            by = "sampid_assigned_for_bioinf")
  #filter to genetic lat/long df to keep just these samples (aka ones we kept in popgen)
  genetic.post <- genetic.post %>% filter(run_acc_sra %in% geneticsampskept$run_acc_sra)
  
  #merge into one dataset to make aes easier
  gbif.clean.plot <- gbif.clean %>% 
    dplyr::rename("run_acc_sra"="sampid", "x"="decimalLongitude", "y"="decimalLatitude") %>% 
    dplyr::select(run_acc_sra, x, y) %>% mutate(dataset = "gbif")
  genetic.post.plot <- genetic.post %>% 
    dplyr::select(run_acc_sra, x, y) %>% mutate(dataset = "genetic")
  points <- rbind(gbif.clean.plot, genetic.post.plot)
  
  #get range mean lat of GBIF points
  meanlat.sp <- meanlat %>% filter(link == link.i)
  meanlat.sp <- meanlat.sp$meanlat.gbif
  
  #set bounds for inset map aka zoomed in map of sample locations
  maxlat <- max(points$y) + 10
  minlat <- min(points$y) - 30
  maxlon <- max(points$x) + 10
  minlon <- min(points$x) - 30
  if (maxlat > 90) {maxlat = 90}
  if (minlat < -90) {minlat = -90}
  if (maxlon > 180) {maxlon = 180}
  if (minlon < -180) {minlon = -180}
  
  #world with box for inset -----
  inset.plot <- ggplot() + 
    #build map
    geom_sf(data = land, fill = "gray90", colour = "gray90", linewidth = 0.2) +
    geom_sf(data = coasts, fill = "transparent", colour = "gray80", linewidth = 0.2) +
    geom_sf(data = gls, fill = "white", colour = "gray80", linewidth = 0.3) +
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
    geom_sf(data = coasts, fill = "transparent", colour = "gray80", linewidth = 0.2) +
    geom_sf(data = gls, fill = "white", colour = "gray80", linewidth = 0.3) +
    #add mean lat of GBIF points
    geom_hline(yintercept = meanlat.sp, colour = "black", linewidth = 0.4, linetype = 2) +
    #add points (last two values of hex aka 60 specify alpha/opacity)
    geom_point(data = points, aes(x=x, y=y, colour=dataset), size = 1) +
    scale_colour_manual(values = c("#0000FF60", "magenta"), 
                        labels = c("Occurrence reports (GBIF)", "Genetic samples (INSDC)")) +
    #manually make a legend for points and lines
    labs(title = paste0("Dataset: ", run_name),
         subtitle = paste0("N = ",n_samples," genetic samples; dashed line shows mean range latitude (occurrence records)"),
         x = "Longitude",
         y = "Latitude") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = c(0.88, 0.08),
          legend.background = element_blank(),
          legend.box.background = element_rect(fill = "white"),
          legend.margin = margin(0, 10, 5, 10),
          legend.key = element_blank(),
          legend.title = element_blank()
          ) +
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
     gbif.clean,genetic.post,
     geneticsampskept,initPars,plot,
     sampkeyFile,meanlat.sp,n_samples,wmfile,
     gbif.clean.plot,genetic.post.plot,points,
     main.plot,inset.plot,
     link.i,dataset)
  
}


#end



