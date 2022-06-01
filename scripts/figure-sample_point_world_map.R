#idea: make a map of the world and plot NCBI current winners marine data on it

#load libraries -------
library(rgdal)      
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)

rm(list=ls())
gc()



# load files for world map --------

load("data/earth_map_objects.RData")

## This will load 5 R objects:
##   xbl.X & lbl.Y are two data.frames that contain labels for graticule lines
##   NE_box is a SpatialPolygonsDataFrame object and represents a bounding box for Earth 
##   NE_countries is a SpatialPolygonsDataFrame object representing countries 
##   NE_graticules is a SpatialLinesDataFrame object that represents 10 dg latitude lines and 20 dg longitude lines
##   NOTE: data downloaded from http://www.naturalearthdata.com/


### CHOOSE YOUR OWN ADVENTURE (2 options) ###

# # *** OPTION 1 *** : build inputs from scratch aka raw files (works for RT's desktop) ----------
# 
# # read in lat/long points we want to plot on map
# #datasets to plot
# #get master spreadsheet from Drive
# using <- googledrive::drive_ls(path = "working_datasheets", pattern = "working_list_marine_projects_with_10indivs-12-4-2020")
# using <- googlesheets4::range_read(using, sheet = 1) %>% as.data.frame() %>% mutate(run_name = paste("bioprj_",link,sep=""))
# #keep cols we want and only datasets that are "in"
# using <- using %>% dplyr::select(organism_biosamp,run_name,link,keepinrunning_YN) %>%
#   filter(grepl("yes|Yes",keepinrunning_YN)) %>% mutate(y = "full list")
# 
# #get lat/longs
# datasets <- list.files(path = "/Users/rachel/Desktop/DivDiv/mapping_gbif_and_genetic_sample_points/lat_long_tables_per_dataset",
#                        pattern = "lat_long_table", full.names = TRUE)
# sites <- data.frame("link"=NA, "run_acc_sra"=NA, "lat"=NA, "long"=NA, "coordinateUncertaintyInMeters"=NA, "origin"=NA)
# for ( d in datasets ) {
#   out <- read.delim(d)
#   sites <- rbind(sites,out)
# }
# sites <- sites %>% filter(is.na(link)==F)
# write.csv(sites, "data/all_divdiv_lat_long_tables_for_genetic_samps.csv", row.names = FALSE)
# 
# #get list of all datasets
# sites.all <- sites %>% filter(link %in% using$link)
# 
# #get list of datasets with wishart values
# havewishart <- read.csv("data/master_df.csv") %>% filter(is.na(s.wish)==FALSE) %>% dplyr::select(link) %>% mutate(link = gsub("bioprj_","",link))
# sites.pg <- sites %>% filter(link %in% havewishart$link)



# *** OPTION 2 *** : load pre=saved files (works for anyone with /divdiv repo) ---------

#get all divdiv lat/longs
sites <- read.csv("data/all_divdiv_lat_long_tables_for_genetic_samps.csv")

#add clade groups on for coloring

#set taxonomic/clade colors (to match phy: figure-phylogeny_of_species.R)
cladeCols <- RColorBrewer::brewer.pal(10,"Paired")

load("data/phylo/divdiv_phy_from_timetreebeta5.Robj")

#bony fishes
des <- ape::getMRCA(phy,c("Anguilla rostrata","Sebastes diaconus")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)] # get rid of internal nodes and just keep tip "nodes"
df1 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "bony fishes")
#chondrichthyes
des <- ape::getMRCA(phy,c("Bathyraja parmifera","Sphyrna tiburo")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)]
df2 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "chondrichthyes")
#sauropsida
des <- ape::getMRCA(phy,c("Pygoscelis papua","Caretta caretta")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)]
df3 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "sauropsida")
#mammals
des <- ape::getMRCA(phy,c("Phocoena sinus","Halichoerus grypus atlantica")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)]
df4 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "mammals")
#echinoderms
des <- ape::getMRCA(phy,c("Apostichopus californicus","Pisaster ochraceus")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)]
df5 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "echinoderms")
#molluscs
des <- ape::getMRCA(phy,c("Bathymodiolus platifrons","Nautilus pompilius")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)]
df6 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "molluscs")
#crustacea
des <- ape::getMRCA(phy,c("Callinectes sapidus","Isocladus armatus")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)]
df7 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "crustacea")
#cnidarians
des <- ape::getMRCA(phy,c("Galaxea horrescens","Ectopleura larynx")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)]
df8 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "cnidarians")
#vascular plants
des <- ape::getMRCA(phy,c("Rhizophora mangle","Laguncularia racemosa")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)]
df9 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "vascular plants")
#ochrophyta
des <- ape::getMRCA(phy,c("Fucus vesiculosus","Sargassum muticum")) %>% phytools::getDescendants(phy,.)
des <- des[des<ape::Ntip(phy)]
df10 <- phy$tip.label[des] %>% as.data.frame() %>% rename("species"=".") %>% mutate(taxclade = "ochrophyta")
taxcladelbs <- rbind(df1,df2,df3,df4,df5,df6,df7,df8,df9,df10)

#get list of datasets with wishart values
havewishart <- read.csv("data/master_df.csv") %>% filter(is.na(s.wish)==FALSE) %>% dplyr::select(link) %>% mutate(link = gsub("bioprj_","",link))
sites.pg <- sites %>% filter(link %in% havewishart$link)



# set up projections for map and our points --------------
## give the PORJ.4 string for Eckert IV projection
PROJ <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" 
## or use the short form "+proj=eck4"

## re-project the shapefiles (warnings about CRS object has no comment okay)
NE_countries.prj  <- spTransform(NE_countries, CRSobj = PROJ)
NE_graticules.prj <- spTransform(NE_graticules, CRSobj = PROJ)
NE_box.prj        <- spTransform(NE_box, CRSobj = PROJ)

## project long-lat coordinates columns for data frames 
## (two extra columns with projected XY are created)
prj.coord <- project(cbind(lbl.Y$lon, lbl.Y$lat), proj = PROJ)
lbl.Y.prj <- cbind(prj.coord, lbl.Y)
names(lbl.Y.prj)[1:2] <- c("X.prj","Y.prj")

prj.coord <- project(cbind(lbl.X$lon, lbl.X$lat), proj = PROJ)
lbl.X.prj <- cbind(prj.coord, lbl.X)
names(lbl.X.prj)[1:2] <- c("X.prj","Y.prj")

#project sites to map
sites.df  = sites.pg
prj.coord <- project(cbind(sites.df$long, sites.df$lat), proj = PROJ)
sites.df.prj <- cbind(prj.coord, sites.df)
names(sites.df.prj)[1:2] <- c("X.prj","Y.prj")
#add tax coloring on
sites.df.prj <- merge(sites.df.prj %>% tidyr::separate(., link, into = c("garbage","species.temp"), sep="_") %>% mutate(species = gsub("-"," ",species.temp)) %>% dplyr::select(-species.temp,-garbage), 
                      taxcladelbs, 
                      by = "species", all.x = T)

jitter <- position_jitter(width = 10, height = 10)

# plot ----------------------
ggplot() +
  
  coord_fixed(ratio = 1) +
  
  ## add projected bounding box, blue for ocean (#abd9e9)
  geom_polygon(data = NE_box.prj, 
               aes(x = long, y = lat), 
               colour = "black", fill = "white", size = .25) +
  
  ## add graticules
  geom_path(data = NE_graticules.prj,
            aes(long, lat, group = group),
            linetype = "dotted", colour = "grey70", size = .25) +
  
  ## add projected countries, a bit darker grey for terrestrial polygons
  geom_polygon(data = NE_countries.prj, 
               aes(long,lat, group = group), 
               colour = "gray50", fill = "gray90", size = .25) +
  
  ## add locations (points)
  geom_point(data = sites.df.prj, aes(x = X.prj, y = Y.prj, fill = taxclade), shape=21, colour="black", size = 4, alpha = 0.5) +
  scale_fill_manual(values = c("#A6CEE3","#FF7F00","#FDBF6F","#B2DF8A","#E31A1C","#1F78B4","#CAB2D6")) +
  
  ## Set empty theme
  theme_void() + # remove the default background, gridlines & default gray color around legend's symbols
  theme(
    legend.position = "none",
    #legend.text = element_text(size = 20),
    #legend.title = element_blank(),
    #plot.margin = unit(c(t=0, r=3, b=0, l=0), unit="cm"), # adjust margins
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white")
  )

ggsave(filename = "figures/world_map_wishart_pts.pdf", width = 30, height = 15, units = c("cm"))







#graveyard ------

#read in lat/long of points we want to plot onto map - previous/old way
sites.raw <- read.csv("marine/current_winners/current_winners_list-sequence_level-11-21-19.csv", header = T, stringsAsFactors = F) %>% 
  dplyr::select(-X)
str(sites.raw)
sites <- sites.raw %>% dplyr::select(organism_biosamp,family_worms,lat_long_biosamp,biosample_acc_sra,project_acc_bioprj)
dim(sites)
#keep just records with lat/long (so we can plot them)
sites <- sites %>% filter(is.na(lat_long_biosamp) == F) %>% filter(lat_long_biosamp != "unknown") %>% filter(lat_long_biosamp != "not given in NCBI")
dim(sites)
#convert lat/long from N/S/E/W to +/- (so they will plot nice)
sites <- sites %>% separate(.,lat_long_biosamp,into = c("lat","lat_dir","lon","lon_dir"), sep = " ", remove = F)
sites %>% group_by(lat_dir) %>% summarise(n=n())
sites <- sites %>% filter(lat_dir %in% c("N","S"))
sites %>% group_by(lat_dir) %>% summarise(n=n())
sites %>% group_by(lon_dir) %>% summarise(n=n())
sites$lat <- as.numeric(sites$lat)
sites$lon <- as.numeric(sites$lon)
#north and east are + ; south and west are - 
sites <- sites %>% mutate(lat = ifelse(lat_dir == "N", lat, lat*-1)) %>% 
  mutate(lon = ifelse(lon_dir == "E", lon, lon*-1))
sites <- sites %>% dplyr::select(lat,lon,everything())




#assign tax categories
df.tax <- df.1 %>% mutate(tax_category = ifelse(kingdom_ncbitaxonomy.via.taxize == "Fungi", "Fungi", ".")) %>% 
  mutate(tax_category = ifelse(kingdom_ncbitaxonomy.via.taxize == "Viridiplantae", "Other Plants", tax_category)) %>%
  mutate(tax_category = ifelse(class_ncbitaxonomy.via.taxize == "Cycadopsida", "Gymnosperms", tax_category)) %>% 
  mutate(tax_category = ifelse(class_ncbitaxonomy.via.taxize == "Gnetopsida", "Gymnosperms", tax_category)) %>% 
  mutate(tax_category = ifelse(class_ncbitaxonomy.via.taxize == "Pinopsida", "Gymnosperms", tax_category)) %>%
  mutate(tax_category = ifelse(class_ncbitaxonomy.via.taxize == "Magnoliopsida", "Angiosperms", tax_category)) %>% 
  mutate(tax_category = ifelse(class_ncbitaxonomy.via.taxize == "Aves", "Birds", tax_category)) %>% 
  mutate(tax_category = ifelse(class_ncbitaxonomy.via.taxize == "Mammalia", "Mammals", tax_category)) %>% 
  mutate(tax_category = ifelse(class_ncbitaxonomy.via.taxize %in% c("Actinopteri","Chondrichthyes","Coelacanthimorpha"), "Fish", tax_category)) %>% 
  mutate(tax_category = ifelse(class_ncbitaxonomy.via.taxize %in% c("Lepidosauria","Amphibia"), "Amphibians and Reptiles", tax_category)) %>%
  mutate(tax_category = ifelse(phylum_ncbitaxonomy.via.taxize == "Arthropoda", "Arthropods", tax_category)) %>% 
  mutate(tax_category = ifelse(phylum_ncbitaxonomy.via.taxize == "Mollusca", "Mollusks", tax_category)) %>%
  mutate(tax_category = ifelse(division_taxonomy == "Mammals" & is.na(tax_category)==T, "Mammals", tax_category)) %>% 
  mutate(tax_category = ifelse(division_taxonomy == "Vertebrates" & tax_category==".", "Other Vertebrates", tax_category)) %>% 
  mutate(tax_category = ifelse(division_taxonomy == "Invertebrates" & tax_category==".", "Invertebrates", tax_category)) %>% mutate(tax_category = ifelse(kingdom_ncbitaxonomy.via.taxize == "Fungi" & is.na(tax_category)==T, "Fungi", tax_category)) %>% 
  mutate(tax_category = ifelse(division_taxonomy == "Vertebrates" & is.na(tax_category)==T, "Other Vertebrates", tax_category)) %>% 
  mutate(tax_category = ifelse(division_taxonomy == "Invertebrates" & is.na(tax_category)==T, "Invertebrates", tax_category)) %>% 
  mutate(tax_category = ifelse(division_taxonomy == "Plants and Fungi" & is.na(tax_category)==T, "Other", tax_category)) %>% 
  mutate(tax_category = ifelse(tax_category == "Invertebrates", "Other Invertebrates", tax_category)) %>% 
  dplyr::select(link,tax_category,class_ncbitaxonomy.via.taxize,phylum_ncbitaxonomy.via.taxize,order_ncbitaxonomy.via.taxize,family_ncbitaxonomy.via.taxize,contains("worm")) %>% distinct()

sites <- merge(df, df.tax, by = "link", all.x = T) %>% dplyr::select(lat,lon,everything()) %>% filter(is.na(link)==F)

sites$tax_category[sites$link=="PRJNA311981_Engraulis-encrasicolus"] = "Fish"
sites$tax_category[sites$link=="PRJNA394157_Exaiptasia-pallida"] = "Other Invertebrates"
sites$tax_category[sites$link=="PRJNA473288_Ectocarpus-siliculosus"] = "Other"
sites$tax_category[sites$link=="PRJNA506239_Crassostrea-virginica"] = "Mollusks"
sites$tax_category[sites$tax_category=="Other Invertebrates"] = "Other"
sites$tax_category[sites$tax_category=="Other Vertebrates"] = "Other"
