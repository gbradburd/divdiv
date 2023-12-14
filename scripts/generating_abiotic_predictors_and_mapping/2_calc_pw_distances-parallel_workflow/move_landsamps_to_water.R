#idea: use our custom function that uses some marmap functions to relocate sample points on land 
#to nearest point in ocean

#load libraries -----
suppressMessages(library(marmap, quietly = TRUE, warn.conflicts = FALSE))
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyr, quietly = TRUE, warn.conflicts = FALSE)
library(ggplot2, quietly = TRUE, warn.conflicts = FALSE)
library(stringr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyfast, quietly = TRUE, warn.conflicts = FALSE)
suppressMessages(library(gdistance, quietly = TRUE, warn.conflicts = FALSE)) #not needed directly, but loaded when loading bathy objs
library(raster, quietly = TRUE, warn.conflicts = FALSE) #not needed directly, but loaded when loading bathy objs
library(sp, quietly = TRUE, warn.conflicts = FALSE) #not needed directly, but loaded when loading bathy objs


#print session info and args for HPCC metadata/log files
print("----------------------------------------------------------------------------------------")
print("----------------------------------------------------------------------------------------")
print("we just opened R and loaded all libraries")
print("session info:")
sessionInfo()
#get and print a vector of all arguments passed from "-export" in submit file
args <- commandArgs(trailingOnly = TRUE)
print("----------------------------------------------------------------------------------------")
print("----------------------------------------------------------------------------------------")
print("we are printing the contents of the vector args to show all variables passed from hpcc to R enviro")
print("arguments passed:")
cat(args, sep = "\n")

#define some variables (pulled in from bash script)
#run_name = "bioprj_PRJNA553831_Acanthephyra-purpurea"
#indir = getwd()
#outdir = getwd()
run_name = args[1] #dataset name
indir = args[3] #indir
outdir = args[4] #outdir
start.row = as.numeric(args[5]) #row of landsamps to start on
stop.row = as.numeric(args[6]) #row of landsamps to end on
#useantimeridian = "FALSE"

#source our functions
source(paste0(indir,"/map_functions.R"))

#load (already fetched) bathymetric maps of world (one -180 to 180 longitude and another 0 - 360 longitude)
#bathymap <- getNOAA.bathy(lat1 = -90, lon1 =-180, lat2 = 90, lon2 = 180, resolution = 4, keep = FALSE, antimeridian = FALSE, path = NULL)
#save(bathymap, file = "/Users/rachel/Desktop/DivDiv/mapping_GBIF_and_genetic_sample_points/world_NOAA_bathy_res4.Robj")
#load bathymap of world
load(paste(indir,"/world_NOAA_bathy_res4.Robj",sep=""))

#bathymap.antimeridian <- getNOAA.bathy(lon1 = -180, lon2 = 180, lat1 = -90, lat2 = 90, resolution = 4, keep = FALSE, antimeridian = TRUE, path = NULL)
#save(bathymap.antimeridian, file = "/Users/rachel/Desktop/DivDiv/mapping_GBIF_and_genetic_sample_points/world_NOAA_bathy_res4-antimeridian.Robj")
#load antimeridian bathymap of world
load(paste(indir,"/world_NOAA_bathy_res4-antimeridian.Robj",sep=""))


#define resistance matrix
#min.depth=0, max.depth=NULL says just avoid land
#resistance.matrix <- trans.mat(bathymap, min.depth=0, max.depth=NULL)
#save(resistance.matrix, file = "/Users/rachel/Desktop/DivDiv/mapping_GBIF_and_genetic_sample_points/world_NOAA_bathy_res4-resistancematrix-mindepth0.Robj")
#load resistance matrix of whole world
load(paste(indir,"/world_NOAA_bathy_res4-resistancematrix-mindepth0.Robj",sep=""))

#resistance.matrix.antimeridian <- trans.mat(bathymap.antimeridian, min.depth=0, max.depth=NULL)
#save(resistance.matrix.antimeridian, file = "/Users/rachel/Desktop/DivDiv/mapping_GBIF_and_genetic_sample_points/world_NOAA_bathy_res4-resistancematrix-mindepth0-antimeridian.Robj")
#load resistance matrix of whole world
load(paste(indir,"/world_NOAA_bathy_res4-resistancematrix-mindepth0-antimeridian.Robj",sep=""))


#print objects loaded in R (so we can see functions and variables were correctly loaded)
print("----------------------------------------------------------------------------------------")
print("----------------------------------------------------------------------------------------")
print("so far we have loaded libraries and printed args")
print("objects currently loaded in R:")
ls()



# ************************************************************
# ************************************************************
# load data -----------------------

#get x,y points for gbif and genetic samples on land currently
landsamps <- read.csv(paste(indir,"/landsamps-gbif_and_genetic-",run_name,".csv",sep="")) %>% 
  arrange(ordered_id) %>% dplyr::select(-ordered_id,-X)

#filter landsamps to run just samples for this job
if (nrow(landsamps)>0) {
  landsamps.sub <- landsamps[start.row:stop.row,]
  print(paste0("Moving points ", start.row, " thru ", stop.row, " from land to water"))
}

if (nrow(landsamps)==0) {
  print("There were no points on land to move, ending job")
}


  
# ************************************************************
# ************************************************************
# move land points to water (if needed) -----------------------

if (nrow(landsamps)>0) {
  #find points on land and move them to water
  landsamps.sub <- move_land_points_to_water(landsamps = landsamps.sub, 
                                         bathymap = bathymap, bathymap.antimeridian = bathymap.antimeridian,
                                         starting.search.radius = 1, bathyres = 4,
                                         resistance.matrix = resistance.matrix)
  
  #check if any landpoints are still on land (excluding dropped points that broke dist2isobath)
  stillonland <- get.depth(mat = bathymap, x = landsamps.sub$end.lon, y = landsamps.sub$end.lat, locator = F) %>% cbind(., landsamps.sub) %>% filter(depth > 0) %>% mutate(check = ifelse(initial.end.lon == 0 & initial.end.lat == 0 & end.lon == 0 & end.lat == 0, "toss","keep")) %>% filter(check == "keep") %>% dplyr::select(-check)
  stillonland.anti <- get.depth(mat = bathymap.antimeridian, x = ifelse(landsamps.sub$end.lon < 0, landsamps.sub$end.lon + 360, landsamps.sub$end.lon), y = landsamps.sub$end.lat, locator = F) %>% cbind(., landsamps.sub) %>% filter(depth > 0) %>% mutate(check = ifelse(initial.end.lon == 0 & initial.end.lat == 0 & end.lon == 0 & end.lat == 0, "toss","keep")) %>% filter(check == "keep") %>% dplyr::select(-check)
  stillonland <- rbind(stillonland, stillonland.anti)
  if (nrow(stillonland)==0) {
    print("Good - All points have a depth =< 0")
  }  else {
    print("ERROR - some points still have a positive depth")
    print(stillonland)
  }
  
}
  
 
 
#*****************************************************************
#*****************************************************************
# save -------

out <- format(landsamps.sub, digits = 22)
write.csv(out, paste(outdir,"/moved_landsamps.",start.row,"-",stop.row,".",run_name,".csv",sep=""), row.names = FALSE)

  
  
  
