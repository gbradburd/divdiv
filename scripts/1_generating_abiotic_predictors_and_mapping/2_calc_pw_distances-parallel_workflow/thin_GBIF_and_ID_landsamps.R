#idea: thin GBIF points if there are more than X present
# find GBIF and/or genetic sample points that are on land

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
# load data --------------------------------------------------

#get species name (to plug into function to get GBIF table)
species = run_name %>% gsub("bioprj_","",.) %>% str_split(., "_", simplify = T) %>% .[2] %>% gsub("-"," ",.)

#get x,y points for gbif and genetic samples
gbif <- read.csv(paste(indir,"/GBIFLocations_",species,".csv",sep="")) %>% 
  rename("x" = "decimalLongitude") %>% rename("y" = "decimalLatitude") %>%
  rename("run_acc_sra" = "sampid") %>% 
  dplyr::select(x,y,run_acc_sra) %>% arrange(run_acc_sra)

genetic <- read.delim(paste(indir,"/lat_long_table-",run_name,".txt",sep="")) %>% 
  rename("x" = "long") %>% rename("y" = "lat") %>% 
  dplyr::select(x,y,run_acc_sra) %>% arrange(run_acc_sra)



# ************************************************************
# ************************************************************
# thin (GBIF) datasets with >15,000 pts  ---------------------
# ! note ! this is slightly moving the GBIF pts - bc it is picking a midpoint of a raster cell
# so lat/long and sample names will not match between original pts downloaded from GBIF (e.g., GBIFLocations_Caretta caretta.csv)
# and GBIF pts post thinning (e.g., gbif-coordsinputformovinglandsamps-bioprj_PRJNA559677_Caretta-caretta.csv)

if ( nrow(gbif)>15000 ) {
  print(paste0("THINNING GBIF DATASET FROM ", nrow(gbif), " points to fewer points"))
  #from Bruce
  #rasterize data
  coord.nms<-c('x','y')
  dat <- gbif %>% dplyr::select(x,y)
  ncells<-setNames(c(625,625),coord.nms) #higher numbers here mean finer scale raster grid
  rans<-lapply(dat,range)
  breaks<-lapply(coord.nms,function(ii) seq(rans[[ii]][1],
                                            rans[[ii]][2]+(rans[[ii]][2]-rans[[ii]][1])/1e6,
                                            length.out=ncells[ii]+1))
  breaks<-setNames(data.frame(breaks),coord.nms)
  midpts<-data.frame(lapply(breaks,function(ii) (ii[-length(ii)]+ii[-1])/2))
  ras.dat<-lapply(colnames(dat),function(ii) findInterval(dat[,ii],
                                                          breaks[,ii]))
  ras.dat<-setNames(data.frame(ras.dat),coord.nms)
  ras.dat<-setNames(data.frame(lapply(ras.dat,as.factor)),coord.nms)
  #build "thinned" dist matrix based off unique cell IDs
  ras.dat$cellID<-droplevels(Reduce(':',ras.dat))
  unique.cells<-levels(ras.dat$cellID)
  n.unique.cells<-length(unique.cells)
  levels(ras.dat$cellID)<-seq_len(n.unique.cells)
  parsed.unique.cells<-do.call(rbind,lapply(strsplit(unique.cells,':'),as.numeric))
  for(i in seq_along(coord.nms)){
    parsed.unique.cells[,i]<-midpts[,i][parsed.unique.cells[,i]]
  }
  gbif <- parsed.unique.cells %>% as.data.frame() %>% rename("x"="V1", "y"="V2") %>% 
    mutate(run_acc_sra = paste0("fakethinnedgbif",row.names(.)))
  print(paste0("thinned dataset from ", nrow(dat), " to ", nrow(gbif), " points"))
} else {
  print("GBIF dataset has <15,000 points so no thinning applied")
}



# ************************************************************
# ************************************************************
# identify samples that are currently on land ----------------

print("****************************************")
print("identifying samples that are currently on land")

#merge genetic and gbif data into one df
points <- rbind(gbif %>% dplyr::select(x,y,run_acc_sra), genetic %>% dplyr::select(x,y,run_acc_sra))

#check if any points are on land
#(merging based on x/y coords doesn't work bc sometimes trailing 0's get added and then they don't match anymore)
landsamps <- data.frame("run_acc_sra"=NA, "lon"=NA, "lat"=NA, "depth"=NA)
for(i in 1:nrow(points)) {
  #print(paste0("starting point ", i, " of ", nrow(points)))
  apoint <- points[i,]
  sampid <- apoint %>% dplyr::select(run_acc_sra)
  apoint <- apoint %>% dplyr::select(x,y)
  apoint <- cbind(sampid, get.depth(mat = bathymap, x = apoint, locator = FALSE))
  landsamps <- rbind(landsamps, apoint)
}
landsamps <- landsamps %>% filter(depth > 0) %>% dplyr::select(run_acc_sra)

landsamps.anti <- data.frame("run_acc_sra"=NA, "lon"=NA, "lat"=NA, "depth"=NA)
for(i in 1:nrow(points)) {
  #print(paste0("starting point ", i, " of ", nrow(points)))
  apoint <- points[i,]
  sampid <- apoint %>% dplyr::select(run_acc_sra)
  apoint <- apoint %>% dplyr::select(x,y)
  apoint <- cbind(sampid, get.depth(mat = bathymap.antimeridian, x = apoint %>% mutate(x = ifelse(x < 0, x + 360, x)), locator = FALSE))
  apoint <- apoint %>% mutate(lon = ifelse(lon > 180, lon - 360, lon))
  landsamps.anti <- rbind(landsamps.anti, apoint)
}
landsamps.anti <- landsamps.anti %>% filter(depth > 0) %>% dplyr::select(run_acc_sra)

landsamps <- rbind(landsamps, landsamps.anti) %>% distinct()
landsamps <- points %>% filter(run_acc_sra %in% landsamps$run_acc_sra)
rm(landsamps.anti)
landsamps <- landsamps %>% mutate(ordered_id = 1:nrow(.))

print(paste0("there are this many points on land: ", nrow(landsamps)))




#*****************************************************************
#*****************************************************************
# save (thinned) list of GBIF points and landsamps ---------------

write.csv(gbif, paste(outdir,"/gbif-coordsinputformovinglandsamps-",run_name,".csv",sep=""))
write.csv(landsamps, paste(outdir,"/landsamps-gbif_and_genetic-",run_name,".csv",sep=""))
