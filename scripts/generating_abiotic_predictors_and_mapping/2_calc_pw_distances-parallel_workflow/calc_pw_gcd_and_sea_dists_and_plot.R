#idea: finish rest of marmap workflow aka 
#check that all samples are on land now,
# calc distances between all pairs of GBIF points and all pairs of genetic points
# do great circle distance and 
# distance via sea aka avoiding land
# get map with lat (y), long (x), and depth (z) info from NOAA,
# plot GBIF occurrence points and sampled genetic (NCBI) points on it



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
# load data -----------------------

#get species name (to plug into function to get GBIF table)
species = run_name %>% gsub("bioprj_","",.) %>% str_split(., "_", simplify = T) %>% .[2] %>% gsub("-"," ",.)

#get x,y points for gbif and genetic samples
gbif <- read.csv(paste(indir,"/gbif-coordsinputformovinglandsamps-",run_name,".csv",sep="")) %>% 
  dplyr::select(-X) %>% arrange(run_acc_sra)

genetic <- read.delim(paste(indir,"/lat_long_table-",run_name,".txt",sep="")) %>% 
  rename("x" = "long") %>% rename("y" = "lat") %>% 
  dplyr::select(x,y,run_acc_sra) %>% arrange(run_acc_sra)



# ************************************************************
# ************************************************************
# gather and merge all landsamp dataset together -------------

list.landsamp.files <- list.files(path = indir, pattern = "moved_landsamps", full.names = TRUE)
print("list.landsamp.file contains:")
print(list.landsamp.files)
landsamps <- data.frame("run_acc_sra"=NA, "x"=NA, "y"=NA, "initial.end.lon"=NA, "initial.end.lat"=NA,
                        "end.lon"=NA, "end.lat"=NA, "total_dist_moved_km"=NA, "wiggle_dist_moved_km"=NA)
for (file in list.landsamp.files) {
  temp <- read.csv(file, header = T)
  landsamps <- rbind(landsamps, temp)
}
landsamps <- landsamps %>% filter(is.na(run_acc_sra)==FALSE)

#merge back onto gbif and genetic original dfs
gbif <- merge(gbif, landsamps %>% dplyr::select(-x,-y), by = "run_acc_sra", all.x = T)
genetic <- merge(genetic, landsamps %>% dplyr::select(-x,-y), by = "run_acc_sra", all.x = T)
#drop any points we weren't able to move due to issues with dist2isobath
gbif <- gbif %>% mutate(check = ifelse(initial.end.lon != 0 & initial.end.lat != 0 & end.lon != 0 & end.lat != 0 | is.na(initial.end.lon) == T, "keep","toss")) %>% filter(check == "keep") %>% dplyr::select(-check)
genetic <- genetic %>% mutate(check = ifelse(initial.end.lon != 0 & initial.end.lat != 0 & end.lon != 0 & end.lat != 0 | is.na(initial.end.lon) == T, "keep","toss")) %>% filter(check == "keep") %>% dplyr::select(-check)
landsamps <- landsamps %>% mutate(check = ifelse(initial.end.lon != 0 & initial.end.lat != 0 & end.lon != 0 & end.lat != 0 | is.na(initial.end.lon) == T, "keep","toss")) %>% filter(check == "keep") %>% dplyr::select(-check)

#sub in new water locations
gbif <- gbif %>% rename("x.orig" = "x") %>% rename("y.orig" = "y") %>% 
  mutate(x = ifelse(is.na(end.lon)==F, end.lon, x.orig)) %>% 
  mutate(y = ifelse(is.na(end.lat)==F, end.lat, y.orig))
genetic <- genetic %>% rename("x.orig" = "x") %>% rename("y.orig" = "y") %>% 
  mutate(x = ifelse(is.na(end.lon)==F, end.lon, x.orig)) %>% 
  mutate(y = ifelse(is.na(end.lat)==F, end.lat, y.orig))

gbif <- gbif %>% arrange(run_acc_sra)
genetic <- genetic %>% arrange(run_acc_sra)



# *********************************************************************************************
# *********************************************************************************************
# calculate pairwise great circle distances between points moved from land to water ----------

#returns lower triangle distance matrix in km btw points from decimal degrees long,lat df
#will work if longitudes are bewteen -180 to 180 or 0 to 360

print("****************************************")
print("calcing pairwise great circle distances between points moved to water")

finallandsamps <- rbind(gbif %>% filter(is.na(end.lon)==F) %>% dplyr::select(x, y, run_acc_sra), 
                        genetic %>% filter(is.na(end.lon)==F) %>% dplyr::select(x, y, run_acc_sra))

# moved landsamps
quantile = 100
gcd.landsamps <- pwGCD(points = finallandsamps, quantile = quantile)
names(gcd.landsamps)[names(gcd.landsamps) == "max.quant.gcd"] <- paste("max.",quantile,".gcd",".landsamps", sep="")
names(gcd.landsamps)[names(gcd.landsamps) == "pw.gcd"] <- paste("pw.gcd",".landsamps", sep="")
names(gcd.landsamps)[names(gcd.landsamps) == "max.quant.gcdpr"] <- paste("max.",quantile,".gcdpr",".landsamps", sep="")

#one final check that all points are in water now (excluding points dropped bc they broke dist2isobath)
points <- rbind(gbif %>% dplyr::select(x,y,run_acc_sra), genetic %>% dplyr::select(x,y,run_acc_sra))
stillonland <- get.depth(mat = bathymap, x = points$x, y = points$y, locator = F) %>% cbind(., points$run_acc_sra) %>% filter(depth > 0)
stillonland.anti <- get.depth(mat = bathymap.antimeridian, x = ifelse(points$x < 0, points$x + 360, points$x), y = points$y, locator = F) %>% cbind(., points$run_acc_sra) %>% filter(depth > 0)
stillonland <- rbind(stillonland, stillonland.anti)
if (nrow(stillonland)==0) {
  print("Good - All points have a depth =< 0")
}  else {
  print("ERROR - some points still have a positive depth")
  print(stillonland)
}



# ************************************************************
# ************************************************************
# calculate pairwise great circle distances -----------

#returns lower triangle distance matrix in km btw points from decimal degrees long,lat df
#will work if longitudes are bewteen -180 to 180 or 0 to 360

print("****************************************")
print("calcing pairwise great circle distances")

gbif <- gbif %>% dplyr::select(x,y,run_acc_sra,everything())
genetic <- genetic %>% dplyr::select(x,y,run_acc_sra,everything())

# GBIF occurrences - X% quantile
quantile = 100
gcd.gbif <- pwGCD(points = gbif, quantile = quantile)
pw.gcd.gbif <- gcd.gbif$pw.gcd #save local temp copy for plotting later, but don't write out into saved list
gcd.gbif <- within(gcd.gbif, rm(pw.gcd)) #don't save pw dist matrix for gbif points (but we do for genetic pts)
names(gcd.gbif)[names(gcd.gbif) == "max.quant.gcd"] <- paste("max.",quantile,".gcd",".gbif", sep="")
names(gcd.gbif)[names(gcd.gbif) == "max.quant.gcdpr"] <- paste("max.",quantile,".gcdpr",".gbif", sep="")
gcd.gbif.1 <- gcd.gbif

quantile = 95
gcd.gbif <- pwGCD(points = gbif, quantile = quantile)
pw.gcd.gbif <- gcd.gbif$pw.gcd #save local temp copy for plotting later, but don't write out into saved list
gcd.gbif <- within(gcd.gbif, rm(pw.gcd)) #don't save pw dist matrix for gbif points (but we do for genetic pts)
names(gcd.gbif)[names(gcd.gbif) == "max.quant.gcd"] <- paste("max.",quantile,".gcd",".gbif", sep="")
names(gcd.gbif)[names(gcd.gbif) == "max.quant.gcdpr"] <- paste("max.",quantile,".gcdpr",".gbif", sep="")
gcd.gbif.2 <- gcd.gbif

# genetic samples
quantile = 100
gcd.genetic <- pwGCD(points = genetic, quantile = quantile)
names(gcd.genetic)[names(gcd.genetic) == "max.quant.gcd"] <- paste("max.",quantile,".gcd",".genetic", sep="")
names(gcd.genetic)[names(gcd.genetic) == "pw.gcd"] <- paste("pw.gcd",".genetic", sep="")
names(gcd.genetic)[names(gcd.genetic) == "max.quant.gcdpr"] <- paste("max.",quantile,".gcdpr",".genetic", sep="")



# ***************************************************************************************
# ***************************************************************************************
# calculate pairwise sea distances (aka shortest path btw 2 pts avoiding land) -----------

#returns lower triangle distance matrix in km btw points from decimal degrees long,lat df

print("****************************************")
print("calcing pairwise via sea distances")

# GBIF occurrences - X% quantile
quantile = 100
sea.gbif <- pwseadist(points = gbif, resistance.matrix = resistance.matrix, quantile = quantile, resistance.matrix.antimeridian = resistance.matrix.antimeridian)
pw.seadist.gbif <- sea.gbif$pw.seadist
sea.gbif <- within(sea.gbif, rm(pw.seadist)) #don't save pw dist matrix for gbif points (but we do for genetic pts)
names(sea.gbif)[names(sea.gbif) == "max.quant.sea"] <- paste("max.",quantile,".sea",".gbif", sep="")
names(sea.gbif)[names(sea.gbif) == "max.quant.seapr"] <- paste("max.",quantile,".seapr",".gbif", sep="")
names(sea.gbif)[names(sea.gbif) == "max.quant.seapath"] <- paste("max.",quantile,".seapath",".gbif", sep="")
names(sea.gbif)[names(sea.gbif) == "max.quant.seapath.antimeridian"] <- paste("max.",quantile,".seapath.antimeridian",".gbif", sep="")
sea.gbif.1 <- sea.gbif

quantile = 95
sea.gbif <- pwseadist(points = gbif, resistance.matrix = resistance.matrix, quantile = quantile, resistance.matrix.antimeridian = resistance.matrix.antimeridian)
sea.gbif <- within(sea.gbif, rm(pw.seadist)) #don't save pw dist matrix for gbif points (but we do for genetic pts)
names(sea.gbif)[names(sea.gbif) == "max.quant.sea"] <- paste("max.",quantile,".sea",".gbif", sep="")
names(sea.gbif)[names(sea.gbif) == "max.quant.seapr"] <- paste("max.",quantile,".seapr",".gbif", sep="")
names(sea.gbif)[names(sea.gbif) == "max.quant.seapath"] <- paste("max.",quantile,".seapath",".gbif", sep="")
names(sea.gbif)[names(sea.gbif) == "max.quant.seapath.antimeridian"] <- paste("max.",quantile,".seapath.antimeridian",".gbif", sep="")
sea.gbif.2 <- sea.gbif

# genetic samples
quantile = 100
sea.genetic <- pwseadist(points = genetic, resistance.matrix = resistance.matrix, quantile = quantile, resistance.matrix.antimeridian = resistance.matrix.antimeridian)
names(sea.genetic)[names(sea.genetic) == "max.quant.sea"] <- paste("max.",quantile,".sea",".genetic", sep="")
names(sea.genetic)[names(sea.genetic) == "pw.seadist"] <- paste("pw.seadist",".genetic", sep="")
names(sea.genetic)[names(sea.genetic) == "max.quant.seapr"] <- paste("max.",quantile,".seapr",".genetic", sep="")
names(sea.genetic)[names(sea.genetic) == "max.quant.seapath"] <- paste("max.",quantile,".seapath",".genetic", sep="")
names(sea.genetic)[names(sea.genetic) == "max.quant.seapath.antimeridian"] <- paste("max.",quantile,".seapath.antimeridian",".genetic", sep="")



#*****************************************************************
#*****************************************************************
# save -------

max_and_pw_dists <- c(gcd.gbif.1, gcd.gbif.2, gcd.genetic, sea.gbif.1, sea.gbif.2, sea.genetic, gcd.landsamps)

save(max_and_pw_dists,file=paste(outdir,"/max_and_pw_dists",run_name,".Robj",sep = ""))
write.csv(gbif, paste(outdir,"/gbif-coordspostmarmap-",run_name,".csv",sep=""))
write.csv(genetic, paste(outdir,"/genetic-coordspostmarmap-",run_name,".csv",sep=""))

#save (big) pairwise dist matricies for gbif data as separate objects
pw.gbif <- list("pw.seadist.gbif" = pw.seadist.gbif,
                "pw.gcd.gbif" = pw.gcd.gbif)
save(pw.gbif,file=paste(outdir,"/GBIF_pw_dists",run_name,".Robj",sep = ""))



#*****************************************************************
#*****************************************************************
# plot data on map ------------------

load(paste(outdir,"/max_and_pw_dists",run_name,".Robj",sep = ""))
gbif <- read.csv(paste(outdir,"/gbif-coordspostmarmap-",run_name,".csv",sep=""))
genetic <- read.csv(paste(outdir,"/genetic-coordspostmarmap-",run_name,".csv",sep=""))

# plot one full world map plus other zoomed in views

bathmap.simple <- bathymap
bathmap.simple[bathmap.simple > 0] <- 5000
bathmap.simple[bathmap.simple <= 0] <- -5000

bathymap.antimeridian.simple <- bathymap.antimeridian
bathymap.antimeridian.simple[bathymap.antimeridian.simple > 0] <- 5000
bathymap.antimeridian.simple[bathymap.antimeridian.simple <= 0] <- -5000


#get zoomed bathymap (just for non antimeridian case)
maxlat <- max(gbif$y.orig, genetic$y.orig, 
              max_and_pw_dists$max.95.seapath.gbif$y, max_and_pw_dists$max.100.seapr.genetic$y,
              max_and_pw_dists$max.95.gcdpr.gbif$y, max_and_pw_dists$max.100.gcdpr.genetic$y) + 0.1
minlat <- min(gbif$y.orig, genetic$y.orig, 
              max_and_pw_dists$max.95.seapath.gbif$y, max_and_pw_dists$max.100.seapr.genetic$y,
              max_and_pw_dists$max.95.gcdpr.gbif$y, max_and_pw_dists$max.100.gcdpr.genetic$y) - 0.1
maxlon <- max(gbif$x.orig, genetic$x.orig, 
              max_and_pw_dists$max.95.seapath.gbif$x, max_and_pw_dists$max.100.seapr.genetic$x,
              max_and_pw_dists$max.95.gcdpr.gbif$x, max_and_pw_dists$max.100.gcdpr.genetic$x) + 0.1
minlon <- min(gbif$x.orig, genetic$x.orig, 
              max_and_pw_dists$max.95.seapath.gbif$x, max_and_pw_dists$max.100.seapr.genetic$x,
              max_and_pw_dists$max.95.gcdpr.gbif$x, max_and_pw_dists$max.100.gcdpr.genetic$x) - 0.1

if (maxlat > 90) {maxlat = 90}
if (minlat < -90) {minlat = -90}
if (maxlon > 180) {maxlon = 180}
if (minlon < -180) {minlon = -180}

bathymapZOOMED <- getNOAA.bathy(lat1 = minlat, lon1 = minlon, lat2 = maxlat, lon2 = maxlon, resolution = 4, keep = FALSE, antimeridian = FALSE, path = NULL)

bathymapZOOMED.simple <- bathymapZOOMED
bathymapZOOMED.simple[bathymapZOOMED.simple > 0] <- 5000
bathymapZOOMED.simple[bathymapZOOMED.simple <= 0] <- -5000




pdf(file = paste(outdir,"/distance-plots-",run_name,".pdf", sep = ""), width = 14, height = 10)

#everything - standard aka not antimeridian -----
plot1 <- ggplot() + 
  #build map
  geom_raster(data = bathmap.simple, aes(x=x, y=y, fill=z)) +
  geom_contour(data = bathymap, aes(x=x, y=y, z=z),
               breaks=0, #contour for land
               colour="black", size=0.3) +
  scale_fill_gradient2(low="#deebff", mid="#daebfb", high="gray90", midpoint = 0) +
  #add points
  geom_point(data = gbif, aes(x=x.orig, y=y.orig, color = "blue"), alpha = 0.35, size = 1) +
  geom_point(data = genetic, aes(x=x.orig, y=y.orig, color = "magenta"), alpha = 1, size = 1) +
  #add max sea dists
  geom_path(data = max_and_pw_dists$max.95.seapath.gbif, aes(x=x, y=y, color = "chartreuse2")) +
  geom_path(data = max_and_pw_dists$max.100.seapath.genetic, aes(x=x, y=y, colour = "#ff9f0f"), linetype = "dashed") +
  #add max GCD's
  geom_line(data = max_and_pw_dists$max.95.gcdpr.gbif, aes(x=x, y=y, color = "forestgreen")) +
  geom_line(data = max_and_pw_dists$max.100.gcdpr.genetic, aes(x=x, y=y, color = "red"), linetype = "dashed") +
  #manually make a legend for points and lines
  scale_color_identity(guide = "legend",
                       name = "",
                       breaks = c("blue", "magenta", 
                                  "chartreuse2", "#ff9f0f",
                                  "forestgreen", "red"),
                       labels = c("Occurrence reports (GBIF)", "Genetic samples (INSDC)",
                                  "Max 95% sea dist. GBIF","Max sea dist. INSDC",
                                  "Max 95% GCD GBIF", "Max GCD INSDC")
  ) +
  labs(title = "Full world map, normal meridian",
       subtitle = paste0("Dataset: ", run_name,
                         "\nMax 95% GCD GBIF = ", round(max_and_pw_dists$max.95.gcd.gbif,0), " km",
                         "\nMax 95% sea dist. GBIF = ", round(max_and_pw_dists$max.95.sea.gbif,0), " km",
                         "\nNOTE - drawn distance paths may not reflect actual max dist (bc of antimeridian)"),
       x = "Longitude",
       y = "Latitude") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(fill = "none") +
  coord_fixed()
print(plot1)
#

#everything - antimeridian -----
#note - sea path calced using a resistance/trans matrix from 0 - 360 are also on scale 0 - 360
gbif.anti <- gbif %>% mutate(x.orig = ifelse(x.orig<0,x.orig+360,x.orig))
genetic.anti <- genetic %>% mutate(x.orig = ifelse(x.orig<0,x.orig+360,x.orig))
max.95.gcdpr.gbif.anti <- max_and_pw_dists$max.95.gcdpr.gbif %>% mutate(x = ifelse(x<0,x+360,x))
max.100.gcdpr.genetic.anti <- max_and_pw_dists$max.100.gcdpr.genetic %>% mutate(x = ifelse(x<0,x+360,x))

plot1.5 <- ggplot() + 
  #build map
  geom_raster(data = bathymap.antimeridian.simple, aes(x=x, y=y, fill=z)) +
  geom_contour(data = bathymap.antimeridian, aes(x=x, y=y, z=z),
               breaks=0, #contour for land
               colour="black", size=0.3) +
  scale_fill_gradient2(low="#deebff", mid="#daebfb", high="gray90", midpoint = 0) +
  #add points
  geom_point(data = gbif.anti, aes(x=x.orig, y=y.orig, color = "blue"), alpha = 0.35, size = 1) +
  geom_point(data = genetic.anti, aes(x=x.orig, y=y.orig, color = "magenta"), alpha = 1, size = 1) +
  #add max sea dists
  geom_path(data = max_and_pw_dists$max.95.seapath.antimeridian.gbif, aes(x=x, y=y, color = "chartreuse2")) +
  geom_path(data = max_and_pw_dists$max.100.seapath.antimeridian.genetic, aes(x=x, y=y, colour = "#ff9f0f"), linetype = "dashed") +
  #add max GCD's
  geom_line(data = max.95.gcdpr.gbif.anti, aes(x=x, y=y, color = "forestgreen")) +
  geom_line(data = max.100.gcdpr.genetic.anti, aes(x=x, y=y, color = "red"), linetype = "dashed") +
  #manually make a legend for points and lines
  scale_color_identity(guide = "legend",
                       name = "",
                       breaks = c("blue", "magenta", 
                                  "chartreuse2", "#ff9f0f",
                                  "forestgreen", "red"),
                       labels = c("Occurrence reports (GBIF)", "Genetic samples (INSDC)",
                                  "Max 95% sea dist. GBIF","Max sea dist. INSDC",
                                  "Max 95% GCD GBIF", "Max GCD INSDC")
  ) +
  labs(title = "Full world map, ANTI meridian",
       subtitle = paste0("Dataset: ", run_name,
                         "\nMax 95% GCD GBIF = ", round(max_and_pw_dists$max.95.gcd.gbif,0), " km",
                         "\nMax 95% sea dist. GBIF = ", round(max_and_pw_dists$max.95.sea.gbif,0), " km",
                         "\nNOTE - drawn distance paths may not reflect actual max dist (bc of antimeridian)"),
       x = "Longitude",
       y = "Latitude") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(fill = "none") +
  coord_fixed()
print(plot1.5)
#

# just sea dists ZOOMED ------
plot2 <- ggplot() + 
  #build map
  geom_raster(data = bathymapZOOMED.simple, aes(x=x, y=y, fill=z)) +
  geom_contour(data = bathymapZOOMED, aes(x=x, y=y, z=z),
               breaks=0, #contour for land
               colour="black", size=0.3) +
  scale_fill_gradient2(low="#deebff", mid="#daebfb", high="gray90", midpoint = 0) +
  #add points
  geom_point(data = gbif, aes(x=x.orig, y=y.orig, color = "blue"), alpha = 0.35, size = 1) +
  geom_point(data = genetic, aes(x=x.orig, y=y.orig, color = "magenta"), alpha = 1, size = 1) +
  #add max sea dists
  geom_path(data = max_and_pw_dists$max.95.seapath.gbif, aes(x=x, y=y, color = "chartreuse2")) +
  geom_path(data = max_and_pw_dists$max.100.seapath.genetic, aes(x=x, y=y, colour = "#ff9f0f"), linetype = "dashed") +
  #manually make a legend for points and lines
  scale_color_identity(guide = "legend",
                       name = "",
                       breaks = c("blue", "magenta", 
                                  "chartreuse2", "#ff9f0f",
                                  "forestgreen", "red"),
                       labels = c("Occurrence reports (GBIF)", "Genetic samples (INSDC)",
                                  "Max 95% sea dist. GBIF","Max sea dist. INSDC",
                                  "Max 95% GCD GBIF", "Max GCD INSDC")
  ) +
  labs(title = "Zoomed map, normal meridian, SEA DISTS",
       subtitle = paste0("Dataset: ", run_name,
                         "\nMax 95% sea dist. GBIF = ", round(max_and_pw_dists$max.95.sea.gbif,0), " km",
                         "\nNOTE - drawn distance paths may not reflect actual max dist (bc of antimeridian)"),
       x = "Longitude",
       y = "Latitude") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(fill = "none") +
  coord_fixed()
print(plot2)

# just GCDs ZOOMED --------
plot3 <- ggplot() + 
  #build map
  geom_raster(data = bathymapZOOMED.simple, aes(x=x, y=y, fill=z)) +
  geom_contour(data = bathymapZOOMED, aes(x=x, y=y, z=z),
               breaks=0, #contour for land
               colour="black", size=0.3) +
  scale_fill_gradient2(low="#deebff", mid="#daebfb", high="gray90", midpoint = 0) +
  #add points
  geom_point(data = gbif, aes(x=x.orig, y=y.orig, color = "blue"), alpha = 0.35, size = 1) +
  geom_point(data = genetic, aes(x=x.orig, y=y.orig, color = "magenta"), alpha = 1, size = 1) +
  #add max GCD's
  geom_line(data = max_and_pw_dists$max.95.gcdpr.gbif, aes(x=x, y=y, color = "forestgreen")) +
  geom_line(data = max_and_pw_dists$max.100.gcdpr.genetic, aes(x=x, y=y, color = "red"), linetype = "dashed") +
  #manually make a legend for points and lines
  scale_color_identity(guide = "legend",
                       name = "",
                       breaks = c("blue", "magenta", 
                                  "chartreuse2", "#ff9f0f",
                                  "forestgreen", "red"),
                       labels = c("Occurrence reports (GBIF)", "Genetic samples (INSDC)",
                                  "Max 95% sea dist. GBIF","Max sea dist. INSDC",
                                  "Max 95% GCD GBIF", "Max GCD INSDC")
  ) +
  labs(title = "Zoomed map, normal meridian, GREAT CIRCLE DISTS",
       subtitle = paste0("Dataset: ", run_name,
                         "\nMax 95% GCD GBIF = ", round(max_and_pw_dists$max.95.gcd.gbif,0), " km",
                         "\nNOTE - drawn distance paths may not reflect actual max dist (bc of antimeridian)"),
       x = "Longitude",
       y = "Latitude") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(fill = "none") +
  coord_fixed()
print(plot3)

#just lines ZOOMED --------
plot4 <- ggplot() + 
  #build map
  geom_raster(data = bathymapZOOMED.simple, aes(x=x, y=y, fill=z)) +
  geom_contour(data = bathymapZOOMED, aes(x=x, y=y, z=z),
               breaks=0, #contour for land
               colour="black", size=0.3) +
  scale_fill_gradient2(low="#deebff", mid="#daebfb", high="gray90", midpoint = 0) +
  #add max sea dists
  geom_path(data = max_and_pw_dists$max.95.seapath.gbif, aes(x=x, y=y, color = "chartreuse2"), size = 1.3) +
  geom_path(data = max_and_pw_dists$max.100.seapath.genetic, aes(x=x, y=y, colour = "#ff9f0f"), size = 1.3, linetype = "dashed") +
  #add max GCD's
  geom_line(data = max_and_pw_dists$max.95.gcdpr.gbif, aes(x=x, y=y, color = "forestgreen"), size = 1.3) +
  geom_line(data = max_and_pw_dists$max.100.gcdpr.genetic, aes(x=x, y=y, color = "red"), size = 1.3, linetype = "dashed") +
  #manually make a legend for points and lines
  scale_color_identity(guide = "legend",
                       name = "",
                       breaks = c("blue", "magenta", 
                                  "chartreuse2", "#ff9f0f",
                                  "forestgreen", "red"),
                       labels = c("Occurrence reports (GBIF)", "Genetic samples (INSDC)",
                                  "Max 95% sea dist. GBIF","Max sea dist. INSDC",
                                  "Max 95% GCD GBIF", "Max GCD INSDC")
  ) +
  labs(title = "Zoomed map, normal meridian, DISTS",
       subtitle = paste0("Dataset: ", run_name,
                         "\nMax 95% GCD GBIF = ", round(max_and_pw_dists$max.95.gcd.gbif,0), " km",
                         "\nMax 95% sea dist. GBIF = ", round(max_and_pw_dists$max.95.sea.gbif,0), " km",
                         "\nNOTE - drawn distance paths may not reflect actual max dist (bc of antimeridian)"),
       x = "Longitude",
       y = "Latitude") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(fill = "none") +
  coord_fixed()
print(plot4)


if (nrow(landsamps)>0) {
  #original points vs. points moved from land to water ------
  pointinspect.gbif <- rbind(gbif %>% dplyr::select(run_acc_sra,x,y) %>% mutate(type = "movedtowater"), 
                             gbif %>% dplyr::select(run_acc_sra,x.orig,y.orig) %>% rename("x" = "x.orig", "y" = "y.orig") %>% 
                               mutate(type = "original")) %>% 
    mutate(dataset = "gbif")
  pointinspect.genetic <- rbind(genetic %>% dplyr::select(run_acc_sra,x,y) %>% mutate(type = "movedtowater"), 
                                genetic %>% dplyr::select(run_acc_sra,x.orig,y.orig) %>% rename("x" = "x.orig", "y" = "y.orig") %>% 
                                  mutate(type = "original")) %>% 
    mutate(dataset = "genetic")
  pointinspect <- rbind(pointinspect.gbif, pointinspect.genetic)
  
  plot5 <- ggplot() +
    geom_raster(data = bathymapZOOMED.simple, aes(x=x, y=y, fill=z)) +
    geom_contour(data = bathymapZOOMED, aes(x=x, y=y, z=z),
                 breaks=0, #contour for land
                 colour="black", size=0.05) +
    scale_fill_gradient2(low="white", mid="#daebfb", high="gray97", midpoint = 0) +
    geom_point(data = pointinspect, aes(x=x, y=y, colour=type, shape=dataset), alpha = 0.15, size = 1) +
    scale_color_manual(values = c("red","blue")) +
    labs(title = "Zoomed map, normal meridian, MOVING POINTS FROM LAND TO WATER",
         subtitle = paste0("Dataset: ", run_name),
         x = "Longitude",
         y = "Latitude") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    guides(fill = "none") +
    coord_fixed()
  
  print(plot5)
  
  #points AND lines --------
  plot6 <- ggplot() +
    geom_raster(data = bathymapZOOMED.simple, aes(x=x, y=y, fill=z)) +
    geom_contour(data = bathymapZOOMED, aes(x=x, y=y, z=z),
                 breaks=0, #contour for land
                 colour="black", size=0.05) +
    scale_fill_gradient2(low="white", mid="#daebfb", high="gray97", midpoint = 0) +
    geom_point(data = pointinspect, aes(x=x, y=y, colour=type, shape=dataset), alpha = 0.15, size = 1) +
    scale_color_manual(values = c("red","blue")) +
    geom_line(data = pointinspect, aes(x=x, y=y, group = run_acc_sra), size = 0.25) +
    labs(title = "Zoomed map, normal meridian, MOVING POINTS FROM LAND TO WATER - CONNECTOR LINES",
         subtitle = paste0("Dataset: ", run_name),
         x = "Longitude",
         y = "Latitude") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    guides(fill = "none") +
    coord_fixed()
  
  print(plot6)
  
  #correlation btwn original GCD and GCD after land points moved to water --------
  landsampcorr.moved <- max_and_pw_dists$pw.gcd.landsamps %>% as.data.frame() %>% mutate(samp = rownames(.)) %>% 
    dplyr::select(samp, everything()) %>% pivot_longer(., names_to = "samp2", values_to = "pw.gcd.moved", cols = 2:ncol(.)) %>% 
    mutate(type = ifelse(grepl("gbif", samp), "gbif", "genetic")) %>% arrange(samp,samp2)
  
  landsampcorr.allgenetic <- max_and_pw_dists$pw.gcd.genetic %>% as.data.frame() %>% mutate(samp = rownames(.)) %>% 
    dplyr::select(samp, everything()) %>% pivot_longer(., names_to = "samp2", values_to = "pw.gcd.og", cols = 2:ncol(.))
  
  landsampcorr.allgbif <- pw.gcd.gbif %>% as.data.frame() %>% mutate(samp = rownames(.)) %>% 
    dplyr::select(samp, everything()) %>% pivot_longer(., names_to = "samp2", values_to = "pw.gcd.og", cols = 2:ncol(.))
  
  landsampcorr.og <- rbind(landsampcorr.allgenetic, landsampcorr.allgbif) %>% arrange(samp,samp2) #this has more points than landsampcorr.moved bc it includes all points, not just those moved from land
  
  landsampcorr <- merge(landsampcorr.moved, landsampcorr.og, by = c("samp","samp2"), all.x = T) %>% 
    filter(is.na(pw.gcd.og)==F) %>%  #above gives pairwise distances btwn all samples, including gbif to genetic, which we don't want, so filter those out here (the NAs, bc we didn't calc og pw gcd betweeen gbif and genetic)
    filter(samp != samp2) %>% mutate(pw.gcd.difference = pw.gcd.og - pw.gcd.moved)
  
  rm(landsampcorr.allgbif,landsampcorr.allgenetic,landsampcorr.og,landsampcorr.moved)
  
  plot7 <- ggplot() +
    geom_point(data = landsampcorr, aes(x=pw.gcd.og, y=pw.gcd.moved, colour=type), alpha = 0.15, size = 2.5) +
    geom_abline(slope = 1, intercept = 0, size = 0.5, colour = "black") +
    scale_color_manual(values = c("orange","green")) +
    labs(title = "Correlation between pw GCD for points moved from land to water",
         subtitle = paste0("Dataset: ", run_name),
         x = "pairwise great circle distance - original point",
         y = "pairwise great circle distance - moved to water") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    guides(fill = "none") +
    coord_fixed() +
    facet_grid(~type)
  
  print(plot7)
  
  #histogram of difference between original pairwise GCD and new pairwise GCD after land points moved to water --------
  plot8 <- ggplot() +
    geom_histogram(data = landsampcorr, aes(x=pw.gcd.difference, fill = type), binwidth = 1, colour = "black", size=0.2) +
    labs(title = "Histogram of difference between pw GCD before and after moving points to water (just landsamps plotted)",
         subtitle = paste0("Dataset: ", run_name),
         x = "Difference in pairwise great circle distance before and after moving land points to water (km)",
         y = "Frequency") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    guides(fill = "none") +
    facet_wrap(~type, scales = "free")
  
  print(plot8)
  
  #histogram of distances points moved from land to water --------
  plot9 <- ggplot() +
    geom_histogram(data = landsamps, aes(x=total_dist_moved_km), binwidth = 1, colour = "black", fill = "gray70") +
    labs(title = "Histogram of total distance land points moved to be in water (GBIF and genetic points lumped)",
         subtitle = paste0("Dataset: ", run_name),
         x = "Distance moved (km)",
         y = "Frequency") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    guides(fill = "none") +
    coord_fixed()
  
  print(plot9)
  
}


dev.off()

#end

