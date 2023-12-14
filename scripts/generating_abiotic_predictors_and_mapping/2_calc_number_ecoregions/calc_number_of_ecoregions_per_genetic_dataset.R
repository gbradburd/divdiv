#idea: calculate how many marine ecoregions are in range of each species 
#aka for all GBIF points for a species, how many total ecoregions total do they fall into


#load libraries
library(rgdal, warn.conflicts = FALSE, quietly = TRUE)
library(sp, warn.conflicts = FALSE, quietly = TRUE)
library(raster, warn.conflicts = FALSE, quietly = TRUE)
library(stringr, warn.conflicts = FALSE, quietly = TRUE)
library(dplyr, warn.conflicts = FALSE, quietly = TRUE)
library(tidyr, warn.conflicts = FALSE, quietly = TRUE)



#get variables for running on hpcc
args <- commandArgs(trailingOnly = TRUE)
indir = args[1]
outdir = args[1]
species = args[2]
run_name = args[3]
print(paste0("indir is ",indir))

#get marine ecoregions
#as shape file
# ecor.prj.meow <- rgdal::readOGR(dsn= "/Users/rachel/Desktop/DivDiv/mapping_GBIF_and_genetic_sample_points/Marine_Ecoregions_Of_the_World_(MEOW)-shp/", 
#                    layer="Marine_Ecoregions_Of_the_World__MEOW_",
#                    verbose=T)
ecor.prj.meow <- rgdal::readOGR(dsn= paste0(indir,"/Marine_Ecoregions_Of_the_World_(MEOW)-shp/"), 
                                layer="Marine_Ecoregions_Of_the_World__MEOW_",
                                verbose=T)
#get proj.4 code for ecoregions shape files and save it so we can project rest of data into same coord system
raster::crs(ecor.prj.meow)
#PROJ <- "+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +wktext +no_defs"

#get countries of world shape obj
# countries <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sp")
# #save for hpcc
# save(countries, file = "/Users/rachel/Desktop/DivDiv/mapping_GBIF_and_genetic_sample_points/countries.outlines.Robj")
# load("/Users/rachel/Desktop/DivDiv/mapping_GBIF_and_genetic_sample_points/countries.outlines.Robj")
load(paste0(indir,"/countries.outlines.Robj"))
countries.prj <- sp::spTransform(countries, crs(ecor.prj.meow))

#get Costello et al. 2017 marine realms
ecor.prj.costello <- rgdal::readOGR(dsn= paste0(indir,"/MarineRealmsShapeFile/"), 
                                    layer="MarineRealms",
                                    verbose=T)
ecor.prj.costello <- sp::spTransform(ecor.prj.costello, crs(ecor.prj.meow))


#get genetic points
print(paste0("starting to process dataset: ", run_name))
points <- read.delim(paste0(indir,"/lat_long_table-",run_name,".txt")) %>% dplyr::rename("decimalLatitude"="lat", "decimalLongitude"="long")
#make sure there are no points left with latitude = 90 (which break projection stuff aka spTransform())
#should be all gone already anyways, just a double check
points <- points %>% filter(decimalLatitude < 90)
#project to same coord sys as ecoregions
sp::coordinates(points) <- ~decimalLongitude+decimalLatitude
raster::crs(points)
sp::proj4string(points) <- sp::CRS("+proj=longlat +datum=WGS84")
raster::crs(points)
points.prj <- sp::spTransform(points, raster::crs(ecor.prj.meow))
#ggplot can't handle SpatialPointsDataFrame so save a normal df for visualization
points.prj.df <- as.data.frame(points.prj) %>% rename("X.prj" = "decimalLongitude", "Y.prj" = "decimalLatitude")


#check that genetic points, ecoregions (and countries) are all in same projection now
sp::identicalCRS(points.prj, ecor.prj.meow)
sp::identicalCRS(points.prj, countries.prj)
sp::identicalCRS(ecor.prj.meow, ecor.prj.costello)

#get the ecoregions that each point falls within - MEOW ----------
out.meow <- sp::over(points.prj, ecor.prj.meow, returnList = TRUE)

#wrangle list into a df
out.df.meow <- data.frame("FID"=NA, "ECO_CODE"=NA, "ECOREGION"=NA, "PROV_CODE"=NA, "PROVINCE"=NA, 
                          "RLM_CODE"=NA, "REALM"=NA, "ALT_CODE"=NA, "ECO_CODE_X"=NA, "Lat_Zone"=NA, 
                          "SHAPE_Leng"=NA, "SHAPE_Area"=NA, "point_i"=NA)
for (i in 1:length(out.meow)) {
  apoint <- out.meow[i] %>% as.data.frame() %>% mutate(i = i)
  colnames(apoint) <- c("FID", "ECO_CODE", "ECOREGION", "PROV_CODE", "PROVINCE", 
                        "RLM_CODE", "REALM", "ALT_CODE", "ECO_CODE_X", "Lat_Zone", 
                        "SHAPE_Leng", "SHAPE_Area", "point_i")
  out.df.meow <- rbind(apoint,out.df.meow)
}
out.df.meow <- out.df.meow %>% filter(is.na(i)==F)


#get the ecoregions that each point falls within - Costello ----------
out.costello <- sp::over(points.prj, ecor.prj.costello, returnList = TRUE)

#wrangle list into a df
out.df.costello <- data.frame("Realm"=NA, "point_i"=NA)
for (i in 1:length(out.costello)) {
  apoint <- out.costello[i] %>% as.data.frame() %>% mutate(point_i = i)
  out.df.costello <- rbind(apoint,out.df.costello)
}
out.df.costello <- out.df.costello %>% filter(is.na(i)==F) %>% rename("Realm_Costello" = "Realm")

out.df <- merge(out.df.meow, out.df.costello, by = "point_i", all.x = T, all.y = T) %>% mutate(link = run_name)

write.csv(out.df, paste0(outdir,"/number_of_ecoregions_.genetic.",run_name,".csv"), row.names = F)

#end


