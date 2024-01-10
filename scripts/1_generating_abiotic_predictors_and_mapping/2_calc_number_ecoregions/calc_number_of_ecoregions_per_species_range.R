#idea: calculate how many marine ecoregions are in range of each species 
#aka for all GBIF points for a species, how many total ecoregions total do they fall into


#load libraries
library(rgdal, warn.conflicts = FALSE, quietly = TRUE)
library(sp, warn.conflicts = FALSE, quietly = TRUE)
library(raster, warn.conflicts = FALSE, quietly = TRUE)
library(stringr, warn.conflicts = FALSE, quietly = TRUE)
library(dplyr, warn.conflicts = FALSE, quietly = TRUE)
library(tidyr, warn.conflicts = FALSE, quietly = TRUE)


#rm(list = ls())
#gc()

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


  #get GBIF points
  print(paste0("starting to process species: ", species))
  points <- read.csv(paste0(indir,"/GBIFLocations_",species,".csv"))
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
  
  
  #check that GBIF points, ecoregions (and countries) are all in same projection now
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

  out.df <- merge(out.df.meow, out.df.costello, by = "point_i", all.x = T, all.y = T) %>% mutate(species = species)
  
speciesnospace = gsub(" ","-",species)
write.csv(out.df, paste0(outdir,"/number_of_ecoregions_.gbif.",speciesnospace,".csv"), row.names = F)

#end




# 
# # 
# # plot MEOWs
# ggplot() +
# 
#   geom_polygon(data = countries.prj,
#                aes(long, lat, group = group),
#                colour = "gray50", fill = "white", size = 0.25) +
#   geom_polygon(data = ecor.prj,
#                aes(long,lat, group = group),
#                size = 0.5, color = "transparent", fill = "cyan1", alpha = 0.2) +
#   geom_point(data = points.prj.df, aes(X.prj, Y.prj), colour = "red", size = 0.74) +
#   #geom_text(data = points.prj.df, aes(X.prj, Y.prj, label = rownames(points.prj.df)), size = 2) +
#   geom_polygon(data = ecor.prj,
#                aes(long,lat, group = group),
#                size = 0.5, color = "black", fill = "transparent", alpha = 0.2) +
#   coord_sf() +
#   lims(y = c(-15000000,23000000))
# 
# 
# # plot Costello marine realms
# ggplot() +
#   
#   geom_polygon(data = countries.prj,
#                aes(long, lat, group = group),
#                colour = "gray50", fill = "white", size = 0.25) +
#   geom_polygon(data = ecor.prj.costello,
#                aes(long,lat, group = group),
#                size = 0.5, color = "transparent", fill = "cyan1", alpha = 0.2) +
#   geom_point(data = points.prj.df, aes(X.prj, Y.prj), colour = "red", size = 0.74) +
#   #geom_text(data = points.prj.df, aes(X.prj, Y.prj, label = rownames(points.prj.df)), size = 2) +
#   geom_polygon(data = ecor.prj.costello,
#                aes(long,lat, group = group),
#                size = 0.5, color = "black", fill = "transparent", alpha = 0.2) +
#   coord_sf() +
#   lims(y = c(-15000000,23000000))
# #
# 
# #plot MEOW ecoregions and Costello relams
# ggplot() +
#   geom_polygon(data = countries.prj,
#                aes(long, lat, group = group),
#                colour = "gray50", fill = "white", size = 0.25) +
#   geom_polygon(data = ecor.prj.costello,
#                aes(long,lat, group = group),
#                size = 0.5, color = "black", fill = "cyan1", alpha = 0.2) +
#   geom_polygon(data = ecor.prj.meow,
#                aes(long,lat, group = group),
#                size = 0.5, color = "red", fill = "red", alpha = 0.1) +
#   coord_sf() +
#   lims(y = c(-24000000,23000000))
# #
# 
# #plot MEOW realms and Costello relams
# ecor.prj.meow <- rgdal::readOGR(dsn= "/Users/rachel/Desktop/DivDiv/mapping_GBIF_and_genetic_sample_points/Marine_Ecoregions_Of_the_World_(MEOW)-shp/",
#                                 layer="Marine_Ecoregions_Of_the_World__MEOW_",
#                                 verbose=T)
# #create df of shape obj with data slots, so we can color by them
# ecor.prj.meow@data$id <- rownames(ecor.prj.meow@data)
# fort <- ggplot2::fortify(ecor.prj.meow, region = "id")
# df <- merge(fort, ecor.prj.meow@data,
#             by = "id") %>% mutate(RLM_CODE = as.factor(RLM_CODE)) %>% mutate(PROV_CODE = as.factor(PROV_CODE))
# 
# ggplot() +
#   geom_polygon(data = countries.prj,
#                aes(long, lat, group = group),
#                colour = "gray50", fill = "white", size = 0.25) +
#   geom_polygon(data = ecor.prj.costello,
#                aes(long,lat, group = group),
#                size = 0.5, color = "transparent", fill = "cyan1", alpha = 0.1) +
#   geom_polygon(data = df,
#                aes(long,lat, group = group, fill = RLM_CODE),
#                size = 0.5, alpha = 0.6, colour = "transparent") +
#   geom_polygon(data = ecor.prj.costello,
#                aes(long,lat, group = group),
#                size = 0.5, color = "black", fill = "transparent") +
#   coord_sf() +
#   lims(y = c(-24000000,23000000)) +
#   theme(panel.grid = element_blank())
# 
# ggplot() +
#   geom_polygon(data = countries.prj,
#                aes(long, lat, group = group),
#                colour = "gray50", fill = "white", size = 0.25) +
#   geom_polygon(data = ecor.prj.costello,
#                aes(long,lat, group = group),
#                size = 0.5, color = "transparent", fill = "cyan1", alpha = 0.1) +
#   geom_polygon(data = df,
#                aes(long,lat, group = group, fill = PROV_CODE),
#                size = 0.5, alpha = 0.6, colour = "transparent") +
#   geom_polygon(data = ecor.prj.costello,
#                aes(long,lat, group = group),
#                size = 0.5, color = "black", fill = "transparent") +
#   coord_sf() +
#   lims(y = c(-24000000,23000000)) +
#   theme(panel.grid = element_blank())
# #



#graveyard --------

#prj.coord <- rgdal::project(cbind(points$decimalLongitude, points$decimalLatitude), proj = PROJ)
#points.prj <- cbind(prj.coord, points)
#names(points.prj)[1:2] <- c("X.prj","Y.prj")


#as simple features
#ecor.simple <- sf::read_sf("/Users/rachel/Desktop/DivDiv/mapping_GBIF_and_genetic_sample_points/Marine_Ecoregions_Of_the_World_(MEOW)-shp/Marine_Ecoregions_Of_the_World__MEOW_.shp")


# #no contries
# ggplot() + 
#   geom_sf(data = ecor.simple, size = 0.5, color = "black", fill = "cyan1", alpha = 0.2) +
#   geom_point(data = sites.df.prj, aes(x = X.prj, y = Y.prj), colour = "red", size = 1.5) +
#   coord_sf()
# 
# 
# #using simple object for ecoregions
# ggplot() +
#   
#   geom_polygon(data = countries.prj, 
#                aes(long,lat, group = group), 
#                colour = "gray50", fill = "white", size = 0.25) +
#   geom_sf(data = ecor.simple, size = 0.5, color = "transparent", fill = "cyan1", alpha = 0.2) +
#   geom_point(data = sites.df.prj, aes(x = X.prj, y = Y.prj), colour = "red", size = 1.5) +
#   geom_sf(data = ecor.simple, size = 0.5, color = "black", fill = "transparent") +
#   coord_sf() +
#   lims(y = c(-14000000,22000000))



# #old code for converting shape obj to df (for plotting)
# df <- data.frame("V1"=NA,"V2"=NA,"i.sub"=NA,"i"=NA)
# for (i in 1:length(ecor.prj.meow)){
#   print(paste0("i is ",i))
#   out <- sapply(
#     ecor.prj.meow@polygons[[i]]@Polygons,
#     slot,
#     "coords")
#   if (is.list(out)==FALSE) {
#     out <- cbind(out[1:(length(out)/2)],out[((length(out)/2)+1):length(out)]) %>% as.data.frame() %>% mutate(i.sub = 1) %>% mutate(i=i)
#     df <- rbind(df,out)
#   } else {
#     for (i.sub in 1:length(out)) {
#       out.sub <- out[[i.sub]] %>% as.data.frame()
#       out.sub <- rbind(out.sub,out.sub[1,]) %>% mutate(i.sub = i.sub) %>% mutate(i = i)
#       df <- rbind(df,out.sub)
#     }
#   }
# }
# df <- df %>% filter(is.na(i)==FALSE) %>% mutate(group = as.factor(paste0(i,i.sub,sep = ".")))
# labels <- ecor.prj.meow@data
# df <- merge(df, labels, by.x = c("i"), by.y = c("FID"), all.x = T) %>%
#   mutate(RLM_CODE = as.factor(RLM_CODE)) %>%
#   mutate(PROV_CODE = as.factor(PROV_CODE))

