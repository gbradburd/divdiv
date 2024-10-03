#idea: make a spilhaus map of the world and plot marine ecoregion polygons on it

#https://spatialreference.org - to look up projection codes, like PROJ.4 codes

#one of most active threads about implementing the spilhaus projection: 
#https://github.com/OSGeo/PROJ/issues/1851

#github for only implemntation in R I am currently aware of / have found online:
#https://github.com/rtlemos/spilhaus/tree/main




#load libraries -------
library(ggplot2)
library(ncdf4)
library(dplyr)
library(tidyr)

rm(list=ls())
gc()


#load Spilhaus projection functions
source("spilhaus_functions.R") 



# make spilhaus base map -----------------
#create Splihaus map grid (table of x, y points)
spilhaus_df = make_spilhaus_xy_gridpoints(spilhaus_res=1500)
x_lo = min(spilhaus_df$x); x_hi = max(spilhaus_df$x)
y_lo = min(spilhaus_df$y); y_hi = max(spilhaus_df$y)

#convert grid to lat/lon
lonlat = from_spilhaus_xy_to_lonlat(spilhaus_df$x, spilhaus_df$y)

#download marine data to locate where ocean is 
#create mask aka pull ocean shape out of square grid
da = download_land_mask()
spilhaus_df$z = extract_mask(da, lonlat)
spilhaus_df$l = is.na(spilhaus_df$z)

# pretty up projected layers
#not exactly sure what this function does - basically makes things pretty for plotting/does some cleanup
pretty_spilhaus_df = pretify_spilhaus_df(spilhaus_df)




# get marine ecoregions ------------------
ecor <- terra::vect("Marine_Ecoregions_Of_the_World__MEOW_.shp")
#reproject into norm lat/long and pull out coords
PROJ.latlon <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
ecor <- terra::project(ecor, PROJ.latlon)
#pull out lat/longs
ecor <- terra::geom(ecor) %>% as.data.frame() %>% mutate(grp = paste0(geom,"_",part))
#do they look right in current proj? (yes)
ggplot() + geom_polygon(data = ecor, aes(x=x, y=y, group=grp), fill = "transparent", colour = "orange", alpha = 0.2)
#convert lat/longs to spilhaus and add grouping variable back on
ecor.prj = from_lonlat_to_spilhaus_xy(ecor$x, ecor$y)
ecor.prj <- cbind(ecor.prj,ecor)
head(ecor.prj)

#hmmm, not great looking
ggplot() + geom_path(data = ecor.prj, aes(x=spilhaus_x, y=spilhaus_y, group=grp), colour = "orange", alpha = 1)




#find the polygons that are problems using some math
ecor.prj %>% group_by(grp) %>% mutate(range.x = max(spilhaus_x)-min(spilhaus_x)) %>% 
  mutate(range.y = max(spilhaus_y)-min(spilhaus_y)) %>% 
  dplyr::select(geom,grp,range.x,range.y) %>% 
  distinct() %>%
  pivot_longer(., names_to = "dim", values_to = "range", cols = range.x:range.y) %>% 
  ggplot() + geom_histogram(aes(x=range)) + facet_wrap(.~dim, ncol = 1)

#pull out the prob polygons
probs <- ecor.prj %>% group_by(grp) %>% mutate(range.x = max(spilhaus_x)-min(spilhaus_x)) %>% 
  mutate(range.y = max(spilhaus_y)-min(spilhaus_y)) %>% 
  dplyr::select(geom,grp,range.x,range.y) %>% 
  distinct() %>%
  pivot_longer(., names_to = "dim", values_to = "range", cols = range.x:range.y) %>% 
  filter(range > 5e6)

#viz where prob polygons are
probs.plot <- ecor.prj %>% filter(grp %in% probs$grp)

ggplot() + 
  geom_path(data = ecor.prj, aes(x=spilhaus_x, y=spilhaus_y, group=grp), colour = "orange", alpha = 1) +
  geom_path(data = probs.plot, aes(x=spilhaus_x, y=spilhaus_y, group=grp), colour = "red", alpha = 1)

ggplot() + 
  geom_polygon(data = ecor, aes(x=x, y=y, group=grp), fill = "transparent", colour = "orange", alpha = 0.2) +
  geom_polygon(data = ecor %>% filter(grp %in% probs$grp), aes(x=x, y=y, group=grp), fill = "red", colour = "red", alpha = 0.2)



#look at just one of probs
head(probs)
prob = "34_1"
prob = "37_1"
prob = "42_1"

#what it looks like
ggplot() + 
  geom_path(data = ecor.prj %>% filter(!grp %in% probs$grp), aes(x=spilhaus_x, y=spilhaus_y, group=grp), colour = "orange", alpha = 1) +
  geom_path(data = ecor.prj %>% filter(grp == prob), aes(x=spilhaus_x, y=spilhaus_y, group=grp), colour = "red", alpha = 1)
#what it should look like
ggplot() + 
  geom_polygon(data = ecor, aes(x=x, y=y, group=grp), fill = "transparent", colour = "orange", alpha = 0.2) +
  geom_polygon(data = ecor %>% filter(grp == prob), aes(x=x, y=y, group=grp), fill = "red", colour = "red", alpha = 0.2)
#single out the specific prob points in polygon
ggplot() + 
  geom_path(data = ecor.prj %>% filter(!grp %in% probs$grp), aes(x=spilhaus_x, y=spilhaus_y, group=grp), colour = "orange", alpha = 1) +
  geom_path(data = ecor.prj %>% filter(grp == prob), aes(x=spilhaus_x, y=spilhaus_y, group=grp), colour = "red", alpha = 1) +
  geom_point(data = ecor.prj %>% filter(grp == prob) %>% filter(spilhaus_x < -1e7), aes(x=spilhaus_x, y=spilhaus_y, group=grp), colour = "purple")
#try to prettify / does this solve anything ??
test <- ecor.prj %>% filter(grp == prob) %>% filter(spilhaus_x < -1e7) %>%
  rename("lon"="x", "lat" = "y", "x" = "spilhaus_x", "y" = "spilhaus_y") %>% 
  mutate(l = FALSE) %>%
  mutate(z = grp) %>% 
  dplyr::select(x, y, z, l)
test.pretty = pretify_spilhaus_df(test)
#points do move, but not to the right place
ggplot() + 
  geom_path(data = ecor.prj %>% filter(!grp %in% probs$grp), aes(x=spilhaus_x, y=spilhaus_y, group=grp), colour = "orange", alpha = 1) +
  geom_path(data = ecor.prj %>% filter(grp == prob), aes(x=spilhaus_x, y=spilhaus_y, group=grp), colour = "red", alpha = 1) +
  geom_point(data = ecor.prj %>% filter(grp == prob) %>% filter(spilhaus_x < -1e7), aes(x=spilhaus_x, y=spilhaus_y, group=grp), colour = "purple") +
  geom_point(data = test.pretty, aes(x=x, y=y), colour = "blue")
#


