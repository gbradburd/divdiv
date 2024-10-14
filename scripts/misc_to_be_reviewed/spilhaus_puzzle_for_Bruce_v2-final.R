#Bruce code

rm(list=ls())
gc()
library(ncdf4)

#load Spilhaus projection functions
source("scripts/spilhaus_functions.R") 

#get marine ecoregions
#original V Bruce worked with - ecor's stop at coastline
#individuals parts of ecoregions split out--aggregating them
ecor<-terra::aggregate(terra::vect("/Users/rachel/data0/meow_ecos_expl_clipped_expl.shp"),
                       by="ECOREGION",fun="modal",count=FALSE)
                       
#V we used in divdiv - where ecor's buffered and extend inland some
ecor<-terra::vect("/Users/rachel/divdiv/data/abiotic/input_and_working/Marine_Ecoregions_Of_the_World_(MEOW)-shp/Marine_Ecoregions_Of_the_World__MEOW_.shp")
ecor<-terra::project(ecor,"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
ggplot() + tidyterra::geom_spatvector(data = ecor) + theme_bw()

#get coastlines
NE_coastlines <- terra::vect(rnaturalearth::ne_coastline(scale = "medium", returnclass = 'sf'))

#gridlines
gridlines<-terra::vect(gen.gridlines(c(-180,-90),c(180,90),c(36,18),c(3000,1800)),
                       "lines")

#landmask
landmask<-download_land_mask()
lands<-terra::rast(t(ncvar_get(landmask)))
terra::ext(lands)<-c(-180,180,-90,90)
lands<-lands==1
lands<-terra::as.polygons(lands)[2]

#preprocess (split regions/lines along seam, etc.) and pull out geom info
#width controls the width of the seam (ensures split polygons end up on correct side)
#might need to be tweaked in some cases
#found that as low as 100 m worked for buffered V of data (that extends beyond coasts)! 
#Changed function default accordingly
#width of 2.5km aka width=2500 needed for V of ecor's that stop at coastline
ecor.coords<-as.data.frame(terra::geom(preproc.split(ecor,width=100)))
coast.coords<-as.data.frame(terra::geom(preproc.split(NE_coastlines,width=100)))
grid.coords<-as.data.frame(terra::geom(preproc.split(gridlines,width=100)))
lands.coords<-as.data.frame(terra::geom(preproc.split(lands,width=100)))
ggplot() + geom_polygon(data=ecor.coords %>% mutate(id = paste0(geom,"_",part)), aes(x=x, y=y, group=id), fill = "gray80", colour = "black", alpha = 1) + coord_fixed()

#convert lat/longs to spilhaus
ecor.coords[,c("lon","lat")]<-from_lonlat_to_spilhaus_xy(ecor.coords$x,
                                                         ecor.coords$y)
coast.coords[,c("lon","lat")]<-from_lonlat_to_spilhaus_xy(coast.coords$x,
                                                          coast.coords$y)
grid.coords[,c("lon","lat")]<-from_lonlat_to_spilhaus_xy(grid.coords$x,
                                                         grid.coords$y)
lands.coords[,c("lon","lat")]<-from_lonlat_to_spilhaus_xy(lands.coords$x,
                                                          lands.coords$y)
ggplot() + geom_polygon(data=ecor.coords %>% mutate(id = paste0(geom,"_",part)), aes(x=lon, y=lat, group=id), fill = "gray80", colour = "black", alpha = 1) + coord_fixed()


#convert back to SpatVector in Spilhaus coordinates
ecor.prj<-terra::vect(as.matrix(ecor.coords[,c("geom","part","lon","lat","hole")]),
                      "polygons",atts=terra::values(ecor))
coast.prj<-terra::vect(as.matrix(coast.coords[,c("geom","part","lon","lat","hole")]),
                       "lines",atts=terra::values(NE_coastlines))
grid.prj<-terra::vect(as.matrix(grid.coords[,c("geom","part","lon","lat","hole")]),
                      "lines",atts=terra::values(gridlines))
lands.prj<-terra::vect(as.matrix(lands.coords[,c("geom","part","lon","lat","hole")]),
                       "polygons",atts=terra::values(lands))
ggplot() + tidyterra::geom_spatvector(data = ecor.prj) + theme_bw()


#getting rid of Pacific seams (only really necessary for polygons, not lines)
#three parameters to potentially tweak:
# - intersection.check.tol --> Controls the thickness of line used to check
#     for overlap with "seam" of lonlat projection (180 longitude line).
#     Haven't played around with too much, but 1000 seems fine.
# - snap.tol --> Controls the radius at which borders are "snapped" together
#     (because imprecision of transform makes it so that borders along the
#     seam aren't perfectly touching and/or overlapping).
#     Lower might be ideal, but 5000 was the lowest I could get it to work for
#     all cases!
# - exp.fudge.fac --> Controls how much regions are buffered (i.e., expanded)
#     following snapping--this was necessary in some cases to get things
#     working correctly.
#     Obviously not ideal since it slightly distorts the geometry, but otherwise
#     you need to use some crazy high snap.tol values that ends up making things
#     much worse!
#     100 was the lowest value I could get away with
#all units are in the Spilhaus projection coordinate units, which obviously
# isn't ideal for interpretation...
#the values below are the default values for the function
ecor.prj<-fix.pacific.seams(ecor.prj,
                            intersection.check.tol=1000,
                            snap.tol=5000,
                            exp.fudge.fac=100)
ggplot() + tidyterra::geom_spatvector(data = ecor.prj) + theme_bw()


#fix landmask corners...
lim<-11825474
wd1<-0.575*lim
wd2<-0.15*lim
corners<-list(cbind(c(-lim,-lim+wd1,-lim+wd1,-lim),
                    c(-lim,-lim,-lim+wd1,-lim+wd1)),
              cbind(c(lim-wd1,lim,lim,lim-wd1),
                    c(lim-wd1,lim-wd1,lim,lim)),
              cbind(c(-lim,-lim+wd2,-lim+wd2,-lim),
                    c(lim-wd2,lim-wd2,lim,lim)),
              cbind(c(lim-wd2,lim,lim,lim-wd2),
                    c(-lim,-lim,-lim+wd2,-lim+wd2)))
lands.prj<-terra::aggregate(terra::union(terra::buffer(lands.prj,0),
                                         terra::vect(corners,"polygons")))
#use to erase land gridlines
grid.prj<-terra::erase(grid.prj,lands.prj)

####NEW STUFF 9/28####
##A COUPLE DIFFERENT SOLUTIONS FOR DEALING WITH DISTORTION IN TOPRIGHT CORNER##

#use landmask to erase inland parts of marine ecoregions
ecor.prj.1<-terra::erase(ecor.prj,lands.prj)
ggplot() + tidyterra::geom_spatvector(data = ecor.prj.1) + theme_bw()

#use above corner-padding trick to fix corner of offending ecoregion
corners<-cbind(c(lim-wd1,lim,lim,lim-wd2),
               c(lim,lim,lim-wd1,lim-wd2))
terra::plot(ecor.prj)
polygon(corners)
prob<-which(terra::values(ecor.prj)$ECOREGION=="East China Sea")
tmp<-terra::aggregate(terra::union(ecor.prj[prob],
                                   terra::vect(corners,"polygon")))
terra::values(tmp)<-terra::values(ecor.prj)[prob,]
ecor.prj<-rbind(ecor.prj,tmp)
ecor.prj.2<-ecor.prj[-prob]
ggplot() + tidyterra::geom_spatvector(data = ecor.prj.2) + theme_bw()
ggplot() + tidyterra::geom_spatvector(data = lands.prj, fill = "gray80", colour = "black") + 
  tidyterra::geom_spatvector(data = ecor.prj.2, fill = "blue", alpha = 0.2) + theme_bw()

ggplot() + tidyterra::geom_spatvector(data = ecor) + theme_bw()
prob.poly<-subset(ecor, ecor$ECOREGION=="East China Sea")
ggplot() + tidyterra::geom_spatvector(data = ecor) + tidyterra::geom_spatvector(data = prob.poly, fill = "red") + theme_bw()

#I'm most familiar with base plot, but there is a tidyterra package for ggplot and terra!
#I would recommend it since terra has a lot of built-in functionality for dealing with multi-part polygons, holes, etc.
terra::plot(ecor.prj,
            col=hcl.colors(20,"dark2",alpha=0.5),
            border=hcl.colors(20,"set2",alpha=0.5))
terra::plot(lands.prj,add=TRUE,col="gray")
terra::plot(grid.prj,add=TRUE)


#just curious to look at how spilhaus grid looks on typical projection
wd<-0.05*lim
gridlines<-terra::vect(gen.gridlines(c(-lim+wd,-lim+wd),c(lim-wd,lim-wd),c(20,20),c(2000,2000)),
                       "lines")

#map -180/180 longitude line onto spilhaus projection
seam.lonlat<-cbind(180,seq(-90,90,length.out=1000))
seam.xy<-from_lonlat_to_spilhaus_xy(seam.lonlat[,1],seam.lonlat[,2])
#split out separate components (thankfully one half always has positive y coords, the other always negative)
seam.xy<-lapply(split(seam.xy,seam.xy[,2]>0),matrix,ncol=2)
#fatten up to make sure intersection test gets everything
seam<-terra::buffer(terra::vect(seam.xy,"lines"),100)

gridlines<-terra::erase(gridlines,seam)

grid.coords2<-as.data.frame(terra::geom(gridlines))
grid.coords2[,c("lon","lat")]<-from_spilhaus_xy_to_lonlat(grid.coords2$x,
                                                          grid.coords2$y)
grid.prj2<-terra::vect(as.matrix(grid.coords2[,c("geom","part","lon","lat","hole")]),
                       "lines",atts=terra::values(gridlines))
maps::map()
terra::plot(grid.prj2,add=TRUE)
