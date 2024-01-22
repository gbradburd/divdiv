

#function - calculate pairwise great circle distances
# points is a df with x, y, and run_acc_sra of the points you want to calculate GCDs between
# note x needs to be the first column and y needs to be the second column
# quantile is a number that says what quantile to chop pairwise dist distribution at
# e.g. quantile = 95 means chop off 2.5% in each tail and then get max dist
# quantile = 100 means take absolute max dist (no truncating of pairwise dist. distribution)

pwGCD <- function(points, quantile) {

  #calc quantile values
  lower = ((100 - quantile)/2)/100
  upper = (100 - (lower*100))/100
  
  #calc GCD
  #quantile(obj_of_class_dist) gives same result as long list of pairwise dists with self-self 0's excluded
  gcd <- fields::rdist.earth(points,miles=FALSE)
  rownames(gcd) <- points$run_acc_sra
  colnames(gcd) <- points$run_acc_sra 
  
  #find the max for X quantile
  max.quant.gcd = stats::quantile(gcd, c(lower, upper)) %>% .[2]
  
  #select the first pair that represent the X% quantile max dist so we can plot it on map
  max.quant.gcdpr <- gcd %>% as.data.frame() %>% mutate(samp = row.names(.)) %>% 
    dplyr::select(samp,everything()) %>% 
    pivot_longer(., 2:ncol(.), names_to = "samp2", values_to = "dist_km") %>% 
    mutate(distfrommax = abs(dist_km - max.quant.gcd)) %>% 
    filter(distfrommax == min(distfrommax)) %>% 
    .[1,]
  samp1 <- points %>% filter(run_acc_sra == max.quant.gcdpr$samp)
  samp2 <- points %>% filter(run_acc_sra == max.quant.gcdpr$samp2)
  max.quant.gcdpr  <- rbind(samp1,samp2)
  
  #save
  gcd_dists <- list("max.quant.gcd" = max.quant.gcd,
                    "pw.gcd" = gcd,
                    "max.quant.gcdpr" = max.quant.gcdpr)
  
  return(gcd_dists)
  
}



#function - move points on land to water
# landsamps is a df with x, y, and run_acc_sra of the points you want to 
# move to water (that already have been checked to have a depth of > 0)
# in points: x is decimal longitude, y is decimal latitude, and run_acc_sra is a unique sample name for the point

# bathymap which is class bathy and output from marmap's getNOAA.bathy() that will be used in downstream analyses 
#(warning that isobath lines change depending on resolution and extent downloaded)

# bathymap.antimeridian which is class bathy and output from marmap's getNOAA.bathy() that will be used in downstream analyses 
# note extent and resolution should be same between bathymap and bathymap.antimeridian, just use.antimerdian = TRUE (instead of FALSE like for bathymap)

# starting.search.radius which is a numerical value that defines search radius (in decimal degrees) to start with to find nearest water cell if 
# dist2isobath doesn't return a water cell to start

# bathyres which is a numerical value in arcminutes that matches resolution passed to getNOAA.bathy() for bathy map objects being used

# outputs landsamps df which has run_acc_sra aka sample name, x which is new decimal longitude in water, 
# y which is new decimal latitude in water, depth at point, 
# intial.end.lon and initial.end.lat which are output of marmap's dist2isobath() that should be in water (decimal degrees),
# end.lon and end.lat which are same as intial.end.lon and initial.end.lat when dist2isobath() returns a point in water
# and wiggled points (aka output of dist2isobath() moved to nearest negative grid cell when output of dist2isobath() did not end up in water),
# total_dist_moved_km which is GCD (in km) between original sample point and final end point in water,
# and wiggle_dist_moved_km which is GCD (in km) between output of dist2isobath() and final point output after adding random normal deviate
# until reached N tries or point landed in water

move_land_points_to_water <- function(landsamps, bathymap, bathymap.antimeridian, starting.search.radius, bathyres, resistance.matrix) {

  #find nearest point on coast for each landsamp
  watersamps <- data.frame("run_acc_sra" = NA, "initial.end.lon" = NA, "initial.end.lat" = NA, "end.lon" = NA, "end.lat" = NA)

  for (i in 1:nrow(landsamps)) {  
  
    #duplicate row (bc function breaks with only one point)
    landpoint <- rbind(landsamps[i,c("x","y")],landsamps[i,c("x","y")])
    samp = landsamps[i,c("run_acc_sra")]
  
    #calc great circle distance from this point to isobath X 
    #return distance in meters
    out <- try(dist2isobath(mat=bathymap, x=landpoint, isobath=0), silent = T)
    
    if ("try-error" %in% class(out)) {
      #if dist2isobath() breaks, temporarily assign 0's to all end point lat/lons, then filter them out outside of function
      print(paste0("Dropping this point bc it breaks marmap::dist2isobath - ", samp))
      waterpoint <- data.frame(run_acc_sra = samp, initial.end.lon = 0, initial.end.lat = 0, end.lon = 0, end.lat = 0)
    } else {
      
      #move the end.point from dist2isobath() to center of nearest grid cell where depth is < 0 if it's not already
      end.point <- out[1,4:5]
      end.point.initial <- end.point %>% rename("initial.end.lon" = "end.lon", "initial.end.lat" = "end.lat")
      ocean_test_point <- data.frame(x = 131.921260, y = -49.299319)
      while(get.depth(mat = bathymap, x = end.point, locator = FALSE)$depth > 0 | 
         get.depth(mat = bathymap.antimeridian, x = cbind(ifelse(end.point$end.lon < 0,end.point$end.lon + 360 , end.point$end.lon),end.point$end.lat), locator = FALSE)$depth > 0 |
         lc.dist(trans = resistance.matrix, loc = rbind(end.point %>% dplyr::rename("x"="end.lon", "y"="end.lat"),ocean_test_point), res ="dist") %>% as.matrix() %>% .[1,2] > 50000
      ){
        
        print(paste0("normal depth greater than 0: ", get.depth(mat = bathymap, x = end.point, locator = FALSE)$depth > 0))
        print(paste0("anti depth greater than 0: ", get.depth(mat = bathymap.antimeridian, x = cbind(ifelse(end.point$end.lon < 0,end.point$end.lon + 360 , end.point$end.lon),end.point$end.lat), locator = FALSE)$depth > 0))
        print(paste0("dist to atlantis bonkers: ", lc.dist(trans = resistance.matrix, loc = rbind(end.point %>% dplyr::rename("x"="end.lon", "y"="end.lat"),ocean_test_point), res ="dist") %>% as.matrix() %>% .[1,2] > 50000))
        print(paste0("end.point is: ", end.point))
        
        #Set search radius (in decimal degrees)
        starting.search.radius.using = starting.search.radius
        #set resolution of bathymap downloaded (arc minutes)
        bathyres = bathyres
        
        #Assign current point
        x <- end.point[,1]
        y <- end.point[,2]
        bathydata.alreadytried <- data.frame("end.lon"=NA, "end.lat"=NA, "z"=NA, "link"=NA)
        for (j in 1:20){ #exit after 20 loops aka rounds of expanding search radius
          #Subset points within search radius from giant dataset
          print(paste0("for testing j is ",j," and start radius is ", starting.search.radius.using))
          bathymap.sub <- bathymap[(as.numeric(rownames(bathymap))<(x+starting.search.radius.using))&(as.numeric(rownames(bathymap))>(x-starting.search.radius.using)),(as.numeric(colnames(bathymap))>(y-starting.search.radius.using))&(as.numeric(colnames(bathymap))<(y+starting.search.radius.using))]
          #Turn subsetted bathy object into dataframe
          bathydata <- expand.grid("end.lon"=as.numeric(rownames(bathymap.sub)),"end.lat"=as.numeric(colnames(bathymap.sub)))
          bathydata$z <- as.vector(bathymap.sub)
          #this is where we pull out already tested cells if we want to add that code
          if (j == 1) {
            bathydata.saved.full <- bathydata
          } else {
            bathydata <- bathydata %>% mutate(link = paste0(end.lon,"_",end.lat)) %>% filter(!link %in% bathydata.alreadytried$link) %>% dplyr::select(-link)
            bathydata.saved.full <- bathydata
          }
          
          #Remove object to save memory
          rm("bathymap.sub")
          #Find only negative z
          bathydata <- bathydata[bathydata$z<0,]
          
          #if we find at least some negative cells
          if(nrow(bathydata)>0){
            #Calculate pairwise distance to each negative grid cell (GCD bc/if longlat = TRUE)
            dist.vector <- sp::spDistsN1(as.matrix(bathydata[,c(1,2)]), c(x,y), longlat = TRUE)
            bathydata <- cbind(bathydata,dist.vector) %>% arrange(dist.vector)
            
            #Pick the nearest negative grid cell that doesn't give a bonkers distance to ocean check point
            #and has negative depth for antimeridian too
            #keeps running until either antimerid is negative AND dist to ocean check isn't bonkers, or when end of df of potential pts is reached
            nearestcellrow = 1
            end.point <- bathydata[nearestcellrow,c(1,2)]
            while(nearestcellrow < (nrow(bathydata)) & 
                  (lc.dist(trans = resistance.matrix, loc = rbind(end.point %>% dplyr::rename("x"="end.lon", "y"="end.lat"),ocean_test_point), res ="dist") %>% as.matrix() %>% .[1,2] > 50000 | 
                  get.depth(mat = bathymap.antimeridian, x = cbind(ifelse(end.point$end.lon < 0,end.point$end.lon + 360 , end.point$end.lon),end.point$end.lat), locator = FALSE)$depth > 0)) {
              print(paste0("for testing nearestcellrow is ", nearestcellrow))
              nearestcellrow = nearestcellrow+1
              end.point <- bathydata[nearestcellrow,c(1,2)]
            }
            #check if we were able to find an end.point that works from returned possible list aka bathydata
            if(lc.dist(trans = resistance.matrix, loc = rbind(end.point %>% dplyr::rename("x"="end.lon", "y"="end.lat"),ocean_test_point), res ="dist") %>% as.matrix() %>% .[1,2] > 50000) {
              arewedone=FALSE
            } else {
              arewedone=TRUE
            }
            print(paste0("are we done is: ", arewedone))
          } else {
            arewedone=FALSE #didn't find any negative grid cells, try again in wider search radius
          } #done trying to pick an end.point from possible list

          if(nrow(bathydata) > 0 & arewedone==TRUE){
            starting.search.radius.using <- starting.search.radius
            break
          }
          starting.search.radius.using <- starting.search.radius.using+1 #Increase search radius and search again when nothing is found
          bathydata.saved.full <- bathydata.saved.full %>% mutate(link = paste0(end.lon,"_",end.lat))
          bathydata.alreadytried <- rbind(bathydata.alreadytried,bathydata.saved.full) %>% filter(is.na(link)==F)
        } #done with iteratively expanding search radius
      } #done with moving end.point to suitable cell if dist.to.isobath didn't work on first try
      
      #check that point is indeed in water for both bathymap and antimeridian bathymap
      depthout1 <- get.depth(mat = bathymap, x = end.point$end.lon, y = end.point$end.lat, locator = F)$depth
      depthout2 <- get.depth(mat = bathymap.antimeridian, x = ifelse(end.point$end.lon < 0, end.point$end.lon + 360, end.point$end.lon), y = end.point$end.lat, locator = F)$depth
      if (depthout1>0 | depthout2>0){
        print("ERROR! - non negative final depths for sample ", samp)
        print(paste0("Final end.point depth for bathymap is ", depthout1, " and bathymap.antimeridian is ",depthout2))
      }
      #save new end.point that's in water
      waterpoint <- cbind(data.frame(run_acc_sra = samp), end.point.initial, end.point)  
    } #done with processing one point for which we got an isobath dist back for
  
    watersamps <- rbind(watersamps, waterpoint) 
    
    print(paste0("done with landpoint ", i, " of ", nrow(landsamps)))
    
  } #totally done with one landsamp point
  watersamps <- watersamps %>% filter(is.na(run_acc_sra)==F)

  #merge back onto original lat/lons for landsamps
  landsamps <- merge(landsamps, watersamps, by = "run_acc_sra", all.x = T) %>% arrange(run_acc_sra)
  #calc distance points were moved
  #total from start to finish
  distmoved.out1 <- sp::spDists(x = as.matrix(landsamps[,c("x","y")]), 
                              y = as.matrix(landsamps[,c("end.lon","end.lat")]), 
                              longlat = TRUE, segments = FALSE, diagonal = TRUE) %>% 
    as.data.frame() %>% rename("total_dist_moved_km" = ".") %>% 
    cbind(landsamps$run_acc_sra, .) %>% rename("run_acc_sra" = "landsamps$run_acc_sra")
  #from original dist2isobath() output to final with depth < 0
  distmoved.out2 <- sp::spDists(x = as.matrix(landsamps[,c("initial.end.lon","initial.end.lat")]), 
                              y = as.matrix(landsamps[,c("end.lon","end.lat")]), 
                              longlat = TRUE, segments = FALSE, diagonal = TRUE) %>% 
    as.data.frame() %>% rename("wiggle_dist_moved_km" = ".") %>% 
    cbind(landsamps$run_acc_sra, .) %>% rename("run_acc_sra" = "landsamps$run_acc_sra")
  distmoved <- merge(distmoved.out1, distmoved.out2, by = "run_acc_sra", all = T)
  landsamps <- merge(landsamps, distmoved, by = "run_acc_sra", all.x = T) %>% arrange(run_acc_sra)
  
  return(landsamps)
  
} #done with all points aka end of function



#function - calculate pairwise distances via isobath X (aka by sea, avoiding land)
# points is a df with x, y, and run_acc_sra of the points you want to calculate GCDs between
# resistance.matrix is a resistence matrix generated using marmap's trans.mat()
# quantile is a number that says what quantile to chop pairwise dist distribution at
# e.g. quantile = 95 means chop off 2.5% in each tail and then get max dist
# quantile = 100 means take absolute max dist (no truncating of pairwise dist. distribution)

pwseadist <- function(points, resistance.matrix, quantile, resistance.matrix.antimeridian=NULL) {
  
  #calc quantile values
  lower = ((100 - quantile)/2)/100
  upper = (100 - (lower*100))/100
  
    #calc distance via sea - NO ANTIMERIDIAN aka a map with Africa in middle and longs from -180 to 180
    #quantile(obj_of_class_dist) gives same result as long list of pairwise dists with self-self 0's excluded
  #print("here 1")
    seadist <- as.matrix(lc.dist(trans = resistance.matrix, loc = points %>% dplyr::select(x,y), res ="dist"))
    rownames(seadist) <- points$run_acc_sra 
    colnames(seadist) <- points$run_acc_sra
  
    if (is.null(resistance.matrix.antimeridian)==FALSE){
    
    print("I am also running sea distance calculations using a map with antimeridian=T and picking shorter distance btwn antimeridian=T and antimeridian=F")
    
    #calc distance via sea - USING ANTIMERIDIAN aka a map with AK, USA in middle and longs from 0 to 360
    #quantile(obj_of_class_dist) gives same result as long list of pairwise dists with self-self 0's excluded
      points.anti <- points %>% dplyr::select(x,y) %>% mutate(x = ifelse(x<0,x+360,x))
      seadist.anti <- as.matrix(lc.dist(trans = resistance.matrix.antimeridian, 
                                      loc = points.anti, 
                                      res ="dist"))
    rownames(seadist.anti) <- points.anti$run_acc_sra 
    colnames(seadist.anti) <- points.anti$run_acc_sra
    #print("here 2")
    #find and save the smaller distance for each pairwise combination
    #(assumes both seadist and seadist.anti are same dimensions)
    seadist.out <- matrix(NA, nrow = nrow(seadist), ncol = ncol(seadist))
    for (i in 1:ncol(seadist)){
      for (j in 1:ncol(seadist)){
 
        min <- min(seadist[i,j], seadist.anti[i,j])
        seadist.out[i,j] <- min
        
      }
    }
   # print("here 3")
    seadist <- seadist.out
    rownames(seadist) <- points$run_acc_sra 
    colnames(seadist) <- points$run_acc_sra
    
  }
  
  #find the max for X quantile
  max.quant.sea = quantile(seadist, c(lower, upper)) %>% .[2]
  #print("here 4")
  #select the first pair that represent the X% quantile max dist so we can plot it on map
  max.quant.seapr <- seadist %>% as.data.frame() %>% mutate(samp = row.names(.)) %>% 
    dplyr::select(samp,everything()) %>% 
    pivot_longer(., 2:ncol(.), names_to = "samp2", values_to = "dist_km") %>% 
    mutate(distfrommax = abs(dist_km - max.quant.sea)) %>% 
    filter(distfrommax == min(distfrommax)) %>% 
    .[1,]
  samp1 <- points %>% filter(run_acc_sra == max.quant.seapr$samp)
  samp2 <- points %>% filter(run_acc_sra == max.quant.seapr$samp2)
  max.quant.seapr  <- rbind(samp1,samp2)
  #print("here 5")
  #get traced path for the chosen max dist pair
  max.quant.seapath <- lc.dist(trans = resistance.matrix, loc = max.quant.seapr %>% dplyr::select(x,y), res ="path") %>% as.data.frame()
  max.quant.seapath.antimeridian <- NULL #so the save step doesn't break if we don't actually calc this below
  
  if (is.null(resistance.matrix.antimeridian)==FALSE){
    max.quant.seapath.antimeridian <- lc.dist(trans = resistance.matrix.antimeridian, 
                                             loc = max.quant.seapr %>% dplyr::select(x,y) %>% mutate(x = ifelse(x<0,x+360,x)), 
                                             res ="path") %>% as.data.frame()
  }
  #print("here 6")
  #save
  sea_dists <- list("max.quant.sea" = max.quant.sea,
                    "pw.seadist" = seadist,
                    "max.quant.seapr" = max.quant.seapr,
                    "max.quant.seapath" = max.quant.seapath,
                    "max.quant.seapath.antimeridian" = max.quant.seapath.antimeridian)
  
  return(sea_dists)
  
}



#graveyard ----

# THIS FUNCTION TAKES SUCCESSIVELY LARGER RANDOM JUMPS UNTIL A LAND POINT ENDS UP IN WATER 

# NOTE - we left this function (10-15-2021) in favor of version that searchs for nearest non-zero cell to 
#endpoint of disttoisobath if it doesn't already land in water

#function - move points on land to water
# landsamps is a df with x, y, and run_acc_sra of the points you want to 
# move to water (that already have been checked to have a depth of > 0)
# in points, x is decimal longitude, y is decimal latitude, and run_acc_sra is a unique sample name for the point

# bathymap which is class bathy and output from marmap's getNOAA.bathy() for the region the points fall within

#useantimeridian is "TRUE" or "FALSE", if TRUE, add 360 degrees to negative longitudes to flip them over and plot on 
#far righthand side of map (instead of lefthand side with 0 longitude line in middle of map)
#this should match whatever antimeridian in getNOAA.bathy() was set to

# outputs landsamps df which has run_acc_sra aka sample name, x which is new decimal longitude in water, 
# y which is new decimal latitude in water, depth at point, 
# intial.end.lon and initial.end.lat which are output of marmap's dist2isobath() that should be in water (decimal degrees),
# end.lon and end.lat which are same as intial.end.lon and initial.end.lat when dist2isobath() returns a point in water
# and wiggled points (aka output of dist2isobath() plus a random normal deviate) when output of dist2isobath() did not end up in water,
# total_dist_moved_km which is GCD (in km) between original sample point and final end point in water,
# and wiggle_dist_moved_km which is GCD (in km) between output of dist2isobath() and final point output after adding random normal deviate
# until reached N tries or point landed in water
# 
# move_land_points_to_water <- function(landsamps, bathymap, useantimeridian=NULL, bathymap.antimeridian=NULL) {
#   
#   #find nearest point on coast for each landsamp
#   watersamps <- data.frame(run_acc_sra = NA, initial.end.lon = NA, initial.end.lat = NA, end.lon = NA, end.lat = NA)
#   
#   for (i in 1:nrow(landsamps)) {  
#     
#     #duplicate row (bc function breaks with only one point)
#     landpoint <- rbind(landsamps[i,c("x","y")],landsamps[i,c("x","y")])
#     samp = landsamps[i,c("run_acc_sra")]
#     
#     #calc great circle distance from this point to isobath X 
#     #return distance in meters
#     out <- dist2isobath(mat=bathymap, x=landpoint, isobath=0)
#     #move the end.point from dist2isobath() around a tiny bit until the depth is < 0 if it's not already
#     end.point <- out[1,4:5]
#     end.point.initial <- end.point %>% rename("initial.end.lon" = "end.lon", "initial.end.lat" = "end.lat")
#     end.point.guess <- end.point
#     n_tries = 0
#     while(get.depth(mat = bathymap, x = end.point.guess, locator = FALSE)$depth > 0 | 
#           get.depth(mat = bathymap.antimeridian, x = cbind(ifelse(end.point.guess$end.lon < 0,end.point.guess$end.lon + 360 , end.point.guess$end.lon),end.point.guess$end.lat), locator = FALSE)$depth > 0 & 
#           n_tries<100000){
#       end.point.guess <- end.point + c(rnorm(2,0,0.001+n_tries%/%1e3))
#       n_tries <- n_tries + 1
#     }
#     
#     #depthout1 <- get.depth(mat = bathymap, x = end.point.guess$end.lon, y = end.point.guess$end.lat, locator = F)$depth
#     #depthout2 <- get.depth(mat = bathymap.antimeridian, x = ifelse(end.point.guess$end.lon < 0, end.point.guess$end.lon + 360, end.point.guess$end.lon), y = end.point.guess$end.lat, locator = F)$depth
#     #print(paste0("final depths are ", depthout1, " and ",depthout2))
#     
#     #get new point that's in water
#     waterpoint <- cbind(data.frame(run_acc_sra = samp), end.point.initial, end.point.guess)
#     watersamps <- rbind(watersamps, waterpoint) 
#     
#     print(paste0("done with landpoint ", i, " of ", nrow(landsamps)))
#     
#   }
#   watersamps <- watersamps %>% filter(is.na(run_acc_sra)==F)
#   
#   #merge back onto original lat/lons for landsamps
#   landsamps <- merge(landsamps, watersamps, by = "run_acc_sra", all.x = T) %>% arrange(run_acc_sra)
#   #calc distance points were moved
#   #total from start to finish
#   distmoved.out1 <- sp::spDists(x = as.matrix(landsamps[,c("x","y")]), 
#                                 y = as.matrix(landsamps[,c("end.lon","end.lat")]), 
#                                 longlat = TRUE, segments = FALSE, diagonal = TRUE) %>% 
#     as.data.frame() %>% rename("total_dist_moved_km" = ".") %>% 
#     cbind(landsamps$run_acc_sra, .) %>% rename("run_acc_sra" = "landsamps$run_acc_sra")
#   #from original dist2isobath() output to final with depth < 0
#   distmoved.out2 <- sp::spDists(x = as.matrix(landsamps[,c("initial.end.lon","initial.end.lat")]), 
#                                 y = as.matrix(landsamps[,c("end.lon","end.lat")]), 
#                                 longlat = TRUE, segments = FALSE, diagonal = TRUE) %>% 
#     as.data.frame() %>% rename("wiggle_dist_moved_km" = ".") %>% 
#     cbind(landsamps$run_acc_sra, .) %>% rename("run_acc_sra" = "landsamps$run_acc_sra")
#   distmoved <- merge(distmoved.out1, distmoved.out2, by = "run_acc_sra", all = T)
#   landsamps <- merge(landsamps, distmoved, by = "run_acc_sra", all.x = T) %>% arrange(run_acc_sra)
#   
#   return(landsamps)
#   
# }
# 

