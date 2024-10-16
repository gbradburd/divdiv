library(sf)

expand.map<-function(x,
                     amount=0.05,
                     exp.fudge.fac=10){
  if(amount>1|amount<0){
    stop("amount must be between 0 and 1")
  }
  geom.type<-terra::geomtype(x)
  lim<-11825474
  tmp<-as.data.frame(terra::geom(x))
  rotate<-function(df,ang,org=c(0,0)){
    df[,c("lon","lat")]<-list(cos(ang)*(df$x-org[1])-sin(ang)*(df$y-org[2])+org[1],
                              sin(ang)*(df$x-org[1])+cos(ang)*(df$y-org[2])+org[2])
    df
  }
  tmp<-matrix(list(rotate(tmp,pi,c(-lim,lim)),
                   rotate(tmp,pi/2,c(-lim,-lim)),
                   rotate(tmp,pi,c(-lim,-lim)),
                   rotate(tmp,-pi/2,c(lim,lim)),
                   rotate(tmp,0,c(0,0)),
                   rotate(tmp,-pi/2,c(-lim,-lim)),
                   rotate(tmp,pi,c(lim,lim)),
                   rotate(tmp,pi/2,c(lim,lim)),
                   rotate(tmp,pi,c(lim,-lim))),
              3,3)
  out.coords<-do.call(rbind,as.vector(tmp))
  tmp.prj<-terra::vect(as.matrix(out.coords[,c("geom","part","lon","lat","hole")]),
                       geom.type)
  
  nn<-length(x)
  out.prj<-NULL
  for(i in 1:nn){
    if(geom.type=="polygons"){
      tmp<-terra::aggregate(terra::buffer(tmp.prj[(0:8)*nn+i],width=exp.fudge.fac))
    }else{
      tmp<-terra::aggregate(tmp.prj[(0:8)*nn+i])
    }
    terra::values(tmp)<-terra::values(x)[i,]
    if(is.null(out.prj)){
      out.prj<-tmp
    }else{
      out.prj<-rbind(out.prj,tmp)
    }
  }
  
  wd<-amount*2*lim
  rectangle<-cbind(c(-lim-wd,
                     -lim-wd,
                     lim+wd,
                     lim+wd),
                   c(-lim-wd,
                     lim+wd,
                     lim+wd,
                     -lim-wd))
  terra::crop(out.prj,terra::vect(rectangle,"polygons"))
}

patch.corners<-function(x,
                        NA.radius=1e6+10,
                        SA.radius=1e6+10,
                        Asia.radius=1e6+10){
  lim<-11825474
  if(NA.radius>0){
    NA.patch1<-terra::vect(
      cbind(c(-lim,seq(-lim,-lim+NA.radius,length.out=100)),
            c(lim,seq(lim-NA.radius,lim,length.out=100))),
      "polygons")
    inds<-which(terra::is.related(x,NA.patch1,"intersects"))
    if(length(inds)>0){
      for(i in inds){
        tmp<-terra::aggregate(terra::union(x[i],NA.patch1))
        terra::values(tmp)<-terra::values(x)[i,]
        x<-rbind(x,tmp)
      }
      x<-x[-inds]
    }
    NA.patch2<-terra::vect(
      cbind(c(lim,seq(lim,lim-NA.radius,length.out=100)),
            c(-lim,seq(-lim+NA.radius,-lim,length.out=100))),
      "polygons")
    inds<-which(terra::is.related(x,NA.patch2,"intersects"))
    if(length(inds)>0){
      for(i in inds){
        tmp<-terra::aggregate(terra::union(x[i],NA.patch2))
        terra::values(tmp)<-terra::values(x)[i,]
        x<-rbind(x,tmp)
      }
      x<-x[-inds]
    }
  }
  if(SA.radius>0){
    SA.patch<-
      terra::vect(
        cbind(c(-lim,seq(-lim,-lim+SA.radius,length.out=100)),
              c(-lim,seq(-lim+SA.radius,-lim,length.out=100))),
        "polygons")
    inds<-which(terra::is.related(x,SA.patch,"intersects"))
    if(length(inds)>0){
      for(i in inds){
        tmp<-terra::aggregate(terra::union(x[i],SA.patch))
        terra::values(tmp)<-terra::values(x)[i,]
        x<-rbind(x,tmp)
      }
      x<-x[-inds]
    }
  }
  if(Asia.radius>0){
    Asia.patch<-
      terra::vect(
        cbind(c(lim,seq(lim,lim-Asia.radius,length.out=100)),
              c(lim,seq(lim-Asia.radius,lim,length.out=100))),
        "polygons")
    inds<-which(terra::is.related(x,Asia.patch,"intersects"))
    if(length(inds)>0){
      for(i in inds){
        tmp<-terra::aggregate(terra::union(x[i],Asia.patch))
        terra::values(tmp)<-terra::values(x)[i,]
        x<-rbind(x,tmp)
      }
      x<-x[-inds]
    }
  }
  x
}

gen.gridlines<-function(mins,maxs,n,res=n*10){
  mins<-rep(mins,length.out=2)
  maxs<-rep(maxs,length.out=2)
  n<-rep(n,length.out=2)
  res<-rep(res,length.out=2)
  x.seq1<-seq(mins[1],maxs[1],length.out=n[1]+1)
  x.seq2<-seq(mins[1],maxs[1],length.out=res[1]+1)
  y.seq1<-seq(mins[2],maxs[2],length.out=n[2]+1)
  y.seq2<-seq(mins[2],maxs[2],length.out=res[2]+1)
  c(lapply(x.seq1,function(ii) cbind(ii,y.seq2)),
    lapply(y.seq1,function(ii) cbind(x.seq2,ii)))
}

#fix polygon seams resulting from intersecting with -180/180 longitude line
fix.pacific.seams<-function(x,
                            intersection.check.tol=1000,
                            snap.tol=5000,
                            exp.fudge.fac=300){
  #map -180/180 longitude line onto spilhaus projection
  seam.lonlat<-cbind(180,seq(-90,90,length.out=1000))
  seam.xy<-from_lonlat_to_spilhaus_xy(seam.lonlat[,1],seam.lonlat[,2])
  #split out separate components (thankfully one half always has positive y coords, the other always negative)
  seam.xy<-lapply(split(seam.xy,seam.xy[,2]>0),matrix,ncol=2)
  #fatten up to make sure intersection test gets everything
  seam<-terra::buffer(terra::vect(seam.xy,"lines"),intersection.check.tol)
  #find geometries overlapping with seam
  inds<-which(terra::is.related(x,seam,"intersects"))
  
  #for each problem geometry...
  #first diaggregate into individual pieces
  #then snap their borders together (within some tolerance)
  #then buffer the shapes slightly
  #(not ideal, but found this little bit of fudging necessary to cover all cases)
  #then finally reaggregate into a single, dissolving any overlapping borders
  #9/28-->calling aggregate without "by" arguments results in dropping field attributes
  #now fixed by manually re-adding fields to aggregated shapes
  for(i in inds){
    tmp<-terra::aggregate(
      terra::buffer(
        terra::snap(
          terra::disagg(x[i]),
          tol=snap.tol),
        exp.fudge.fac))
    terra::values(tmp)<-terra::values(x)[i,]
    x<-rbind(x,
             tmp)
  }
  #eliminate the old geometry elements (R didn't seem to like directly replacing the indices for whatever reason)
  x<-x[-inds]
  
  x
}

#splits any polygons overlapping with "seams" of Spilhaus projection
preproc.split<-function(x,
                        NA.radius=1e6,
                        SA.radius=1e6,
                        Asia.radius=1e6){
  
  #max coordinate in spilhaus projection
  lim<-11825474
  
  #new idea
  lim.seq<-seq(-lim,lim,length.out=1000)
  seam.xy<-rbind(cbind(-lim,lim.seq),
                 cbind(lim.seq,lim),
                 cbind(lim,-lim.seq),
                 cbind(-lim.seq,-lim))
  seam.lonlat<-from_spilhaus_xy_to_lonlat(seam.xy[,1],seam.xy[,2])
  
  #eliminate the cases that failed (projecting the border can of course make things go wonky)
  seam.lonlat<-seam.lonlat[complete.cases(seam.lonlat),]
  
  #new idea--also remove anything too close to intersection point?
  #where do corners break down?
  #less than 53899,53899
  #less than 675000,0
  #less than 0,815000
  # 815000/675000*53899
  #triangle from 675000,0 to 0,815000 encompasses 53899,53899
  #so removing NAs should be taking care of things...
  #except things change a bit depending on corner, so...
  #1000000 is a safe bet
  # corn.pat<-from_spilhaus_xy_to_lonlat(c(-lim+1000000*cos(seq(-pi/2,0,length.out=100)),
  #                   lim+1000000*cos(seq(pi/2,pi,length.out=100))),
  #             c(lim+1000000*sin(seq(-pi/2,0,length.out=100)),
  #               -lim+1000000*sin(seq(pi/2,pi,length.out=100))))
  #maybe try triangles?
  #yeah, this seems to work better!
  NA.patch<-from_spilhaus_xy_to_lonlat(c(seq(-lim,-lim+NA.radius,length.out=100),
                                         seq(lim,lim-NA.radius,length.out=100)),
                                       c(seq(lim-NA.radius,lim,length.out=100),
                                         seq(-lim+NA.radius,-lim,length.out=100)))
  SA.patch<-from_spilhaus_xy_to_lonlat(seq(-lim,-lim+SA.radius,length.out=100),
                                         seq(-lim+SA.radius,-lim,length.out=100))
  Asia.patch<-from_spilhaus_xy_to_lonlat(seq(lim-Asia.radius,lim,length.out=100),
                                           seq(lim,lim-Asia.radius,length.out=100))
  
  
  
  #split up at seam of typical mercator projection (i.e., meridian in the Pacific)
  seam.lonlat<-lapply(split(seam.lonlat,seam.lonlat[,1]>0),matrix,ncol=2)
  
  #find min pair of first split 
  nn<-nrow(seam.lonlat[[1]])
  min1<-which.min((seam.lonlat[[1]]+seam.lonlat[[1]][c(nn,1:(nn-1)),1])/2)
  min1<-c(min1-2,min1-1)%%nn+1
  # plot(seam.lonlat[[1]],ylim=c(64.9,65.1))
  # points(seam.lonlat[[1]][min1,],pch=16) #perfect
  #and max pair of second split
  nn<-nrow(seam.lonlat[[2]])
  max2<-which.max((seam.lonlat[[2]]+seam.lonlat[[2]][c(nn,1:(nn-1)),1])/2)
  max2<-c(max2-2,max2-1)%%nn+1
  # plot(seam.lonlat[[2]],ylim=c(64.9,65.1))
  # points(seam.lonlat[[2]][max2,],pch=16) #perfect
  tmp.xx<-seam.lonlat[[1]][,1]
  tmp.yy<-seam.lonlat[[1]][,2]
  if(min1[1]<min1[2]){
    tmp.xx<-append(tmp.xx,seam.lonlat[[2]][max2[2],1]-360,min1[1])
    tmp.yy<-append(tmp.yy,seam.lonlat[[2]][max2[2],2],min1[1])
    tmp.xx<-append(tmp.xx,seam.lonlat[[2]][max2[1],1]-360,min1[2])
    tmp.yy<-append(tmp.yy,seam.lonlat[[2]][max2[1],2],min1[2])
  }else{
    tmp.xx<-append(tmp.xx,seam.lonlat[[2]][max2[1],1]-360,min1[2]-1)
    tmp.yy<-append(tmp.yy,seam.lonlat[[2]][max2[1],2],min1[2]-1)
    tmp.xx<-append(tmp.xx,seam.lonlat[[2]][max2[2],1]-360,min1[1]+1)
    tmp.yy<-append(tmp.yy,seam.lonlat[[2]][max2[2],2],min1[1]+1)
  }
  tmp.seam<-cbind(tmp.xx,tmp.yy)
  
  #how to fix self-intersection?
  #check every edge? that would take too long right?
  #nope, not if you vectorize!
  #however, this all ended up being of limited utility compared to just buffering with 0 width
  #didn't allow me to use lower widths...
  #giving up for now, but code may prove useful in future...
  # cross<-function(p1s,p2s){
  #   p1s[,1]*p2s[,2]-p1s[,2]*p2s[,1]
  # }
  # check.intersect<-function(s1s,e1s,s2s,e2s){
  #   
  #   #check for cases including the same points...
  #   tmp.check<-cbind(s1s[,1]==s2s[,1],
  #                    s1s[,1]==e2s[,1],
  #                    e1s[,1]==s2s[,1],
  #                    e1s[,1]==e2s[,1])&
  #     cbind(s1s[,2]==s2s[,2],
  #           s1s[,2]==e2s[,2],
  #           e1s[,2]==s2s[,2],
  #           e1s[,2]==e2s[,2])
  #   tmp.check<-tmp.check[,1]|tmp.check[,2]|tmp.check[,3]|tmp.check[,4]
  #   
  #   e1s<-e1s-s1s
  #   e2s<-e2s-s2s
  #   tmp1s<-s2s-s1s
  #   tmp2s<-cross(e1s,e2s)
  #   ts<-cross(tmp1s,e2s/tmp2s)
  #   us<-cross(tmp1s,e1s/tmp2s)
  #   which(ts>0&ts<1&us>0&us<1&!tmp.check)
  # }
  # tmp.seq<-seq_len(nrow(tmp.seam))
  # inds1<-rep(tmp.seq,rev(tmp.seq)-1)
  # inds2<-do.call(c,lapply(tmp.seq[-nrow(tmp.seam)]+1,seq,to=nrow(tmp.seam)))
  # s1s<-tmp.seam[(inds1-2)%%nrow(tmp.seam)+1,]
  # e1s<-tmp.seam[inds1,]
  # s2s<-tmp.seam[(inds2-2)%%nrow(tmp.seam)+1,]
  # e2s<-tmp.seam[inds2,]
  # intersection<-check.intersect(s1s,e1s,s2s,e2s)
  # #this may be a stupid, hacky way to do it...but...
  # while(length(intersection)>0){
  #   inds<-range(c(inds1[intersection],inds2[intersection]))
  #   tmp.seam[inds[1]:(inds[2]-1),]<-tmp.seam[(inds[2]-1):inds[1],]
  #   s1s<-tmp.seam[(inds1-2)%%nrow(tmp.seam)+1,]
  #   e1s<-tmp.seam[inds1,]
  #   s2s<-tmp.seam[(inds2-2)%%nrow(tmp.seam)+1,]
  #   e2s<-tmp.seam[inds2,]
  #   intersection<-check.intersect(s1s,e1s,s2s,e2s)
  # }
  
  tmp.xx<-seam.lonlat[[2]][,1]
  tmp.yy<-seam.lonlat[[2]][,2]
  if(max2[1]<max2[2]){
    tmp.xx<-append(tmp.xx,seam.lonlat[[1]][min1[2],1]+360,max2[1])
    tmp.yy<-append(tmp.yy,seam.lonlat[[1]][min1[2],2],max2[1])
    tmp.xx<-append(tmp.xx,seam.lonlat[[1]][min1[1],1]+360,max2[2])
    tmp.yy<-append(tmp.yy,seam.lonlat[[1]][min1[1],2],max2[2])
  }else{
    tmp.xx<-append(tmp.xx,seam.lonlat[[1]][min1[1],1]+360,max2[2]-1)
    tmp.yy<-append(tmp.yy,seam.lonlat[[1]][min1[1],2],max2[2]-1)
    tmp.xx<-append(tmp.xx,seam.lonlat[[1]][min1[2],1]+360,max2[1]+1)
    tmp.yy<-append(tmp.yy,seam.lonlat[[1]][min1[2],2],max2[1]+1)
  }
  seam.lonlat[[1]]<-tmp.seam
  seam.lonlat[[2]]<-cbind(tmp.xx,tmp.yy)
  
  seam<-terra::buffer(terra::vect(seam.lonlat,"polygons",
                                  crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"),
                      width=0)
  seam1<-terra::aggregate(terra::union(seam[1],terra::vect(list(NA.patch,SA.patch),"polygons")))
  seam2<-terra::aggregate(terra::union(seam[2],terra::vect(Asia.patch,"polygons")))
  seam<-rbind(seam1,seam2)
  #use seam to erase bits of geometry at border
  terra::erase(x,seam)
}

# auxiliary function 1: fetch 5km resolution land mask
download_land_mask = function() {
  myInfo = rerddap::info('NOAA_DHW_monthly', url='https://coastwatch.pfeg.noaa.gov/erddap/')
  myData = rerddap::griddap(datasetx = myInfo,
                            fields = "mask",
                            latitude = c(-89.975, 89.975),
                            longitude = c(-179.975, 179.975),
                            time = c('2024-01-01', '2024-01-01'))
  da = ncdf4::nc_open(myData$summary$filename)
  return(da)
}


# auxiliary function 2: determine whether arbitrary lon-lat points fall on land
extract_mask = function(da, lonlat) {
  lons = sort(ncdf4::ncvar_get(da, "longitude"))
  lats = sort(ncdf4::ncvar_get(da, "latitude"))
  dlon = lons[2] - lons[1]
  dlat = lats[2] - lats[1]
  ln = as.integer(round((lonlat[,1] - lons[1]) / dlon) + 1)
  la = as.integer(round((lonlat[,2] - lats[1]) / dlat) + 1)
  msk = ncdf4::ncvar_get(da, "mask", start = c(1,1,1), count = c(7200, 3600,1))
  msk = msk[,ncol(msk):1]
  chunk = msk[ln + (la - 1) * dim(msk)[1]]
  chunk[chunk == 1] = NA # in the original mask, 1 denotes land
  chunk[chunk == 4] = 0 # in the original mask, 4 denotes ice
  return(chunk)
}


from_lonlat_to_spilhaus_xy <- function(longitude, latitude){

  # constants (https://github.com/OSGeo/PROJ/issues/1851)
  e = sqrt(0.00669438)
  lat_center_deg = -49.56371678
  lon_center_deg = 66.94970198
  azimuth_deg = 40.17823482

  # parameters derived from constants
  lat_center_rad = lat_center_deg * pi / 180
  lon_center_rad = lon_center_deg * pi / 180
  azimuth_rad = azimuth_deg * pi / 180
  conformal_lat_center = -pi / 2 + 2 * atan(
    tan(pi/4 + lat_center_rad/2) *
      ((1 - e * sin(lat_center_rad)) / (1 + e * sin(lat_center_rad))) ^ (e / 2)
  )
  alpha = -asin(cos(conformal_lat_center) * cos(azimuth_rad))
  lambda_0 = lon_center_rad + atan2(tan(azimuth_rad), -sin(conformal_lat_center))
  beta = pi + atan2(-sin(azimuth_rad), -tan(conformal_lat_center))

  # coordinates in radians
  lon = longitude * pi / 180
  lat = latitude * pi / 180

  # conformal latitude, in radians
  lat_c = -pi / 2 + 2 * atan(
    tan(pi/4 + lat/2) * ((1 - e * sin(lat)) / (1 + e * sin(lat))) ^ (e / 2)
  )

  # transformed lat and lon, in degrees
  lat_s = 180 / pi * asin(sin(alpha) * sin(lat_c) - cos(alpha) * cos(lat_c) * cos(lon - lambda_0))
  lon_s = 180 / pi * (
    beta + atan2(
      cos(lat_c) * sin(lon - lambda_0),
      (sin(alpha) * cos(lat_c) * cos(lon - lambda_0) + cos(alpha) * sin(lat_c))
    )
  )

  # projects transformed coordinates onto plane (Adams World in a Square II)
  adams_ws2 = "+proj=adams_ws2 +no_defs +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m"
  projected = sf_project(from=sf::st_crs(4326), to=adams_ws2, pts=cbind(lon_s, lat_s))
  adams_x = projected[,1]
  adams_y = projected[,2]
  spilhaus_x = -(adams_x + adams_y) / sqrt(2)
  spilhaus_y = (adams_x - adams_y) / sqrt(2)

  return(cbind(spilhaus_x, spilhaus_y)) #, adams_x, adams_y, lon_s, lat_s))
}

make_spilhaus_xy_gridpoints <- function(spilhaus_res=1000) {
  # regular grid of points in Spilhaus map
  extreme = 11825474
  m = seq(-extreme, extreme, len=spilhaus_res)
  spilhaus_df = expand.grid(x=m, y=m)
  return(spilhaus_df)
}

flip1 = function(x) {
  res = as.integer(sqrt(length(x)))
  return(c(t(matrix(nrow=res, ncol=res, x)[,res:1])))
}
flip2 = function(x) {
  res = as.integer(sqrt(length(x)))
  return(c(t(matrix(nrow=res, ncol=res, x))[,res:1]))
}



pretify_spilhaus_df <- function(spilhaus_df) {

    spilhaus_x = spilhaus_df$x
    spilhaus_y = spilhaus_df$y
    spilhaus_z = spilhaus_df$z
    spilhaus_l = spilhaus_df$l

    extreme = 11825474

    # augmented grid points
    aug_x = c(
      spilhaus_x,
      spilhaus_x - 2 * extreme,
      spilhaus_x + 2 * extreme,
      spilhaus_x,
      spilhaus_x
    )
    aug_y = c(
      spilhaus_y,
      spilhaus_y,
      spilhaus_y,
      spilhaus_y - 2 * extreme,
      spilhaus_y + 2 * extreme
    )
    aug_z = c(
      spilhaus_z,
      flip1(spilhaus_z),
      flip1(spilhaus_z),
      flip2(spilhaus_z),
      flip2(spilhaus_z)
    )
    aug_l = c(
      spilhaus_l,
      flip1(spilhaus_l),
      flip1(spilhaus_l),
      flip2(spilhaus_l),
      flip2(spilhaus_l)
    )

    cutpoint = 1.1 * extreme
    keep = !(
      aug_l
      | (aug_x < -cutpoint)
      | (aug_x > cutpoint)
      | (aug_y < -cutpoint)
      | (aug_y > cutpoint)
      | (aug_y > 1.089e7 - 0.176 * aug_x)
      | (aug_x < -0.984e7 - 0.565 * aug_y)
      | (aug_y < -1.378e7 + 0.46 * aug_x)
      | (aug_x > 1.274e7 + 0.172 * aug_y)
      | (aug_y > 1e7 - 0.5 * aug_x)
      | (aug_y > 2.3e7 + aug_x)
      | ((aug_y < 0.29e7) & (aug_x < -1.114e7))
      | ((aug_y < 0.39e7) & (aug_x < -1.17e7))
      | ((aug_y < -1.21e7) & (aug_x > 0.295e7))
      | ((aug_y < -1.2e7) & (aug_x > 0.312e7))
      | ((aug_y < -1.16e7) & (aug_x > 0.4e7))
      | ((aug_y < -1.11e7) & (aug_x > 0.45e7))
    )

    pretty_spilhaus_df = data.frame(
      x=aug_x[keep],
      y=aug_y[keep],
      z=aug_z[keep]
    )
    pretty_spilhaus_df = pretty_spilhaus_df[!duplicated(pretty_spilhaus_df[c(1,2)]),]

    return(pretty_spilhaus_df)
}

from_spilhaus_xy_to_lonlat <- function(spilhaus_x, spilhaus_y) {
  # constants
  e = sqrt(0.00669438)
  lat_center_deg = -49.56371678
  lon_center_deg = 66.94970198
  azimuth_deg = 40.17823482

  # parameters derived from constants
  lat_center_rad = lat_center_deg * pi / 180
  lon_center_rad = lon_center_deg * pi / 180
  azimuth_rad = azimuth_deg * pi / 180
  conformal_lat_center = -pi / 2 + 2 * atan(
    tan(pi/4 + lat_center_rad/2) *
      ((1 - e * sin(lat_center_rad)) / (1 + e * sin(lat_center_rad))) ^ (e / 2)
  )
  alpha = -asin(cos(conformal_lat_center) * cos(azimuth_rad))
  lambda_0 = lon_center_rad + atan2(tan(azimuth_rad), -sin(conformal_lat_center))
  beta = pi + atan2(-sin(azimuth_rad), -tan(conformal_lat_center))

  adams_x = (spilhaus_y - spilhaus_x) * sqrt(2) / 2
  adams_y = - (spilhaus_y + spilhaus_x) * sqrt(2) / 2

  adams_ws2 = "+proj=adams_ws2 +no_defs +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m"
  projection_fun = function(x, y) {
    tryCatch(
      sf_project(from=adams_ws2, to=sf::st_crs(4326), pts=c(x, y)),
      error=function(e) c(NA, NA)
    )
  }
  projected = sf_project(from=adams_ws2, to=sf::st_crs(4326), pts=cbind(adams_x, adams_y), keep = TRUE, warn = FALSE)
  lon_s = projected[,1]
  lat_s = projected[,2]

  #transformed coords in radians
  lon_s_rad = lon_s * pi / 180
  lat_s_rad = lat_s * pi / 180

  # conformal latitude
  lat_c = asin(sin(alpha) * sin(lat_s_rad) + cos(alpha) * cos(lat_s_rad) * cos(lon_s_rad - beta))

  # longitude, in radians
  lon = lambda_0 + atan2(
    cos(lat_s_rad) * sin(lon_s_rad - beta),
    sin(alpha) * cos(lat_s_rad) * cos(lon_s_rad - beta) - cos(alpha) * sin(lat_s_rad)
  )

  # latitude (iterative formula from https://mathworld.wolfram.com/ConformalLatitude.html)
  lat = lat_c
  for (i in 0:9) {
    lat = -0.5 * pi + 2 * atan(
      tan(pi / 4 + lat_c / 2) *
        ((1 + e * sin(lat)) / (1 - e * sin(lat))) ^ (e / 2)
    )
  }

  # coordinates in degrees
  longitude = ((lon * 180 / pi + 180) %% 360) - 180
  latitude = lat * 180 / pi

  return(cbind(longitude, latitude))
}


# auxiliary function to extract netCDF file from Coastwatch
download_sst_data = function(start_date, end_date) {
  myInfo = rerddap::info('NOAA_DHW_monthly', url='https://coastwatch.pfeg.noaa.gov/erddap/')
  myData = rerddap::griddap(datasetx = myInfo,
                            fields = "sea_surface_temperature",
                            latitude = c(-89.975, 89.975),
                            longitude = c(-179.975, 179.975),
                            time = c(start_date, end_date))
  da = ncdf4::nc_open(myData$summary$filename)
  return(da)
}


# auxiliary function to extract a list (length = number of instants) 
#of snapshots (7200x3600 pixels) from the netCDF file
extract_sst_data = function(da, lonlat) {
  lons = sort(ncdf4::ncvar_get(da, "longitude"))
  lats = sort(ncdf4::ncvar_get(da, "latitude"))
  times = ncdf4::ncvar_get(da, "time")
  
  dlon = lons[2] - lons[1]
  dlat = lats[2] - lats[1]
  ln = as.integer(round((lonlat[,1] - lons[1]) / dlon) + 1)
  la = as.integer(round((lonlat[,2] - lats[1]) / dlat) + 1)
  get_chunk = function(timepoint) {
    sst = ncdf4::ncvar_get(da, "sea_surface_temperature",
                           start = c(1,1,timepoint), count = c(7200, 3600, 1))
    sst = sst[,ncol(sst):1]
    chunk = sst[ln + (la - 1) * dim(sst)[1]]
    return(chunk)
  }
  sst_data = mapply(1:length(times), FUN=function(timepoint) {get_chunk(timepoint)})
  return(sst_data)
}
