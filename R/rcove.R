#' A Rcove Function
#'
#' This function allows you to Calculate Dispersal Range of Propagules using Rotated Components of Water Velocity
#' @name rcove
#' @param theta Angle about which to rotate current velocities. Numeric
#' @param velocity Velocity time series. Vector consisting of two columns u,v.
#' @param release_pts Start locations for propagules. sf object
#' @param PD Planktonic propagule duration of species of interest
#' @param CPD Competent propagule duration of species of interest
#' @param land Land to be subtracted from final concatenated polygons
#' @param proj Projected crs in units m.
#' @keywords dispersal
#' @export
#' @examples
#' rcove()

library(rmapshaper)
library(dplyr)
library(tidyr)
library(units)
library(purrr)
library(nngeo)
library(sf)

rcove <- function(theta, velocity, release_pts, PD, CPD, land, proj){

# Prepare velocity time series ====

#Format velocity time series
theta <- set_units(theta, "degrees") #set units of theta
rotation.matrix <- function(a) matrix(data = c(cos(a), -sin(a), sin(a), cos(a)), nrow = 2, ncol = 2, byrow = TRUE)
rotated.velocity <- as.matrix(velocity)%*%rotation.matrix(theta) %>% #Rotate velocity time series about angle theta
  as.data.frame() %>%
  rename("along.shore" = "V1",
         "cross.shore" = "V2")

#Calculate mean of each component of velocity
components_of_velocity <- rotated.velocity %>%
  pivot_longer(cols = c(1,2),
               names_to = "direction",
               values_to = "velocity") %>%
  mutate(sign = case_when(
    velocity > 0 ~ "positive",
    TRUE ~ "negative"
  )) %>%
  group_by(direction, sign) %>%
  summarise(speed = round(mean(abs(velocity)), digits = 4)) %>%
  ungroup() %>%
  mutate(dir = case_when(
    direction == "cross.shore" & sign == "positive" ~ "north",
    direction == "along.shore" & sign == "positive" ~ "east",
    direction == "cross.shore" & sign == "negative" ~ "south",
    TRUE ~ "west"
  )) %>%
  dplyr::select(-direction, -sign)

compassline <- function(theta,h,origin=NA){
  # browser()
  if(theta%%90<45){
    o <- sin(theta%%45*pi/180)*h
    a <- cos(theta%%45*pi/180)*h
  } else {
    o <- sin((45-theta%%45)*pi/180)*h
    a <- cos((45-theta%%45)*pi/180)*h
  }
  if(theta%/%90==0){
    if(theta%%90<45){
      x=o
      y=a
    } else {
      x=a
      y=o
    }
  }else if(theta%/%90==1){
    if(theta%%90>45){
      x=o
      y=-a
    } else {
      x=a
      y=-o
    }
  }else if(theta%/%90==2){
    if(theta%%90<45){
      x=-o
      y=-a
    } else {
      x=-a
      y=-o
    }
  }else if(theta%/%90==3){
    if(theta%%90>45){
      x=-o
      y=a
    } else {
      x=-a
      y=o
    }
  }
  if(!is.na(origin)){
    return(st_sfc(st_linestring(matrix(c(st_coordinates(origin),st_coordinates(origin)+c(x,y)),2,2,byrow = TRUE))))
  } else {
    return(c(x,y))
  }
}

# Calculate dispersal area for PD ====

#Apply rule of thumb for dispersal distance (distance = average velocity X time (PD))
foci <- components_of_velocity %>%
      mutate(distance = set_units(as.numeric(speed)*86400*PD, m))

angles <- set_units(c(270, 90, 180, 360), degrees)
rotated_angles <- as.data.frame(angles - theta) %>%
  rename("angle" = "angles - theta")
foci_lengths <- cbind(foci, rotated_angles)

foci_lengths$angle <- as.numeric(foci_lengths$angle)
foci_lengths$distance <- as.numeric(foci_lengths$distance)

foci_vectors <- release_pts %>%
  mutate(id = as.factor(row.names(.))) %>%
  merge(foci_lengths) %>%
  mutate(cl=pmap(list(angle,distance),
                 function(angle,distance) compassline(angle,distance))) %>%
  mutate(origin=geometry,
         geometry=pmap(list(geometry,cl),
                       function(geometry,cl) st_point(st_coordinates(geometry)+cl)) %>%
           st_sfc(crs=proj)) %>%
      dplyr::select(-origin, -cl)

##calculate new centroids based on the 4 kite points
new_release_pts <- foci_vectors %>%
  dplyr::select(id,dir, geometry) %>%
  rowwise %>%
  mutate(x = st_coordinates(geometry)[1],
         y = st_coordinates(geometry)[2]) %>%
  st_drop_geometry() %>%
  pivot_wider(names_from = dir, values_from = c(x,y)) %>%
  mutate(x_centroid = (x_east + x_west + x_south + x_north)/4,
         y_centroid = (y_east + y_west + y_south + y_north)/4) %>%
  st_as_sf(coords = c(10,11), crs = proj) %>%
  dplyr::select(id, geometry)

##calculates the distance between the new centroids and N,E,S,W
west_foci <- foci_vectors %>%
  dplyr::filter(dir == "west")
east_foci <- foci_vectors %>%
  dplyr::filter(dir == "east")
north_foci <- foci_vectors %>%
  dplyr::filter(dir == "north")
south_foci <- foci_vectors %>%
  dplyr::filter(dir == "south")
ex <- as.numeric(st_distance(east_foci, new_release_pts, by_element = TRUE))
wx <- as.numeric(st_distance(west_foci, new_release_pts, by_element = TRUE))
ny <- as.numeric(st_distance(north_foci, new_release_pts, by_element = TRUE))
sy <- as.numeric(st_distance(south_foci, new_release_pts, by_element = TRUE))


##Loop for concentating four ellipse quadrats starts here
m <- 1

for(m in 1:4){
  if (m == 1){
    ellipse <- st_ellipse(pnt = new_release_pts, ex = wx/2, ey = sy/2, res = 60) %>%
      st_as_sf()
    } else if (m == 2){
      ellipse <- st_ellipse(pnt = new_release_pts, ex = ex/2, ey = sy/2, res = 60) %>%
        st_as_sf()
      } else if (m == 3){
        ellipse <- st_ellipse(pnt = new_release_pts, ex = ex/2, ey = ny/2, res = 60) %>%
          st_as_sf()
        } else {
          ellipse <- st_ellipse(pnt = new_release_pts, ex = wx/2, ey = ny/2, res = 60) %>%
            st_as_sf()
        }

st_crs(ellipse) <- proj

##Rotate ellipse around newrelease_pts
rotated_ellipse <- (st_geometry(ellipse) - st_geometry(new_release_pts)) * rotation.matrix(-theta) + st_geometry(new_release_pts)
st_crs(rotated_ellipse) <- proj

#Convert to linestring and cut into relative quadrat
ellipse_points <- rotated_ellipse %>%
  st_cast("MULTIPOINT", ids = id) %>%
  as.data.frame() %>%
  st_as_sf() %>%
  mutate(id = row_number()) %>%
  st_cast("POINT", ids = id) %>%
  mutate(ellipse_id = rep(seq(1:60), times = nrow(new_release_pts)))

##this works but need to make an ifelse for the loop
if (m == 1){
  cut_sw <- ellipse_points %>%
    dplyr::filter(ellipse_id <= 15)
  } else if (m == 2){
    cut_se <- ellipse_points %>%
      dplyr::filter(ellipse_id > 15 & ellipse_id <=30)
    } else if (m == 3){
      cut_ne <- ellipse_points %>%
        dplyr::filter(ellipse_id > 30 & ellipse_id <=45)
      } else {
        cut_nw <- ellipse_points %>%
          dplyr::filter(ellipse_id > 45)
        }

m = m+1
} #end of creating 4 ellipse quadrats

##Start concatenation
concatenated_ellipse <- rbind(cut_sw, cut_se, cut_ne, cut_nw) %>%
  dplyr::group_by(id) %>%
  dplyr::summarise() %>%
  st_convex_hull() %>%
  ms_erase(land,
           remove_slivers = TRUE) %>%
  st_as_sf() %>%
  mutate(id = row_number()) %>%
  st_cast("POLYGON", ids = id)
concatenated_ellipse_PD <- concatenated_ellipse %>%
  mutate(area = st_area(concatenated_ellipse)) %>%
  group_by(id) %>%
  dplyr::filter(area == max(area)) %>%
  dplyr::select(-area)

# Calculate dispersal area for CPD ====

#Apply rule of thumb for dispersal distance (distance = average velocity X time (CPD))
foci <- components_of_velocity %>%
  mutate(distance = set_units(as.numeric(speed)*86400*CPD, m))

angles <- set_units(c(270, 90, 180, 360), degrees)
rotated_angles <- as.data.frame(angles - theta) %>%
  rename("angle" = "angles - theta")
foci_lengths <- cbind(foci, rotated_angles)

foci_lengths$angle <- as.numeric(foci_lengths$angle)
foci_lengths$distance <- as.numeric(foci_lengths$distance)

foci_vectors <- release_pts %>%
  mutate(id = as.factor(row.names(.))) %>%
  merge(foci_lengths) %>%
  mutate(cl=pmap(list(angle,distance),
                 function(angle,distance) compassline(angle,distance))) %>%
  mutate(origin=geometry,
         geometry=pmap(list(geometry,cl),
                       function(geometry,cl) st_point(st_coordinates(geometry)+cl)) %>%
           st_sfc(crs=proj)) %>%
  dplyr::select(-origin, -cl)

##calculate new centroids based on the 4 kite points
new_release_pts <- foci_vectors %>%
  dplyr::select(id,dir, geometry) %>%
  rowwise %>%
  mutate(x = st_coordinates(geometry)[1],
         y = st_coordinates(geometry)[2]) %>%
  st_drop_geometry() %>%
  pivot_wider(names_from = dir, values_from = c(x,y)) %>%
  mutate(x_centroid = (x_east + x_west + x_south + x_north)/4,
         y_centroid = (y_east + y_west + y_south + y_north)/4) %>%
  st_as_sf(coords = c(10,11), crs = proj) %>%
  dplyr::select(id, geometry)

##calculates the distance between the new centroids and N,E,S,W
west_foci <- foci_vectors %>%
  dplyr::filter(dir == "west")
east_foci <- foci_vectors %>%
  dplyr::filter(dir == "east")
north_foci <- foci_vectors %>%
  dplyr::filter(dir == "north")
south_foci <- foci_vectors %>%
  dplyr::filter(dir == "south")
ex <- as.numeric(st_distance(east_foci, new_release_pts, by_element = TRUE))
wx <- as.numeric(st_distance(west_foci, new_release_pts, by_element = TRUE))
ny <- as.numeric(st_distance(north_foci, new_release_pts, by_element = TRUE))
sy <- as.numeric(st_distance(south_foci, new_release_pts, by_element = TRUE))


##Loop for concentating four ellipse quadrats starts here
m <- 1

for(m in 1:4){
  if (m == 1){
    ellipse <- st_ellipse(pnt = new_release_pts, ex = wx/2, ey = sy/2, res = 60) %>%
      st_as_sf()
  } else if (m == 2){
    ellipse <- st_ellipse(pnt = new_release_pts, ex = ex/2, ey = sy/2, res = 60) %>%
      st_as_sf()
  } else if (m == 3){
    ellipse <- st_ellipse(pnt = new_release_pts, ex = ex/2, ey = ny/2, res = 60) %>%
      st_as_sf()
  } else {
    ellipse <- st_ellipse(pnt = new_release_pts, ex = wx/2, ey = ny/2, res = 60) %>%
      st_as_sf()
  }

  st_crs(ellipse) <- proj

  ##Rotate ellipse around newrelease_pts
  rotated_ellipse <- (st_geometry(ellipse) - st_geometry(new_release_pts)) * rotation.matrix(-theta) + st_geometry(new_release_pts)
  st_crs(rotated_ellipse) <- proj

  #Convert to linestring and cut into relative quadrat
  ellipse_points <- rotated_ellipse %>%
    st_cast("MULTIPOINT", ids = id) %>%
    as.data.frame() %>%
    st_as_sf() %>%
    mutate(id = row_number()) %>%
    st_cast("POINT", ids = id) %>%
    mutate(ellipse_id = rep(seq(1:60), times = nrow(new_release_pts)))

  ##this works but need to make an ifelse for the loop
  if (m == 1){
    cut_sw <- ellipse_points %>%
      dplyr::filter(ellipse_id <= 15)
  } else if (m == 2){
    cut_se <- ellipse_points %>%
      dplyr::filter(ellipse_id > 15 & ellipse_id <=30)
  } else if (m == 3){
    cut_ne <- ellipse_points %>%
      dplyr::filter(ellipse_id > 30 & ellipse_id <=45)
  } else {
    cut_nw <- ellipse_points %>%
      dplyr::filter(ellipse_id > 45)
  }

  m = m+1
} #end of creating 4 ellipse quadrats

##Start concatenation
concatenated_ellipse <- rbind(cut_sw, cut_se, cut_ne, cut_nw) %>%
  dplyr::group_by(id) %>%
  dplyr::summarise() %>%
  st_convex_hull() %>%
  ms_erase(land,
           remove_slivers = TRUE) %>%
  st_as_sf() %>%
  mutate(id = row_number()) %>%
  st_cast("POLYGON", ids = id)
concatenated_ellipse_CPD <- concatenated_ellipse %>%
  mutate(area = st_area(concatenated_ellipse)) %>%
  group_by(id) %>%
  dplyr::filter(area == max(area)) %>%
  dplyr::select(-area)

final_ellipse <<- st_difference(concatenated_ellipse_PD, concatenated_ellipse_CPD)

}
