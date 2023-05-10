#' A Rcove Function
#'
#' This function allows you to Calculate Dispersal Range of Propagules using Rotated Components of Water Velocity
#' @name rcove
#' @param theta Angle about which to rotate current velocities. Degrees. Numeric
#' @param velocity Velocity time series. Vector consisting of two columns u,v.
#' @param release_pts Start locations for propagules. sf points
#' @param PD Planktonic propagule duration of species of interest. Numeric
#' @param CPD Competent propagule duration of species of interest. Numeric
#' @param land Land to be subtracted from final concatenated polygons. sf polygon
#' @param proj Projected crs in units m.
#' @keywords dispersal
#' @import purrr
#' @import units
#' @import tidyr
#' @import dplyr
#' @import rmapshaper
#' @import nngeo
#' @import sf
#' @export
#' @examples
#' rcove()

# library(rmapshaper)
# library(dplyr)
# library(tidyr)
# library(units)
# library(purrr)
# library(nngeo)
# library(sf)

rcove <- function(theta, velocity, release_pts, PD, CPD, land, proj){

# Prepare velocity time series ====

#Format velocity time series
theta <- units::set_units(theta, "degrees") #set units of theta
release_pts <- sf::st_transform(release_pts, proj)
land <- sf::st_transform(land, proj)
rotation.matrix <- function(a) matrix(data = c(cos(a), -sin(a), sin(a), cos(a)), nrow = 2, ncol = 2, byrow = TRUE)
rotated.velocity <- as.matrix(velocity)%*%rotation.matrix(theta) %>% #Rotate velocity time series about angle theta
  as.data.frame() %>%
  dplyr::rename("along.shore" = "V1",
         "cross.shore" = "V2")

#Calculate mean of each component of velocity
components_of_velocity <- rotated.velocity %>%
  tidyr::pivot_longer(cols = c(1,2),
               names_to = "direction",
               values_to = "velocity") %>%
  dplyr::mutate(sign = dplyr::case_when(
    velocity > 0 ~ "positive",
    TRUE ~ "negative"
  )) %>%
  dplyr::group_by(direction, sign) %>%
  dplyr::summarise(speed = round(mean(abs(velocity)), digits = 4)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(dir = dplyr::case_when(
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
    return(sf::st_sfc(sf::st_linestring(matrix(c(sf::st_coordinates(origin),sf::st_coordinates(origin)+c(x,y)),2,2,byrow = TRUE))))
  } else {
    return(c(x,y))
  }
}

# Calculate dispersal area for PD ====

#Apply rule of thumb for dispersal distance (distance = average velocity X time (PD))
foci <- components_of_velocity %>%
  dplyr::mutate(distance = units::set_units(as.numeric(speed)*86400*PD, m))

angles <- units::set_units(c(270, 90, 180, 360), degrees)
rotated_angles <- as.data.frame(angles - theta) %>%
  dplyr::rename("angle" = "angles - theta")
foci_lengths <- cbind(foci, rotated_angles)

foci_lengths$angle <- as.numeric(foci_lengths$angle)
foci_lengths$distance <- as.numeric(foci_lengths$distance)

foci_vectors <- release_pts %>%
  dplyr::mutate(id = as.factor(row.names(.))) %>%
  merge(foci_lengths) %>%
  dplyr::mutate(cl=purrr::pmap(list(angle,distance),
                 function(angle,distance) compassline(angle,distance))) %>%
  dplyr::mutate(origin=geometry,
         geometry=purrr::pmap(list(geometry,cl),
                       function(geometry,cl) sf::st_point(sf::st_coordinates(geometry)+cl)) %>%
           sf::st_sfc(crs=proj)) %>%
      dplyr::select(-origin, -cl)

##calculate new centroids based on the 4 kite points
new_release_pts <- foci_vectors %>%
  dplyr::select(id,dir, geometry) %>%
  dplyr::rowwise %>%
  dplyr::mutate(x = sf::st_coordinates(geometry)[1],
         y = sf::st_coordinates(geometry)[2]) %>%
  sf::st_drop_geometry() %>%
  tidyr::pivot_wider(names_from = dir, values_from = c(x,y)) %>%
  dplyr::mutate(x_centroid = (x_east + x_west + x_south + x_north)/4,
         y_centroid = (y_east + y_west + y_south + y_north)/4) %>%
  sf::st_as_sf(coords = c(10,11), crs = proj) %>%
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
ex <- as.numeric(sf::st_distance(east_foci, new_release_pts, by_element = TRUE))
wx <- as.numeric(sf::st_distance(west_foci, new_release_pts, by_element = TRUE))
ny <- as.numeric(sf::st_distance(north_foci, new_release_pts, by_element = TRUE))
sy <- as.numeric(sf::st_distance(south_foci, new_release_pts, by_element = TRUE))


##Loop for concentating four ellipse quadrats starts here
m <- 1

for(m in 1:4){
  if (m == 1){
    ellipse <- sf::st_ellipse(pnt = new_release_pts, ex = wx/2, ey = sy/2, res = 60) %>%
      sf::st_as_sf()
    } else if (m == 2){
      ellipse <- sf::st_ellipse(pnt = new_release_pts, ex = ex/2, ey = sy/2, res = 60) %>%
        sf::st_as_sf()
      } else if (m == 3){
        ellipse <- sf::st_ellipse(pnt = new_release_pts, ex = ex/2, ey = ny/2, res = 60) %>%
          sf::st_as_sf()
        } else {
          ellipse <- sf::st_ellipse(pnt = new_release_pts, ex = wx/2, ey = ny/2, res = 60) %>%
            sf::st_as_sf()
        }

sf::st_crs(ellipse) <- proj

##Rotate ellipse around newrelease_pts
rotated_ellipse <- (sf::st_geometry(ellipse) - sf::st_geometry(new_release_pts)) * rotation.matrix(-theta) + sf::st_geometry(new_release_pts)
sf::st_crs(rotated_ellipse) <- proj

#Convert to linestring and cut into relative quadrat
ellipse_points <- rotated_ellipse %>%
  sf::st_cast("MULTIPOINT", ids = id) %>%
  as.data.frame() %>%
  sf::st_as_sf() %>%
  dplyr::mutate(id = row_number()) %>%
  sf::st_cast("POINT", ids = id) %>%
  dplyr::mutate(ellipse_id = rep(seq(1:60), times = nrow(new_release_pts)))

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
concatenated_ellipse_PD <- rbind(cut_sw, cut_se, cut_ne, cut_nw) %>%
  dplyr::group_by(id) %>%
  dplyr::summarise() %>%
  sf::st_convex_hull()

# Calculate dispersal area for CPD ====

#Apply rule of thumb for dispersal distance (distance = average velocity X time (CPD))
foci <- components_of_velocity %>%
  dplyr::mutate(distance = set_units(as.numeric(speed)*86400*CPD, m))

angles <- set_units(c(270, 90, 180, 360), degrees)
rotated_angles <- as.data.frame(angles - theta) %>%
  dplyr::rename("angle" = "angles - theta")
foci_lengths <- cbind(foci, rotated_angles)

foci_lengths$angle <- as.numeric(foci_lengths$angle)
foci_lengths$distance <- as.numeric(foci_lengths$distance)

foci_vectors <- release_pts %>%
  dplyr::mutate(id = as.factor(row.names(.))) %>%
  merge(foci_lengths) %>%
  dplyr::mutate(cl=purrr::pmap(list(angle,distance),
                 function(angle,distance) compassline(angle,distance))) %>%
  dplyr::mutate(origin=geometry,
         geometry=purrr::pmap(list(geometry,cl),
                       function(geometry,cl) st_point(st_coordinates(geometry)+cl)) %>%
           st_sfc(crs=proj)) %>%
  dplyr::select(-origin, -cl)

##calculate new centroids based on the 4 kite points
new_release_pts <- foci_vectors %>%
  dplyr::select(id,dir, geometry) %>%
  dplyr::rowwise %>%
  dplyr::mutate(x = sf::st_coordinates(geometry)[1],
         y = sf::st_coordinates(geometry)[2]) %>%
  sf::st_drop_geometry() %>%
  tidyr::pivot_wider(names_from = dir, values_from = c(x,y)) %>%
  dplyr::mutate(x_centroid = (x_east + x_west + x_south + x_north)/4,
         y_centroid = (y_east + y_west + y_south + y_north)/4) %>%
  sf::st_as_sf(coords = c(10,11), crs = proj) %>%
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
ex <- as.numeric(sf::st_distance(east_foci, new_release_pts, by_element = TRUE))
wx <- as.numeric(sf::st_distance(west_foci, new_release_pts, by_element = TRUE))
ny <- as.numeric(sf::st_distance(north_foci, new_release_pts, by_element = TRUE))
sy <- as.numeric(sf::st_distance(south_foci, new_release_pts, by_element = TRUE))


##Loop for concentating four ellipse quadrats starts here
m <- 1

for(m in 1:4){
  if (m == 1){
    ellipse <- sf::st_ellipse(pnt = new_release_pts, ex = wx/2, ey = sy/2, res = 60) %>%
      sf::st_as_sf()
  } else if (m == 2){
    ellipse <- sf::st_ellipse(pnt = new_release_pts, ex = ex/2, ey = sy/2, res = 60) %>%
      sf::st_as_sf()
  } else if (m == 3){
    ellipse <- sf::st_ellipse(pnt = new_release_pts, ex = ex/2, ey = ny/2, res = 60) %>%
      sf::st_as_sf()
  } else {
    ellipse <- sf::st_ellipse(pnt = new_release_pts, ex = wx/2, ey = ny/2, res = 60) %>%
      sf::st_as_sf()
  }

  st_crs(ellipse) <- proj

  ##Rotate ellipse around newrelease_pts
  rotated_ellipse <- (sf::st_geometry(ellipse) - sf::st_geometry(new_release_pts)) * rotation.matrix(-theta) + sf::st_geometry(new_release_pts)
  sf::st_crs(rotated_ellipse) <- proj

  #Convert to linestring and cut into relative quadrat
  ellipse_points <- rotated_ellipse %>%
    sf::st_cast("MULTIPOINT", ids = id) %>%
    as.data.frame() %>%
    sf::st_as_sf() %>%
    dplyr::mutate(id = row_number()) %>%
    sf::st_cast("POINT", ids = id) %>%
    dplyr::mutate(ellipse_id = rep(seq(1:60), times = nrow(new_release_pts)))

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
concatenated_ellipse_CPD <- rbind(cut_sw, cut_se, cut_ne, cut_nw) %>%
  dplyr::group_by(id) %>%
  dplyr::summarise() %>%
  sf::st_convex_hull()

CPD_ellipse <- st_difference(concatenated_ellipse_PD, concatenated_ellipse_CPD) |>
  dplyr::filter(id == id.1) |>
  rmapshaper::ms_erase(sf::st_simplify(land, dTolerance = 1000),
           remove_slivers = TRUE)

return(CPD_ellipse)

}
