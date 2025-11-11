rad2deg <- function(rad){
  deg = 180 * rad / pi
  return(deg)
}

#'
deg2rad <- function(deg){
  rad = pi * deg / 180
  return(rad)
}

#'
project_degrees <- function(lat, lon, destination_crs){

  latlon <- vect(data.frame(lat = lat, lon = lon),
                 geom = c('lat', 'lon'),
                 crs = 'epsg:4326')

  latlon_projected <- project(latlon, destination_crs)
  latlon_projected <- crds(latlon_projected)

  x <- latlon_projected[ ,1]
  y <- latlon_projected[ ,2]

  return(c(x, y))
}
