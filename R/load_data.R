#' Load and crop data from data file
#'
#' @param filename string relative path to data file
#' @param limits numeric vector extent of region to be cropped (xmin, xmax, ymin, ymax)
#' @return data SpatRaster layer of cropped data
#'
load_and_crop_data <- function(filename, subdataset, limits = NULL){
  message('Started: loading/cropping data from file: ', filename)

  data <- load_data(filename, subdataset)

  if (is.vector(extent, mode = 'numeric')){
    crop_data(data, limits)
  }

  message('Started: loading/cropping data from file: ', filename)

  return(data)
}

#'
load_data <- function(filename, subdataset){
  message('Started: loading data from file: ', filename)

  data <- rast(filename, subds = subdataset)

  message('Stopped: loading data from file: ', filename)

  return(data)
}

#'
crop_data <- function(data, limits){
  message('Started: cropping data within limits: ', list(limits))

  crs_data <- crs(data, proj = TRUE)
  limits_projected <- project_degrees(limits[1:2],limits[3:4], crs_data)

  data_cropped <- crop(data, ext(limits_projected))

  data_cropped_projected <- project(data_cropped,"EPSG:4326")

  message('Stopped: cropping data within limits: ', list(limits))

  return(data_cropped_projected)
}

