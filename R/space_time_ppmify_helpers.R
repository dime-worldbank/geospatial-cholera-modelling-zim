#' The get_int_points_exposure_weights function
#'
#' Helper function for space_time_ppmify
#' @param ppmx 
#' @param ppm_cases_points_counts
#' @param exposure_raster
#' @param num_periods
#' @import velox

# install.packages('remotes')
# remotes::install_version('velox', version = '0.1.0')
#library(velox)
library(terra)
get_int_points_exposure_weights <- function(ppmx, ppm_cases_points_counts, exposure_raster, periods){

  # First identify which are the local_cases and integration rows
  # in the ppm object
  ppm_int_points <- ppmx[ppmx$points==0,]

  # extract population within each voronoi polygon around each integration point
  # First remove any pixels with cases as you need to estimate population in pixels without cases
  dd <- deldir::deldir(ppm_int_points$x, ppm_int_points$y)
  tiles <- deldir::tile.list(dd)
  
  polys <- vector(mode = "list", length = length(tiles))
  for (i in seq(along = polys)) {
    pcrds <- cbind(tiles[[i]]$x, tiles[[i]]$y)
    pcrds <- rbind(pcrds, pcrds[1,])
    polys[[i]] <- sp::Polygons(list(sp::Polygon(pcrds)), ID = as.character(i))
  }
  spoly <- sp::SpatialPolygons(polys)
  
  ppm_int_points_period <- NULL
  for(j in periods){
    exposure_raster_non_case_pixels <- exposure_raster
    
    # Remove any 'case' pixels
    ppm_case_points_coords <- ppm_cases_points_counts[ppm_cases_points_counts$period == j, c("x", "y")]
    exposure_raster_non_case_pixels[raster::cellFromXY(exposure_raster_non_case_pixels,
                                                    ppm_case_points_coords)] <- NA

  # Extract from offset raster
  ppm_int_points$exposure <- terra::extract(exposure_raster_non_case_pixels,
                                            spoly,
                                            fun = sum,
                                            na.rm = TRUE)
  
  # bind with previous periods
    ppm_int_points$period <- j
    ppm_int_points_period <- rbind(ppm_int_points_period, ppm_int_points)
  }
  
  # Remove any points with 0 offset
  if(length(which(ppm_int_points_period$exposure <= 0))>0){
  ppm_int_points_period <- ppm_int_points_period[-which(ppm_int_points_period$exposure <= 0),]
  }
  
  return(ppm_int_points_period)
}  


#' The aggregate_points_space_time function
#'
#' Helper function for space_time_ppmify
#' @param points
#' @param ppmx
#' @param num_periods
#' @param date_start_end
#' @param reference_raster
#' @import lubridate
library(lubridate)

# Deal with case points
# First aggregate any cases occuring in the same pixel in the same month
# # loop through each month and identify any cases in the same pixel
aggregate_points_space_time <- function(points, ppmx, reference_raster, num_periods){

    # if(is.null(points$period)){
    #   num_periods <- 1
    # }else{
    #   num_periods <- length(unique(points$period)) 
    # }
      ppm_cases_points <- ppmx[ppmx$points==1,]
      ppm_cases_points_counts <- ppm_cases_points[FALSE,]
      
      # Define date breaks
      #dates <- seq(ymd(date_start_end[1]), ymd(date_start_end[2]), 1)
      # if(num_periods > 1){
      # date_breaks <- c(levels(cut.Date(dates, periods, right=TRUE, include.lowest = TRUE)),
      #                  date_start_end[2])
      # }else{
      #   date_breaks <- date_start_end
      # }

      for(i in unique(points$period)){
        
        if(num_periods == 1){
          cases_period <- ppm_cases_points
        }else{
          #cases_model_period <- as.numeric(cut.Date(ymd(points$date), ymd(periods), include.lowest = TRUE))
          cases_period <- ppm_cases_points[points$period==i,]
        }
        
        if(nrow(cases_period)>0){
          
          cases_period$period <- i
          
          # calculate which cell each case is in
          cases_period$case_pixel <- raster::cellFromXY(reference_raster,
                                                       cbind(cases_period$x, cases_period$y))
          
          # Aggregate number of cases by cell
          cell_counts <- raster::aggregate(cases_period$case_pixel,
                                           by = list(cases_period$case_pixel), length)
          
          # create aggregated ppm data frame for aggregated case counts
          cases_period_trimmed <- cases_period[match(cell_counts$Group.1, cases_period$case_pixel),]
          
          # Aggregate case coordinates to be the centroid of the pixel
          cases_period_trimmed[,c("x","y")] <- sp::coordinates(reference_raster)[cases_period_trimmed$case_pixel,]
          
          # change number of points to be aggregate number
          cases_period_trimmed$points <- cell_counts$x
          
          # take the columns you need
          cases_period_trimmed<- cases_period_trimmed[, 1:5]
          #names(cases_period_trimmed) <- names(ppm_cases_points)
          
          # bind to complete the loop
          ppm_cases_points_counts <- rbind(ppm_cases_points_counts, cases_period_trimmed)
        }
      }
      return(ppm_cases_points_counts)
}



#' extract_dynamic_covariates
#'
#' Helper function for space_time_ppmify
#' @param covariates_dynamic 
#' @param reference_raster
#' @import raster

extract_dynamic_covariates <- function(covariates_dynamic, ppm_df, reference_raster){

  n_layers_per_stack <- as.vector(sapply(covariates_dynamic, nlayers))
  if(!all(n_layers_per_stack %in% length(unique(ppm_df$period)))){
    stop(paste('The number of layers in `covariates_dynamic` does not equal the number of periods'))
  }
  
  for(n_covar in 1:length(covariates_dynamic)){
    
    covariates_dynamic_n <- covariates_dynamic[[n_covar]]

  if(!((extent(reference_raster)==extent(covariates_dynamic_n)) &
       (res(reference_raster)[1]==res(covariates_dynamic_n)[1]))){
      covariates_dynamic_n <- resample(covariates_dynamic_n, reference_raster)
  }
      
      ppm_df_loop <- NULL
      periods <- sort(unique(ppm_df$period ))
      for(layer in 1:nlayers(covariates_dynamic_n)){
        
        ppm_df_sub <- subset(ppm_df, ppm_df$period == periods[layer])
        ppm_df_sub[[ names(covariates_dynamic)[n_covar] ]] <- raster::extract(covariates_dynamic_n[[layer]],
                                                           ppm_df_sub[,c("x", "y")])
        ppm_df_loop <- rbind(ppm_df_loop, ppm_df_sub)
        
      }
      ppm_df <- ppm_df_loop
    }
  return(ppm_df)
}  

