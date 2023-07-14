# install.packages('hydroGF')
library(hydroGOF)
validate_fit <- function(pred_inc_stack,
                         observations_sf,
                         exposure){
  
  # Check total predicted numbers per timeslice
  pred_per_timeslice <- floor(cellStats((pred_inc_stack/1000) *  exposure, sum))
  obs_per_timeslice <- table(observations_sf$timeslice)
  obs_pred_per_time_plot <- ggplot() +
    geom_line(aes(c(0,max(pred_per_timeslice, obs_per_timeslice)), 
                  c(0,max(pred_per_timeslice, obs_per_timeslice))), color = "blue") +
    geom_point(aes(pred_per_timeslice, obs_per_timeslice),
               size=2) +
    xlab("Predicted number of cases per time period") +
    ylab("Observed number of cases per time period")
  
  # Check r2 between observed and predicted at pixel level (per timeslice)
  observations_stack <- stack()
  for(time in sort(unique(observations_sf$timeslice))){
    cases_raster <- rasterize(as(subset(observations_sf, timeslice==time), "Spatial"), 
                              exposure, 
                              fun='count',
                              field = "timeslice")
    observations_stack <- stack(observations_stack, 
                                cases_raster)
  }
  
  pred_cases_per_pixel <-(pred_inc_stack/1000) *  exposure

  obs_pred_per_time_pixel_plot <- ggplot() +
    geom_line(aes(c(0,max(c(pred_cases_per_pixel[], observations_stack[]),na.rm=T)), 
                  c(0,max(c(pred_cases_per_pixel[], observations_stack[]),na.rm=T))), color = "blue") +
    geom_point(aes(pred_cases_per_pixel[], observations_stack[]),
               size=2, color = rgb(1,0,0.5, 0.2)) +
    xlab("Predicted number of cases per time period") +
    ylab("Observed number of cases per time period")

  RSME <- rmse(as.vector(pred_cases_per_pixel[]), as.vector(observations_stack[]),
               na.rm=T)

  # Return
  return(list(pred_per_timeslice = pred_per_timeslice, 
              obs_per_timeslice = obs_per_timeslice,
              obs_pred_per_time_plot = obs_pred_per_time_plot,
              obs_pred_per_time_pixel_plot = obs_pred_per_time_pixel_plot,
              RSME = RSME))
  
}