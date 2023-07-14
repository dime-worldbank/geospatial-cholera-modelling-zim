
# This script is intended to demonstrate the statistical approach taken
# in Ayling et al. 2022 "A Stitch in Time: The importance of WSS Infrastructure Maintenance on Cholera Risk. 
# A Geospatial Analysis in Harare, Zimbabwe"


# Get worker function
source("R/validate_fit.R")
source("R/rate_density.R")

# Load population raster to use as offset
exposure <- raster("data/population.tif")

# Load covariates
covariates <- stack("data/covariates.tif")
names(covariates) <- c("Sanitation_risk", "Proportion_piped_water", 
                       "Proportion_unimproved_water", "Poverty")
plot(covariates)

# Load dummy case data
cases <- read.csv("data/dummy_cases.csv")
cases_sf <- st_as_sf(cases, coords = c("lng", "lat"))

# install.packages('leaflet')
library(leaflet)
# Visualize cases by timeslice
leaflet() %>% addProviderTiles("CartoDB.Positron") %>% 
  addCircleMarkers(data=cases_sf, color = rainbow(5)[cases_sf$timeslice]) %>%
  addLegend(colors = rainbow(5), labels = 1:5, title = "Week")

# Fit model
# Run function
res_lagged <- rate_density(cases,
                           space_time = TRUE,
                           k = 150,
                           debug = T,
                           prediction_periods = 'each',
                           approx_num_int_points = 1000,
                           covariates = covariates,
                           #covariates_dynamic = list(imported_incidence = imported_incidence),
                           smooth_effects = FALSE,
                           discrete = TRUE,
                           exposure = exposure)

# Get model fit summary
summary(res_lagged$model)

# Plot predicted risk ()
plot(res_lagged$pred_raster, zlim=c(0,0.0432),
     main = paste("timeslice", 1:5))

# Validate model fit per timeslice
validation <- validate_fit(res_lagged$pred_raster,
                           cases_sf,
                           exposure)
validation$obs_pred_per_time_plot
