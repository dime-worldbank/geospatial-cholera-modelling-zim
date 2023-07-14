# Replication Package for: [Geospatial Analysis of Cholera Risk and Water and Sanitation Infrastructure in Harare, Zimbabwe](https://documents.worldbank.org/en/publication/documents-reports/documentdetail/099631202282320091/idu1d0474b1510f5414def1b645112a4f569be12)

## Code
Code for geospatial modeling of cholera risk from Harare, Zimbabwe in 2018-19. 

#### Main Script
* `fit_model.R`: Main script that runs all code

#### Organization
Code is organized into the folder `/R` as follows:

* `get_ppm.R`
* `rate_density.R`
* `space_time_ppmify.R`
* `space_time_ppmify_helpers.R`: 
* `validate_fit.R`

**Note:** The package `velox` in script `space_time_ppmify_helpers.R` needs to be updated to a new package as this previous one was removed from the CRAN repository. Please find the complete note in the following link https://CRAN.R-project.org/package=velox

## Data

### Organization 
Data can be found in the `/data` folder

## To replicate analysis
1. Clone/download this github repository
2. Run `fit_model.R`; this runs all scripts needed to replicate the analysis.
