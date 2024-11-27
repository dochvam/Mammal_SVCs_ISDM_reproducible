###############################################################################
#'
#' setup_code/main_helper_fn.R
#' Author: Ben R. Goldstein
#' Description: a couple misc. functions and objects
#'
###############################################################################

library(lubridate)
library(progress)
library(unmarked)
library(tidyverse)
library(terra)


covars_toscale <- c("Temp_max", "Precipitation", "Pop_den", "Terrain_roughness",
                    "Aridity", "Forest", "Agriculture", "Canopy_height",
                    "Grassland", "year", "longitude", "latitude",
                    "log_roaddist")

deployment_cols_needed <- c("project_deployment_id", "feature_water_source", 
                            "feature_game_trail", "feature_human_trail", 
                            "feature_road", "log_roaddist", "Canopy_height",
                            "year", "longitude", "latitude")

manual_range_bbox <- data.frame(
  species = c(NA),
  longitude_min = c(-95),
  longitude_max = c(Inf),
  latitude_min = c(-Inf),
  latitude_max = c(Inf)
)

make_name_clean <- function(x) {
  gsub(" ", "_", x)
}
