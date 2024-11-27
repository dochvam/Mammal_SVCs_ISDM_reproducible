###############################################################################
#'
#' setup_code/define_fundamentals
#' Author: Ben R. Goldstein
#' Description: define covariates, target species, etc. for non-lineage models
#'
###############################################################################

library(tidyverse)

occ_covars <- c("Pop_den_sqrt_scaled", "Temp_max_scaled", 
                "Precipitation_scaled", "EVI_mean_scaled",
                "Terrain_roughness_scaled", "Agriculture_scaled",
                "Grassland_scaled")

det_covars <- c("yday_scaled", "yday_scaled_sq", "Canopy_height_scaled",
                "log_roaddist_scaled")

covars_toscale <- c("Arable", "Aridity", "Forest", "Grassland",
                    "EVI_mean", "EVI_variability",
                    "Pastureland", "Pop_den_sqrt", "Precipitation", "Agriculture",
                    "Shrubland", "Temp_max", "Terrain_roughness", "Wetlands")


covar_type_df <- data.frame(
  covariate = gsub("_scaled", "", occ_covars),
  type = c("Human (n=2)", "Climate (n=2)", "Climate (n=2)",
           "Landscape (n=3)", "Landscape (n=3)",
           "Human (n=2)", "Landscape (n=3)")
)



var_colors <- c(
  "#a6cee3",
  "#e41a1c",
  "#377eb8",
  "#4daf4a",
  "#ff7f00",
  "#a65628",
  "#f781bf",
  "#ffff33",
  "#984ea3",
  "#dddddd",
  "#999999"
)
names(var_colors) <- c(gsub("_scaled", "", occ_covars), "Out of range",
                       "Competitive")



order_info <- read_csv("intermediate/north_american_terr_mammals.csv") %>% 
  select(`Sci name`, Order) %>% 
  filter(!is.na(Order))


taxon_key <- read_csv("intermediate/target_species.csv") %>% 
  mutate(species = tolower(common_name)) %>% 
  select(-common_name) %>% 
  mutate(sci_name_clean = paste0(taxon_genus, "_", taxon_spec),
         common_name_clean = gsub("[ -']", "_", species)) %>% 
  filter(!is.na(taxon_genus)) %>% 
  left_join(order_info, by = c(sci_name = "Sci name")) %>% 
  mutate(nobs_in_range = NA) %>% 
  filter(sci_name_clean != "Sus_scrofa")
