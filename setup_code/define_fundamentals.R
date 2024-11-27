###############################################################################
#'
#' setup_code/define_fundamentals
#' Author: Ben R. Goldstein
#' Description: define covariates, target species, etc. for lineage models
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

lineage_key <- data.frame(
  species = c(
              "nine-banded armadillo",
              "american black bear",
              "bobcat", 
              "california ground squirrel",
              "douglas's squirrel",
              "eastern chipmunk",
              "eastern gray squirrel",
              "mule deer", 
              "grey fox",
              "red fox",
              "red squirrel",
              "snowshoe hare",
              "mexican flying squirrel",
              "striped skunk",
              "white-tailed antelope squirrel",
              "white-tailed deer"
              ),
  lineage_file = c(
    "lineages_final/Armadillo new phylo/Armadillo new phylo.shp",
    "lineages_final/Black bear new phylo/Black bear new phylo.shp",
    "lineages_final/Bobcat_new_phylo/Bobcat new phylo/Bobcat new phylo.shp",
    "lineages_final/California ground squirrel phylo partial/California ground squirrel phylo partial/California ground squirrel phylo partial.shp",
    "lineages_final/Douglas_squirrel_new_phylo/Douglas squirrel new phylo/Douglas squirrel new phylo.shp",
    "lineages_final/Eastern_chipmunk_phylo_partial/Eastern chipmunk phylo partial/Eastern chipmunk phylo partial.shp",
    "lineages_final/Eastern_gray_squirrel_new_phylo/Eastern gray squirrel new phylo/Eastern gray squirrel new phylo.shp",
    "lineages_final/Mule_deer_new_phylo/Mule deer new phylo/Mule deer new phylo.shp",
    "lineages_final/Northern_gray_fox_phylo_partial/Northern gray fox phylo partial/Northern gray fox phylo partial.shp",
    "lineages_final/Red_fox_new_phylo/Red fox new phylo/Red fox new phylo.shp",
    "lineages_final/Red_squirrel_lineages_iSDM/Red_squirrel_phylo_partial.shp",
    "lineages_final/Snowshoe_hare_new_phylo/Snowshoe hare new phylo/Snowshoe hare new phylo.shp",
    "lineages_final/Southern flying squirrel phylo partial/Southern flying squirrel phylo partial/Southern flying squirrel phylo partial.shp",
    "lineages_final/Striped_skunk_new_phylo/Striped skunk new phylo/Striped skunk new phylo.shp",
    "lineages_final/White-tailed antelope squirrel phylo partial/White-tailed antelope squirrel phylo partial/White-tailed antelope squirrel phylo partial.shp",
    "lineages_final/White-tailed_deer_new_phylo/White-tailed deer new phylo/White-tailed deer new phylo.shp"
  )
) %>% 
  mutate(common_name_clean = gsub("[ -']", "_", species))



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
  filter(!is.na(Order)) %>% 
  filter(species %in% lineage_key$species) %>% 
  filter(species != "white-tailed deer", species != "nine-banded armadillo") %>% 
  mutate(nobs_in_range = NA)
