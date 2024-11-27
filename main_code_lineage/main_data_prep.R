###############################################################################
#'
#' main_code_lineage/main_data_prep.R
#' Author: Ben R. Goldstein
#' Description: Loads in all the data needed for running the model, including
#'    iNat, camera traps, covariates, spatial info, etc.
#'
###############################################################################


library(tidyverse)
library(terra)
library(tidyterra)
library(unmarked)
source("setup_code/main_helper_fn.R")
source("setup_code/define_fundamentals.R")

#### Prep: read in data, establish some key parameters ####
  continent_crs <- '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83'

#### Create a spatial grid covering North America, the "continental grid" ####
  continental_grid_scale2 <- rast("intermediate/lineages/continental_grid_scale2.grd")
  continental_grid_scale3 <- rast("intermediate/lineages/continental_grid_scale3.grd")
  continental_grid_scale4 <- rast("intermediate/lineages/continental_grid_scale4.grd")
  grid_translator <- read_csv("intermediate/lineages/continental_grid_translator.csv")
  cells_to_states <- read_csv("intermediate/lineages/cells_to_states.csv")



# adjacency information for the continental_grid
  continent_adj_info <- readRDS("intermediate/continentAdjInfo.RDS")
  grid_translator_wspec <- read_csv("intermediate/grid_translator_wspecs_lineage.csv")
  adj_info_list <- readRDS("intermediate/lineages/adj_info_byspec.RDS")


#### Map species lineages to the spatial grid ####
  lineage_df <- read_csv("intermediate/lineages/species_lineages.csv")

#### Map ecoregions to the spatial grid ####
  ecoregion_df <- read_csv("intermediate/lineages/ecoregions.csv")

#### Calculate iNaturalist effort, reporting rates in each cell ####
  grid_counts <- read_csv("intermediate/inat_cts_lineage.csv")

#### Retrieve spatial covariates in each cell ####
  covar_values <- read_csv("intermediate/gridcell_covars_lineage.csv")

#### Retrieve camera trap data####
  spec_dat_list <- readRDS("intermediate/ct_datalist_lineage.RDS")
