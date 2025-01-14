###############################################################################
#'
#' main_code_nolineage/model_fn.R
#' Author: Ben R. Goldstein
#' Description: Helper functions containing the core model estimation
#'    code 
#'
###############################################################################

check_collinearity <- function(mtx, threshold = 0.5) {
  temp <- cov(mtx)
  diag(temp) <- 0
  which(temp > threshold, arr.ind = TRUE)
}

make_ecoregion_interaction_dat <- function(covs_df, occ_covars,
                                           ecoregion_final,
                                           ecoregion_levels) {
  
  
  covs_df <- left_join(covs_df, ecoregion_final, by = c("spec_range_gridID"))
  
  ecoregion_vector <- factor(covs_df$ecoregion, levels = ecoregion_levels)
  contrasts(ecoregion_vector) <- contr.sum
  # Levels are first thru. (n-1)th. -1s are the final ecoregion
  model_mtx <- as.matrix(model.matrix(~ecoregion_vector)[, -1])
  
  ecoregion_int_groups <- rep(1:length(occ_covars), each = ncol(model_mtx))
  ecoregion_interaction_data <- matrix(ncol = length(ecoregion_int_groups), 
                                     nrow = nrow(covs_df))
  
  ecoregion_int_names <- paste0(
    expand.grid(occ_covars, 1:(ncol(model_mtx)))[,1], "_L",
    expand.grid(occ_covars, 1:(ncol(model_mtx)))[,2]
  )
  
  for (i in 1:length(occ_covars)) {
    inds <- which(ecoregion_int_groups == i)
    
    this_cov_mtx <- matrix(rep(unlist(covs_df[, occ_covars[i]]),
                               ncol(model_mtx)),
                           nrow = nrow(model_mtx), ncol = ncol(model_mtx))
    
    ecoregion_interaction_data[, inds] <- model_mtx * this_cov_mtx
  }
  return(list(
    df = ecoregion_interaction_data,
    ecoregion_vector = ecoregion_vector,
    ecoregion_int_groups = ecoregion_int_groups
  ))
}


# longform_matrix_to_array <- function(mtx, y.nrow, y.ncol) {
#   
#   rtn = array(NA, dim = c(y.nrow, y.ncol, ncol(detcovs)))
#   
#   for (i in 1:ncol(mtx)) {
#     rtn[,,i] <- matrix(as.numeric(mtx[,i]), nrow = y.nrow, byrow = T)
#   }
#   
#   dimnames(rtn)[[3]] <- colnames(mtx)
#   rtn
# }

# matrix_to_longform <- function(mtx) {
#   as.numeric(t(mtx))
# }

get_inits <- function(modtype, occ_covars, det_covars, spatial_model,
                      this_adj_info = NULL, nEcoregion, 
                      ecoregion_datlist, ndepl) {
  
  if (spatial_model == "SVCs") {
    inits_list <- list(
      # Always
      lambda_intercept = 0.2,
      p_intercept = 0.5,
      lambda_beta_mean = rnorm(length(occ_covars), 0, sd = 0.1),
      p_gamma = rnorm(length(det_covars), 0, sd = 0.1),
      overdisp_inat = 0.1,
      # overdisp_museum = 0.1,
      dummy = 1,
      p_alpha = rnorm(length(occ_covars), 0, sd = 0.1),
      theta0 = -5,
      theta1 = 1,
      # theta0_museum = -5,
      # theta1_museum = 1,
      spat_magnitude_intercept = 1,
      spat_ranef = rnorm(length(this_adj_info$num), 0, 1),
      
      eco_effect =  matrix(rnorm(nEcoregion*length(occ_covars), sd = 0.5),
                           ncol = length(occ_covars)),
      CAR_effect = matrix(
        rnorm(length(this_adj_info$num) * length(occ_covars), 0, 1),
        ncol = length(occ_covars)
        ),
      
      CAR_gamma = rbinom(length(occ_covars), size = 1, prob = 1),
      eco_gamma = rbinom(length(occ_covars), size = 1, prob = 1),
      
      sigma_CAR_0 = rep(0.25, length(occ_covars)),
      sigma_eco_0 = rep(0.25, length(occ_covars)),
      
      sigma_CAR = rep(0.25, length(occ_covars)),
      sigma_eco = rep(0.25, length(occ_covars)),
      
      det_ranef_sd = 0.1,
      det_ranef = runif(ndepl, -0.3, 0.3),

      p_incl_CAR = 0.5,
      p_incl_eco = 0.5,
      logDens = 0
    )
    return(inits_list)
  }
  
  inits_list <- list(
    # Always
    lambda_intercept = 0.2,
    p_intercept = 0.5,
    lambda_beta = rnorm(length(occ_covars), 0, sd = 0.1),
    p_gamma = rnorm(length(det_covars), 0, sd = 0.1),
    overdisp_inat = 0.1,
    # overdisp_museum = 0.1,
    dummy = 1,
    p_alpha = rnorm(length(occ_covars), 0, sd = 0.1)
  )
  

  if (spatial_model == "CAR") {
    inits_list[["spat_magnitude"]] <- 1
    inits_list[["spat_ranef"]] <- rnorm(length(this_adj_info$num), 0, 0.1)
  }
  
  if (modtype == "joint") {
    inits_list[["theta0"]] <- -5
    inits_list[["theta1"]] <- 1
    # inits_list[["theta0_museum"]] <- -5
    # inits_list[["theta1_museum"]] <- 1
  } else if (modtype == "separate") {
    inits_list[["p_alpha"]] <- rnorm(length(occ_covars), 0, sd = 0.1)
    inits_list[["inat_beta"]] <- rnorm(length(occ_covars), 0, sd = 0.1)
    inits_list[["log_inat_intercept"]] <- -5
    # inits_list[["museum_beta"]] <- rnorm(length(occ_covars), 0, sd = 0.1)
    # inits_list[["log_museum_intercept"]] <- -5
    
    if (spatial_model == "CAR") {
      # inits_list[["museum_spat_magnitude"]] <- 0.2
      # inits_list[["museum_spat_ranef"]] <- rnorm(length(this_adj_info$num), 0, 0.1)
      inits_list[["inat_spat_magnitude"]] <- 0.2
      inits_list[["inat_spat_ranef"]] <- rnorm(length(this_adj_info$num), 0, 0.1)
    }
  }
  
  return(inits_list)
}

get_monitors <- function(modtype, inat_disttype, spatial_model) {
  
  if (spatial_model == "SVCs") {
    monitors <- c(
      "lambda_intercept",
      "p_intercept",
      "lambda_beta_mean",
      "p_gamma",
      "p_alpha",
      "brier_score",
      "spat_magnitude_intercept",
      "overdisp_inat",
      
      "det_ranef_sd",

      "sigma_CAR",
      "sigma_eco",
      
      "CAR_gamma",
      "eco_gamma",
      
      "p_incl_CAR", 
      "p_incl_eco",
      
      "theta0", "theta1"#,
      # "theta0_museum", "theta1_museum"
    )
    
    monitors2 <- c(
      "spat_ranef",
      "lambda_beta",
      "CAR_effect",
      "eco_effect"
    )
    
    return(list(
      monitors = monitors, monitors2 = monitors2
    ))
  }
  
  monitors <- c(
    "lambda_intercept",
    "p_intercept",
    "lambda_beta",
    "p_gamma",
    "brier_score"
  )
  
  if (spatial_model == "CAR") {
    monitors <- c(monitors, "spat_magnitude")
  }
  
  if (inat_disttype == "NB") {
    monitors <- c(monitors, "overdisp_inat"#, "overdisp_museum"
                  )
  }
  
  if (modtype == "joint") {
    monitors <- c(monitors, "theta0", "theta1", "p_alpha"#,
                  #"theta0_museum", "theta1_museum"
                  )
  } else if (modtype == "separate") {
    monitors <- c(monitors, "p_alpha", "inat_beta", "log_inat_intercept"#,
                  #"museum_beta", "log_museum_intercept"
                  )
    
    if (spatial_model == "CAR") {
      monitors <- c(monitors, "inat_spat_magnitude"#, "museum_spat_magnitude"
                    )
    }
  }
  
  if (spatial_model == "CAR") {
    monitors2 <- "spat_ranef"
    if (modtype == "separate") {
      monitors2 <- c(monitors2, "inat_spat_ranef"#, "museum_spat_ranef"
                     )
    }
  } else {
    monitors2 <- "dummy"
  }
  
  return(list(
    monitors = monitors, monitors2 = monitors2
  ))
}



#' @param nfolds Number of folds to use for k-fold cross-validation
#' @param folds_by Method to use for selecting holdout data.
#'                 Either "camera" or "subproject"
#' @param modtype 
#' @param seed
#' 
fit_integrated_model_with_CV <- function(nfolds, 
                                         kfold_mode = "kfold",
                                         folds_by = "camera", 
                                         modtype,
                                         inat_disttype = "Pois",
                                         spatial_model = "None",
                                         species,
                                         ni, 
                                         nb, 
                                         nc, 
                                         nt, 
                                         nt2,
                                         seed,
                                         subset_CT = 1,
                                         subset_inat = 1,
                                         suffix = "",
                                         debug_transl = FALSE,
                                         overwrite = FALSE,
                                         return_data_summary = FALSE) {
    
    outfile <- paste0("intermediate/integration_results/final_", 
                      species, "_", modtype, "_", spatial_model, suffix, "_noLineage_samples.RDS")
    if (!overwrite && file.exists(outfile)) {
      write(paste0("Not running ", species, " suffix ", suffix,
                        " at time ", Sys.time()), file = "log.txt", append = TRUE)
      return()
    }
  start_time <- Sys.time()
  # I can disable k-fold validation by setting nfolds to 0. I still construct
  # a holdout dataset so that I don't have to modify the model code, but the
  # holdout data are not actually held out (instead, they contribute to the 
  # model likelihood). In this case, the Brier score is less useful
  if (nfolds == 0) {
    holdout_data <- FALSE
    nfolds <- 5
    kfold_mode <- "single"
    folds_by <- "camera"
  } else {
    holdout_data <- TRUE
  }
  
  # Check that args are valid
  stopifnot(folds_by %in% c("camera", "subproject"))
  stopifnot(kfold_mode %in% c("kfold", "single"))
  stopifnot(modtype %in% c("joint", "separate"))
  stopifnot(inat_disttype %in% c("Pois", "NB"))
  stopifnot(spatial_model %in% c("None", "CAR", "SVCs"))
  
  source("main_code_nolineage/model_code_main.R")
  
  if (spatial_model == "CAR") {
    model_code <- joint_model_CAR
  } else if (spatial_model == "SVCs") {
    model_code <- joint_model_SVCs_asRanefs
  } else {
    model_code <- joint_model_nonspatial
  }
  
  ### Load in and format the data
  if (!return_data_summary) source("main_code_nolineage/main_data_prep.R")
  
  set.seed(seed)
  
  # Check that species / modtype is valid
  stopifnot(species %in% taxon_key$common_name_clean)

  ###### Preparatory data handling for this species ######
  # Rename the species-specific grid ID column to match. This is the ID that 
  #   corresponds to the species' adjacency matrix, i.e. the "spatcell"
  colnames(grid_translator_wspec)[
    grepl(species, colnames(grid_translator_wspec))] <- "spec_range_gridID"
  
  # Get the index for the species in the list of species
  spec_index <- which(taxon_key$common_name_clean == species)
  
  # Get the species detection history from the umf list
  this_CT_umf <- spec_dat_list[[spec_index]]$umf
  
  ###### Process CT data ######
  y.all        <- this_CT_umf@y %>% 
    as.data.frame() %>% 
    mutate(site_ID = row_number()) %>% 
    pivot_longer(cols = paste0(1:ncol(this_CT_umf@y)))
  
  inds_to_keep <- which(!is.na(y.all$value))
  y.all <- y.all[inds_to_keep,]
  
  # Get a df of every scale2 cell in the species' range
  cells_inRange <- grid_translator_wspec %>% 
    filter(!is.na(spec_range_gridID))
  
  ###### Process site-level covariate data ###### 
  
  spec_dat_list[[spec_index]]$coords <- spec_dat_list[[spec_index]]$coords %>% 
    filter(!is.na(scale2_grid_ID))
  
  occ.covs.all <- left_join(spec_dat_list[[spec_index]]$coords, covar_values,
                            by = c("scale2_grid_ID" = "grid_cell")) %>% 
    mutate(site_ID = row_number()) %>% 
    left_join(grid_translator_wspec)
  occ.covs.all$Intercept <- 1
  occ.covs.all$Pop_den_sqrt <- sqrt(occ.covs.all$Pop_den)
  occ.covs.all$Agriculture <- occ.covs.all$Pastureland + occ.covs.all$Arable
  occ.covs.all$subproject_ID <- spec_dat_list[[spec_index]]$umf@siteCovs$subproject_ID
  
  det.covs.all <- this_CT_umf@obsCovs
  det.covs.all <- det.covs.all[inds_to_keep,]
  det.covs.all$Intercept <- 1
  det.covs.all$site_ID <- y.all$site_ID
  
  stopifnot(nrow(det.covs.all) == nrow(y.all))
  stopifnot(nrow(occ.covs.all) == max(y.all$site_ID))
  
  ###### Drop problem sites ######
  
  # drop sites with only one obs. (for NIMBLE, for now)
  sites_to_drop <- count(y.all, site_ID) %>% filter(n == 1)
  
  # browser()
  # Drop sites with missing occupancy covariates
  if (any(rowSums(is.na(occ.covs.all[, gsub("_scaled", "", occ_covars)])) > 0)) {
    sites_to_drop <- sites_to_drop %>% 
      bind_rows(data.frame(
        site_ID = occ.covs.all$site_ID[
          rowSums(is.na(occ.covs.all[, gsub("_scaled", "", occ_covars)])) > 0
        ],
        n = 1
      ))
  }
  

  # Drop sites by subsetting, if subset < 1
  if (subset_CT < 1) {
    testSubset <- sample_frac(occ.covs.all, size = 1-subset_CT)
    
    sites_to_drop <- 
      bind_rows(sites_to_drop, data.frame(site_ID = testSubset$site_ID, n = 1)) %>% 
      distinct()
  }
  
  
  # Execute dropping in occ.covs, det.covs, and y if there are any to drop
  if (nrow(sites_to_drop) > 0) {
    
    new_site_ID_df <- data.frame(
      old_site_ID = occ.covs.all$site_ID,
      included = !(occ.covs.all$site_ID %in% sites_to_drop$site_ID),
      new_site_ID = NA
    )
    new_site_ID_df$new_site_ID[new_site_ID_df$included] <-
      1:sum(new_site_ID_df$included)
    
    occ.covs.all <- occ.covs.all %>% 
      left_join(new_site_ID_df, by = c("site_ID" = "old_site_ID")) %>% 
      filter(!is.na(new_site_ID)) %>% 
      mutate(site_ID = new_site_ID)
    
    y.all <- y.all %>% 
      left_join(new_site_ID_df, by = c("site_ID" = "old_site_ID")) %>% 
      filter(!is.na(new_site_ID)) %>% 
      mutate(site_ID = new_site_ID)
    
    det.covs.all <- det.covs.all %>% 
      left_join(new_site_ID_df, by = c("site_ID" = "old_site_ID")) %>% 
      filter(!is.na(new_site_ID)) %>% 
      mutate(site_ID = new_site_ID)
    
  }
  stopifnot(nrow(det.covs.all) == nrow(y.all))
  stopifnot(nrow(occ.covs.all) == max(y.all$site_ID))
  
  ###### Scale occupancy covariates ######
  scaling_factors <- data.frame(
    covar = covars_toscale, mean = NA, sd = NA
  )
  # 
  # browser()
  covar_values_wdat <- covar_values %>% 
    filter(grid_cell %in% cells_inRange$scale2_grid_ID) %>% 
    na.omit()
  for (i in 1:length(covars_toscale)) {
    scaling_factors$mean[i] <- mean(unlist(covar_values_wdat[, covars_toscale[i]]))
    scaling_factors$sd[i]   <- sd(unlist(covar_values_wdat[, covars_toscale[i]]))
    
    occ.covs.all[[paste0(covars_toscale[i], "_scaled")]] <- as.numeric(
      (occ.covs.all[, covars_toscale[i]] - scaling_factors$mean[i]) / 
        scaling_factors$sd[i]
    )
  }
  
  ###### Scale detection covariates ######
  scaling_factors_det <- data.frame(
    covar = det_covars, mean = NA, sd = NA
  )
  for (i in 1:length(det_covars)) {
    scaling_factors_det$mean[i] <- mean(det.covs.all[, det_covars[i]])
    scaling_factors_det$sd[i]   <- sd(det.covs.all[, det_covars[i]])
    det.covs.all[[det_covars[i]]] <- as.numeric(
      (det.covs.all[, det_covars[i]] - scaling_factors_det$mean[i]) / 
        scaling_factors_det$sd[i]
    )
  }
  
  
  ###### Set up the holdout data ###### 
  if (folds_by == "camera") {
    fold_vec <- sample(1:nfolds, replace = TRUE, size = max(y.all$site_ID))
    
    folds_df <- data.frame(
      camera_ID = sort(unique(y.all$site_ID)),
      fold = fold_vec[sort(unique(y.all$site_ID))]
    )
  } else if (folds_by == "subproject") {
    fold_vec <- sample(1:nfolds, replace = TRUE, size = max(occ.covs.all$subproject_ID))
    
    y.all.temp <- left_join(y.all, occ.covs.all[, c("site_ID", "subproject_ID")],
                            by = "site_ID")
    
    folds_df <- distinct(y.all.temp, site_ID, subproject_ID) %>% 
      rename(camera_ID = site_ID) %>% 
      mutate(
        fold = fold_vec[subproject_ID]
      )
  } else {
    stop('Invalid "folds_by" value')
  }
  
  
  ###### Prepare iNaturalist data and covariate data ######

  covar_values$Agriculture <- covar_values$Arable + covar_values$Pastureland
  covar_values$Pop_den_sqrt <- sqrt(covar_values$Pop_den)
  
  xdat_allS2 <- covar_values %>% 
    rename(scale2_grid_ID = grid_cell) %>% 
    filter(scale2_grid_ID %in% cells_inRange$scale2_grid_ID) %>% 
    left_join(cells_inRange[, c("scale2_grid_ID", "scale3_grid_ID")])
  
  inat_dat <- grid_counts %>% 
    rename(scale3_grid_ID = grid_cell) %>% 
    filter(scale3_grid_ID %in% xdat_allS2$scale3_grid_ID) %>% 
    arrange(scale3_grid_ID) %>% 
    mutate(inat_datcell_ID = as.numeric(as.factor(scale3_grid_ID))) %>% 
    select(scale3_grid_ID, effort = n, inat_datcell_ID, all_of(species))
  
  
  ### Subset the iNat data, if necessary
  if (subset_inat < 1) {
    inat_dat <- sample_frac(inat_dat, size = subset_inat)
  }
  
  ### Drop grid cells that contain NA covariate data
  cells_to_drop <- unique(xdat_allS2$scale3_grid_ID[rowSums(is.na(xdat_allS2)) > 0])
  
  inat_dat <- inat_dat %>% 
    filter(!scale3_grid_ID %in% cells_to_drop) %>% 
    mutate(inat_datcell_ID = as.numeric(as.factor(inat_datcell_ID))) %>% 
    arrange(inat_datcell_ID)
  # xdat_allS2 <- xdat_allS2 %>% 
  #   filter(!scale3_grid_ID %in% cells_to_drop)
  # 
  
  # Update xdat so it only includes necessary s2 cells
  xdat_allS2 <- xdat_allS2 %>%
    left_join(inat_dat[, c("scale3_grid_ID", "inat_datcell_ID")], 
              by = "scale3_grid_ID") %>% 
    filter(!is.na(inat_datcell_ID)) %>% 
    arrange(inat_datcell_ID)
  
  # prepare start/end vectors for inat_dat
  S2toS3_start <- S2toS3_end <- numeric(nrow(inat_dat))
  for (i in 1:length(S2toS3_start)) {
    S2toS3_start[i] <- min(which(xdat_allS2$inat_datcell_ID == i))
    S2toS3_end[i] <-   max(which(xdat_allS2$inat_datcell_ID == i))
  }
  
  if (any(S2toS3_start == S2toS3_end)) {
    warning("At least one S3 cell has only one S2 cell")
  }
  
  # Scale xdat covariates according to CT scaling factors
  for (i in 1:length(covars_toscale)) {
    xdat_allS2[[paste0(covars_toscale[i], "_scaled")]] <- as.numeric(
      (unlist(xdat_allS2[, covars_toscale[i]]) - scaling_factors$mean[i]) / 
        scaling_factors$sd[i]
    )
  }
  
  
  xdat_allS2_mtx <- as.matrix(xdat_allS2[, occ_covars])
  xdat_allS2_backup <- xdat_allS2
  
  ###### [museum data processing] ######
  # xdat_allS2_museum <- covar_values %>% 
  #   rename(scale2_grid_ID = grid_cell) %>% 
  #   filter(scale2_grid_ID %in% cells_inRange$scale2_grid_ID) %>% 
  #   left_join(cells_inRange[, c("scale2_grid_ID", "scale3_grid_ID")])
  # 
  # museum_dat <- grid_counts_museum %>% 
  #   rename(scale3_grid_ID = grid_cell) %>% 
  #   filter(scale3_grid_ID %in% xdat_allS2_museum$scale3_grid_ID) %>% 
  #   arrange(scale3_grid_ID) %>% 
  #   mutate(scale3_new_ID = as.numeric(as.factor(scale3_grid_ID))) %>% 
  #   select(scale3_grid_ID, effort = n, scale3_new_ID, all_of(species))
  # 
  # # browser()
  # 
  # if (subset_inat < 1) {
  #   museum_dat <- sample_frac(museum_dat, size = subset_inat) %>% 
  #     mutate(scale3_new_ID = as.numeric(as.factor(scale3_new_ID)))
  # }
  # 
  # xdat_allS2_museum <- xdat_allS2_museum %>%
  #   left_join(museum_dat[, c("scale3_grid_ID", "scale3_new_ID")], 
  #             by = "scale3_grid_ID") %>% 
  #   filter(!is.na(scale3_new_ID)) %>% 
  #   arrange(scale3_new_ID)
  # 
  # S2toS3_start_museum <- S2toS3_end_museum <- numeric(nrow(museum_dat))
  # for (i in 1:length(S2toS3_start_museum)) {
  #   S2toS3_start_museum[i] <- min(which(xdat_allS2_museum$scale3_new_ID == i))
  #   S2toS3_end_museum[i] <-   max(which(xdat_allS2_museum$scale3_new_ID == i))
  # }
  # 
  # if (any(S2toS3_start_museum == S2toS3_end_museum)) {
  #   warning("At least one S3 cell has only one S2 cell")
  # }
  # 
  # for (i in 1:length(covars_toscale)) {
  #   xdat_allS2_museum[[paste0(covars_toscale[i], "_scaled")]] <- as.numeric(
  #     (unlist(xdat_allS2_museum[, covars_toscale[i]]) - scaling_factors$mean[i]) / 
  #       scaling_factors$sd[i]
  #   )
  # }
  # xdat_allS2_museum_backup <- xdat_allS2_museum
  # xdat_allS2_museum <- xdat_allS2_museum[, occ_covars]
  # xdat_allS2_museum[is.na(xdat_allS2_museum)] <- 0
  # 
  #### Loop over folds.
  # We need to rebuild the model each time because the likelihood structure changes
  # (unavoidable due to diff. number of obs. per camera deployment)
  samples_list <- list()
  summary_list <- list()
  summary_list_spatial <- list()
  
  nfolds_for_loop <- ifelse(kfold_mode == "single", 1, nfolds)
  
  
  ###### Prepare adjacency info ######
  
  # Get the adjacency matrix
  this_adj_info <- continent_adj_info
  # this_adj_info <- adj_info_list[[which(taxon_key$common_name_clean == species)]]
  
  # Get the spatcell for each CT deployment
  spatcell_CT <- occ.covs.all$scale4_grid_ID
  
  # Get the spatcell for each iNat cell
  inat_wspec <- left_join(inat_dat, distinct(grid_translator_wspec[, c("scale3_grid_ID", "scale4_grid_ID", 
                                                                       "spec_range_gridID")]))
  spatcell_inat <- inat_wspec$scale4_grid_ID
  
  ###### Prepare ecoregion info ######
  grid_transl_filtered <- grid_translator_wspec %>% 
    filter(scale4_grid_ID %in% unique(c(spatcell_CT, spatcell_inat))) %>% 
    distinct(scale4_grid_ID, spec_range_gridID)
  
  # Correct ecoregion values to exclude any for which we have no obs.
  ecoregion_updated <- ecoregion_df %>% 
    left_join(grid_transl_filtered, by = "scale4_grid_ID") %>% 
    select(spec_range_gridID, ecoregion) %>% 
    filter(!is.na(spec_range_gridID)) %>% 
    arrange(spec_range_gridID) 
  
  ecoregion_transl <- ecoregion_updated %>% 
    count(ecoregion) %>% 
    arrange(-n) %>% 
    mutate(new_ecoregion = row_number())
  
  ecoregion_updated <- ecoregion_updated %>% 
    left_join(ecoregion_transl, by = "ecoregion") %>% 
    mutate(ecoregion = new_ecoregion) %>% 
    select(-new_ecoregion) %>% select(-n) %>% 
    left_join(grid_transl_filtered)
  
  
  ecoregion_final <- data.frame(
    scale4_grid_ID = 1:length(this_adj_info$num)
  ) %>% 
    left_join(ecoregion_updated, by = "scale4_grid_ID")
  
  stopifnot(nrow(ecoregion_final) == length(this_adj_info$num))
  stopifnot(sum(is.na(ecoregion_final$ecoregion[spatcell_CT]))== 0)
  stopifnot(sum(is.na(ecoregion_final$ecoregion[spatcell_inat])) == 0)
  
  ###### Create the ecoregion interaction covariate data ######
  ecoregion_levels <- ecoregion_final$ecoregion %>% 
    unique() %>% 
    sort()
  
  # browser()
  # Each row matches xdat
  # ecoregion_interaction_dat_occcovs <- make_ecoregion_interaction_dat(covs_df = occ.covs.all,
  #                                                                 occ_covars,
  #                                                                 ecoregion_final,
  #                                                                 ecoregion_levels)
  # Each row matches xdat_allS2
  # ecoregion_interaction_dat_allS2 <- make_ecoregion_interaction_dat(xdat_allS2,
  #                                                               occ_covars,
  #                                                               ecoregion_final,
  #                                                               ecoregion_levels)
  
  # browser()
  #### Loop over folds ####
  for (fold_i in 1:nfolds_for_loop) {
    
    # Set up holdout data for this fold
    holdout_IDs <- folds_df$camera_ID[folds_df$fold == fold_i]
    
    inmod_inds_y   <- which(!(y.all$site_ID %in% holdout_IDs))
    holdout_inds_y <- which(y.all$site_ID %in% holdout_IDs)
    
    y_inmod <- y.all[inmod_inds_y, ]
    y_holdout <- y.all[holdout_inds_y, ]
    
    det_x_inmod <-   det.covs.all[inmod_inds_y, ]
    det_x_holdout <- det.covs.all[holdout_inds_y, ]
    
    # occ_x_inmod <- occ.covs.all[!(occ.covs.all$site_ID %in% holdout_IDs), ]
    # occ_x_holdout <- occ.covs.all[occ.covs.all$site_ID %in% holdout_IDs, ]
    
    # Get the start/end vectors for main data
    siteID_vec <- unique(y_inmod$site_ID)
    start_vec <- end_vec <- numeric(length(siteID_vec))
    
    for (i in 1:length(siteID_vec)) {
      this_id <- siteID_vec[i]
      start_vec[i] <- min(which(y_inmod$site_ID == this_id))
      end_vec[i] <-   max(which(y_inmod$site_ID == this_id))
    }
    
    # Get start/end vectors for holdout data
    siteID_vec_holdout <- unique(y_holdout$site_ID)
    start_vec_holdout <- end_vec_holdout <- numeric(length(siteID_vec_holdout))
    
    for (i in 1:length(siteID_vec_holdout)) {
      this_id <- siteID_vec_holdout[i]
      start_vec_holdout[i] <- min(which(y_holdout$site_ID == this_id))
      end_vec_holdout[i] <-   max(which(y_holdout$site_ID == this_id))
    }
    
    # Check collinearity, for output
    collinearity_CT <- check_collinearity(as.matrix(occ.covs.all[, occ_covars]))
    collinearity_iNat <- check_collinearity(xdat_allS2[, occ_covars])
    
    ecoregion_final$ecoregion[is.na(ecoregion_final$ecoregion)] <-
      max(ecoregion_final$ecoregion, na.rm = TRUE) + 1

    constants_list <- list(
      holdout_data = holdout_data,
      
      nobs = nrow(y_inmod),
      start = start_vec,
      end = end_vec,
      nobs_holdout = nrow(y_holdout),
      start_holdout = start_vec_holdout,
      end_holdout = end_vec_holdout,
      
      modtype = modtype,
      inat_disttype = inat_disttype,
      
      nLamBeta = length(occ_covars),
      nGamma = length(det_covars),
      
      ncameras_indat = length(unique(y_inmod$site_ID)),
      ncameras_holdout = length(unique(y_holdout$site_ID)),
      deployment_ID = y_inmod$site_ID,
      deployment_ID_holdout = y_holdout$site_ID, # Should max at ncameras_all
      ndepl_total = max(c(y_inmod$site_ID, y_holdout$site_ID)),
      
      S2toS3_start = S2toS3_start, 
      S2toS3_end = S2toS3_end,
      nINatCell = nrow(inat_dat),
      
      # 
      # S2toS3_start_museum = S2toS3_start_museum, 
      # S2toS3_end_museum = S2toS3_end_museum,
      # nMuseumCell = nrow(museum_dat),
      
      adj = this_adj_info$adj,
      num = this_adj_info$num,
      adj_len = length(this_adj_info$adj),
      nSpatCell = length(this_adj_info$num),
      spatcell = spatcell_CT,
      spatcell_inat = spatcell_inat,
      # spatcell_museum = spatcell_museum,
      
      inat_effort = as.numeric(unlist(inat_dat[, "effort"])),
      # museum_effort = as.numeric(unlist(museum_dat[, "effort"])),
      
      xdat_allS2  = xdat_allS2_mtx,
      # xdat_allS2_museum  = xdat_allS2_museum[, occ_covars],
      
      ecoregion = ecoregion_final$ecoregion,
      nEcoregion = max(ecoregion_final$ecoregion, na.rm = TRUE) - 1,
      
      # ecoregion_interaction_data = ecoregion_interaction_dat_occcovs$df,
      # xdat_allS2_ecoregionDat = ecoregion_interaction_dat_allS2$df,
      # ecoregion_int_groups = ecoregion_interaction_dat_occcovs$ecoregion_int_groups,
      # nEcoregionInts = length(ecoregion_interaction_dat_occcovs$ecoregion_int_groups),
      
      cc = 0
    )
    
    data_list <- list(
      inat_ct     = as.numeric(unlist(inat_dat[, species])),
      # museum_ct     = as.numeric(unlist(inat_dat[, species])),
      y_holdout =    as.numeric(y_holdout$value),
      y =            as.numeric(y_inmod$value),
      
      xdat =         as.matrix(occ.covs.all[, occ_covars]),
      wdat =         as.matrix(det_x_inmod[, det_covars]),
      wdat_holdout = as.matrix(det_x_holdout[, det_covars])
    )
    
    if (return_data_summary) {
      return(list(
        constants = constants_list,
        data = data_list
      ))
    }
    
    ##### Build and compile the model #####
    build_start_time <- Sys.time()
    # browser()
    mod <- nimbleModel(code = model_code,
                       constants = constants_list,
                       data = data_list,
                       inits = get_inits(modtype, occ_covars, det_covars,
                                         spatial_model, this_adj_info,
                                         nEcoregion = max(ecoregion_updated$ecoregion) + 1,
                                         ndepl = max(c(y_inmod$site_ID, y_holdout$site_ID))),
                       calculate = F)
    cmod <- compileNimble(mod)
    
    # browser()
    # Run, retrieve samples and summary
    
    #### Put together the custom samplers ####
    mcmcConf <- configureMCMC(cmod)
    
    if (spatial_model %in% c("SVCs", "CAR") && modtype != "separate") {
      mcmcConf$removeSamplers(c("theta0", "theta1", "lambda_intercept"))
      # Add a slice sampler for theta0, theta1, and intercept.
      # I'm keeping the base samplers as well, I'm ok with these getting double-sampled
      mcmcConf$addSampler(type = "AF_slice", target = c("theta0", "theta1",
                                                        "lambda_intercept"))
    }
    
    if (spatial_model == "SVCs") {
      configureRJ(mcmcConf,
                  targetNodes = 'sigma_CAR_0',
                  indicatorNodes = 'CAR_gamma',
                  control = list(mean = 0, scale = 1))
      configureRJ(mcmcConf,
                  targetNodes = 'sigma_eco_0',
                  indicatorNodes = 'eco_gamma',
                  control = list(mean = 0, scale = 1))
    }
    
    # Add a custom sampler that will track the log posterior density of the model
    mcmcConf$removeSamplers('logDens') ## remove sampler assigned to 'logDens'
    mcmcConf$addSampler(target = 'logDens', type = 'sumLogProb') ## add our custom sampler
    
    # browser()
    monitors_list <- get_monitors(modtype, inat_disttype, spatial_model)
    
    mcmcConf$setMonitors(c(monitors_list$monitors, "logDens"))
    mcmcConf$setMonitors2(monitors_list$monitors2)
    
    # Build and compile the MCMC
    mcmc <- buildMCMC(mcmcConf)
    cmcmc <- compileNimble(mcmc)
    
    ##### Run the MCMC #####
    # browser()
    mcmc_start_time <- Sys.time()
    samples_list[[fold_i]] <- 
      runMCMC(cmcmc, niter = ni, nburnin = nb, thin = nt, 
              nchains = nc, thin2 = nt2,
              samplesAsCodaMCMC = TRUE, 
              inits = get_inits(modtype, occ_covars, det_covars,
                                spatial_model, this_adj_info,
                                nEcoregion = max(ecoregion_updated$ecoregion) + 1,
                                ndepl = max(c(y_inmod$site_ID, y_holdout$site_ID))))
    
    summary_list[[fold_i]] <- MCMCvis::MCMCsummary(samples_list[[fold_i]]$samples)
    summary_list[[fold_i]]$param <- rownames(summary_list[[fold_i]])
    summary_list[[fold_i]]$fold <- fold_i
    summary_list[[fold_i]]$species <- species
    
    summary_list_spatial[[fold_i]] <- MCMCvis::MCMCsummary(samples_list[[fold_i]]$samples2)
    summary_list_spatial[[fold_i]]$param <- rownames(summary_list_spatial[[fold_i]])
    summary_list_spatial[[fold_i]]$fold <- fold_i
    summary_list_spatial[[fold_i]]$species <- species
  }
  
  # browser()
  
  param_inds <- colnames(samples_list[[1]]$samples[[1]])
  logsc_index <- which(param_inds == "logarithmic_score")
  
  
  # logarithmic_score <- as.numeric(unlist(
  #   lapply(samples_list, function(x) lapply(x$samples, function(y) as.numeric(y[, logsc_index])))
  # )) %>% mean()
  logarithmic_score <- NA
  end_time <- Sys.time()
  
  lambda_covar_indices <- c(grep("lambda_beta", param_inds), 
                            grep("p_alpha", param_inds),
                            grep("inat_beta", param_inds),
                            grep("museum_beta", param_inds),
                            grep("CAR_gamma", param_inds),
                            grep("eco_gamma", param_inds),
                            grep("sigma_CAR", param_inds),
                            grep("sigma_eco", param_inds))
  det_covar_indices <- grep("p_gamma", param_inds)
  
  for (i in 1:length(summary_list)) {
    summary_list[[i]]$parname <- NA
    summary_list[[i]]$parname[lambda_covar_indices] <- gsub("_scaled", "", occ_covars)
    summary_list[[i]]$parname[det_covar_indices] <- gsub("_scaled", "", det_covars)
  }
  
  ###### Return ######
  saveRDS(list(
    time_taken = end_time - start_time,
    build_time_taken = mcmc_start_time - build_start_time,
    mcmc_time_taken = end_time - mcmc_start_time,
    thin = nt, thin2 = nt2,
    samples_list = samples_list,
    summary_list = summary_list,
    spatial_model = spatial_model,
    summary_list_spatial = summary_list_spatial,
    scaling_factors = scaling_factors,
    scaling_factors_det = scaling_factors_det,
    nfolds = nfolds, folds_by = folds_by,
    modtype = modtype, 
    species = species,
    subset = subset,
    subset_inat = subset_inat,
    occ_covars = occ_covars, 
    det_covars = det_covars,
    call = match.call(),
    collinearity_iNat = collinearity_iNat,
    collinearity_CT = collinearity_CT
  ), outfile)
  
}
