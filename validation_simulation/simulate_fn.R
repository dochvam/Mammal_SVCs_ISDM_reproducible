###############################################################################
#'
#' validation_simulation/simulate_fn.R
#' Author: Ben R. Goldstein
#' Description: Code for simulating data for model validation exercise.
#'    Includes a function that replicates the model-building process, following
#'    which the NIMBLE model is used to simulate data
#'
###############################################################################


source("main_code_lineage/model_fn.R")
# source("main_code_lineage/main_data_prep.R")

if (!dir.exists("intermediate/validation_sim")) {
  dir.create("intermediate/validation_sim")
}


run_simulation <- function(scenario, 
                           species,
                           seed,
                           iter,
                           ni = 1000, 
                           nb = 0,
                           nc = 2, 
                           nt = 1, 
                           nt2 = 10) {
  
  stopifnot(scenario %in% 1:3)
  
  mod <- build_model_only(
    nfolds = 0,
    inat_disttype = "NB",
    modtype = "joint",
    species = species,
    spatial_model = "SVCs",
    seed = seed, subset_CT = 1,
    subset_inat = 1
  )

  mod$y[] <- NA
  mod$y_holdout[] <- NA
  mod$inat_ct[] <- NA
  
  # Set true parameter values
  
  ### Intensity intercepts
  mod$p_intercept <- 0.2
  mod$lambda_intercept <- 0.3
  
  ### Indicator variables
  true_spatial_df <- switch(scenario,
                            data.frame(
                              CAR = c(1,1,1,0,0,0,0),
                              eco = c(0,0,0,1,1,0,0),
                              lin = c(0,0,0,0,0,1,1)
                            ),
                            data.frame(
                              CAR = c(1,1,1,0,0,0,0),
                              eco = c(1,1,0,0,0,0,0),
                              lin = c(0,0,0,0,1,1,0)
                            ),
                            data.frame(
                              CAR = c(1,1,1,0,0,0,0),
                              eco = c(0,0,0,0,1,1,0),
                              lin = c(1,1,0,0,0,0,0)
                            )) %>% 
    mutate(param = gsub("_scaled", "", occ_covars)) %>% 
    sample_frac(size = 1, replace = FALSE)
  
  mod$CAR_gamma <- true_spatial_df$CAR
  mod$eco_gamma <- true_spatial_df$eco
  mod$lin_gamma <- true_spatial_df$lin
  
  ### intensity coeffs
  mod$lambda_beta_mean <- runif(length(mod$lambda_beta_mean), -1, 1)
  
  ### intensity spatial fields
  this_adj_info <- continent_adj_info
  
  mod$spat_ranef[] <- simulate_dcar_normal(this_adj_info$adj,
                                           this_adj_info$num,
                                           precision = 1)
  mod$spat_magnitude_intercept <- 0.5
  
  for (i in 1:length(occ_covars)) {
    mod$sigma_CAR_0[] <- 1
    mod$sigma_lin_0[] <- 1
    mod$sigma_eco_0[] <- 1
    
    mod$CAR_effect[, i] <- simulate_dcar_normal(this_adj_info$adj,
                                                 this_adj_info$num,
                                                 precision = 1)
    mod$lin_effect[, i] <- sample(seq(-1.5, 1.5, length.out = nrow(mod$lin_effect)))
    mod$eco_effect[, i] <- sample(seq(-1.5, 1.5, length.out = nrow(mod$eco_effect)))
  }

  ### detection coeffs
  mod$p_alpha <- runif(length(mod$p_alpha), -1, 1)
  mod$p_gamma <- runif(length(mod$p_gamma), -1, 1)
  
  mod$det_ranef_sd <- 0.2
  mod$det_ranef <- rnorm(length(mod$det_ranef), 0, 0.2)
  
  ### iNat parameters
  mod$theta0 <- -5
  mod$theta1 <- 0.5
  mod$overdisp_inat <- 0.1
  
  # Calculate everything deterministic in the model
  mod$calculate(mod$getParents("y", upstream = TRUE))
  mod$calculate(mod$getParents("y_holdout", upstream = TRUE))
  mod$calculate(mod$getParents("inat_ct", upstream = TRUE))
  
  # Simulate y
  mod$simulate("y", includeData = TRUE)
  mod$simulate("y_holdout", includeData = TRUE)
  
  # simulate inat_ct
  mod$simulate("inat_ct", includeData = TRUE)
  
  # Make sure these are still recognized as data
  mod$setData("y")
  mod$setData("y_holdout")
  mod$setData("inat_ct")
  

  cmod <- compileNimble(mod)
  
  mcmcConf <- configureMCMC(cmod)

  mcmcConf$removeSamplers(c("theta0", "theta1", "lambda_intercept"))
  # Add a slice sampler for theta0, theta1, and intercept.
  # I'm keeping the base samplers as well, I'm ok with these getting double-sampled
  mcmcConf$addSampler(type = "AF_slice", target = c("theta0", "theta1",
                                                    "lambda_intercept"))
  
  configureRJ(mcmcConf,
              targetNodes = 'sigma_CAR_0',
              indicatorNodes = 'CAR_gamma',
              control = list(mean = 0, scale = 1))
  configureRJ(mcmcConf,
              targetNodes = 'sigma_lin_0',
              indicatorNodes = 'lin_gamma',
              control = list(mean = 0, scale = 1))
  configureRJ(mcmcConf,
              targetNodes = 'sigma_eco_0',
              indicatorNodes = 'eco_gamma',
              control = list(mean = 0, scale = 1))
  # Add a custom sampler that will track the log posterior density of the model
  mcmcConf$removeSamplers('logDens') ## remove sampler assigned to 'logDens'
  mcmcConf$addSampler(target = 'logDens', type = 'sumLogProb') ## add our custom sampler
  
  
  monitors_list <- get_monitors(modtype, inat_disttype, spatial_model = "SVCs")
  
  mcmcConf$setMonitors(c(monitors_list$monitors, "logDens"))
  mcmcConf$setMonitors2(monitors_list$monitors2)
  
  # Build and compile the MCMC
  mcmc <- buildMCMC(mcmcConf)
  cmcmc <- compileNimble(mcmc)
  
  ##### Run the MCMC #####
  samples_list <- 
    runMCMC(cmcmc, niter = ni, nburnin = nb, thin = nt, 
            nchains = nc, thin2 = nt2,
            samplesAsCodaMCMC = TRUE, 
            inits = get_inits(modtype = "joint", 
                              occ_covars, det_covars,
                              spatial_model = "SVCs", 
                              this_adj_info,
                              ndepl = length(mod$det_ranef),
                              nLineage = nrow(mod$lin_effect),
                              nEcoregion = nrow(mod$eco_effect)))
  
  summary_df <- MCMCvis::MCMCsummary(samples_list$samples)
  
  summary_df <- summary_df %>% 
    mutate(species = species,
           scenario = scenario,
           iter = iter,
           param = rownames(summary_df),
           parnum = parse_number(param))
  
  true_spatial_long <- true_spatial_df %>% 
    mutate(parnum = row_number()) %>% 
    pivot_longer(cols = c("CAR", "lin", "eco")) %>% 
    rename(type = name, truth = value) %>% 
    select(-param)
  
  result_df <- summary_df %>% 
    filter((grepl("_gamma", param) & !grepl("^p_gamma", param)) |
             grepl("logDens", param)) %>% 
    mutate(type = substr(param, 1, 3)) %>% 
    left_join(true_spatial_long) %>% 
    mutate(species = species,
           scenario = scenario,
           )
  
  outfile <- file.path("intermediate", "validation_sim",
                       paste0("VS_", species, "_S", scenario, "_I", iter, ".RDS"))
    
  saveRDS(list(
    result = result_df,
    species = species,
    scenario = scenario,
    iter = iter,
    seed = seed,
    samples = samples_list
  ), file = outfile)
}



build_model_only <- function(nfolds, 
                             kfold_mode = "kfold",
                             folds_by = "camera", 
                             modtype,
                             inat_disttype = "Pois",
                             spatial_model = "None",
                             species,
                             seed,
                             subset_CT = 1,
                             subset_inat = 1) {
  
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
  
  source("main_code_lineage/model_code_main.R")
  
  if (spatial_model == "CAR") {
    model_code <- joint_model_CAR
  } else if (spatial_model == "SVCs") {
    model_code <- joint_model_SVCs_asRanefs
  } else {
    model_code <- joint_model_nonspatial
  }
  
  ### Load in and format the data
  source("main_code_lineage/main_data_prep.R")
  set.seed(seed)
  
  # Check that species / modtype is valid
  stopifnot(species %in% taxon_key$common_name_clean)
  if (spatial_model == 'SVCs') {
    stopifnot(species %in% lineage_key$common_name_clean)
  }
  
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
  
  # If SVC, drop sites without lineages
  scale4_cell_wLineage <- lineage_df$scale4_grid_ID[!is.na(lineage_df[[species]])]
  if (spatial_model == "SVCs") {
    
    sites_without_lineage <- occ.covs.all %>% 
      filter(!(scale4_grid_ID %in% scale4_cell_wLineage)) %>% 
      select(site_ID)
    
    sites_to_drop <- 
      bind_rows(sites_to_drop, sites_without_lineage) %>% 
      distinct()
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
  
  
  ### If SVCs, drop iNat grid cells that don't have lineages
  if (spatial_model == "SVCs") {
    scale3_cells_wLineage <- grid_translator_wspec %>% 
      filter(scale4_grid_ID %in% scale4_cell_wLineage)
    
    inat_dat <- inat_dat %>% 
      filter(scale3_grid_ID %in% scale3_cells_wLineage$scale3_grid_ID)
  }
  
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
  
  # Get the spatcell for each CT deployment
  spatcell_CT <- occ.covs.all$scale4_grid_ID
  
  # Get the spatcell for each iNat cell
  inat_wspec <- left_join(inat_dat, distinct(grid_translator_wspec[, c("scale3_grid_ID",
                                                                       "scale4_grid_ID",
                                                                       "spec_range_gridID")]))
  spatcell_inat <- inat_wspec$scale4_grid_ID
  
  ###### Prepare lineage info ######
  
  # Filter down to only the S4 cells containing data
  grid_transl_filtered <- grid_translator_wspec %>% 
    filter(scale4_grid_ID %in% unique(c(spatcell_CT, spatcell_inat))) %>% 
    distinct(scale4_grid_ID, spec_range_gridID)
  
  # Identify the species we need
  colnames(lineage_df)[colnames(lineage_df) == species] <- "target_lineage"
  
  # Correct lineage values to exclude any for which we have no obs.
  lineage_updated <- lineage_df %>% 
    left_join(grid_transl_filtered) %>% 
    select(scale4_grid_ID, target_lineage) %>% 
    filter(!is.na(target_lineage)) %>% 
    arrange(scale4_grid_ID) 
  
  lineage_transl <- lineage_updated %>% 
    count(target_lineage) %>% 
    arrange(-n) %>% 
    mutate(new_target_lineage = row_number())
  
  lineage_updated <- lineage_updated %>% 
    left_join(lineage_transl, by = "target_lineage") %>% 
    mutate(target_lineage = new_target_lineage) %>% 
    select(-new_target_lineage) %>% select(-n)
  
  lineage_final <- data.frame(
    scale4_grid_ID = 1:length(this_adj_info$num)
  ) %>% 
    left_join(lineage_updated, by = "scale4_grid_ID")
  
  # browser()
  
  stopifnot(nrow(lineage_final) == length(this_adj_info$num))
  stopifnot(sum(is.na(lineage_final$target_lineage[spatcell_CT]))== 0)
  stopifnot(sum(is.na(lineage_final$target_lineage[spatcell_inat])) == 0)
  
  ###### Create the lineage interaction covariate data ######
  xdat_allS2 <- xdat_allS2 %>% 
    left_join(grid_translator_wspec[, c("scale2_grid_ID", "spec_range_gridID")])
  lineage_levels <- lineage_final$target_lineage %>% 
    unique() %>% 
    sort()
  
  # browser()
  # Each row matches xdat
  # lineage_interaction_dat_occcovs <- make_lineage_interaction_dat(covs_df = occ.covs.all,
  #                                                                 occ_covars,
  #                                                                 lineage_final,
  #                                                                 lineage_levels)
  # # Each row matches xdat_allS2
  # lineage_interaction_dat_allS2 <- make_lineage_interaction_dat(xdat_allS2,
  #                                                               occ_covars,
  #                                                               lineage_final,
  #                                                               lineage_levels)
  ###### Prepare ecoregion info ######
  
  # Correct lineage values to exclude any for which we have no obs.
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
  
  fold_i <- 1

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
  lineage_final$target_lineage[is.na(lineage_final$target_lineage)] <-
    max(lineage_final$target_lineage, na.rm = TRUE) + 1
  
  
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
    
    lineage = lineage_final$target_lineage,
    nLineage = max(lineage_final$target_lineage, na.rm = TRUE) - 1,
    
    # lineage_interaction_data = lineage_interaction_dat_occcovs$df,
    # xdat_allS2_lineageDat = lineage_interaction_dat_allS2$df,
    # lineage_int_groups = lineage_interaction_dat_occcovs$lineage_int_groups,
    # nLineageInts = length(lineage_interaction_dat_occcovs$lineage_int_groups),
    
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
  
  ##### Build and compile the model #####
  build_start_time <- Sys.time()
  mod <- nimbleModel(code = model_code,
                     constants = constants_list,
                     data = data_list,
                     inits = get_inits(modtype, occ_covars, det_covars,
                                       spatial_model, this_adj_info,
                                       ndepl = constants_list$ndepl_total,
                                       nLineage = max(lineage_updated$target_lineage) + 1,
                                       nEcoregion = max(ecoregion_updated$ecoregion) + 1),
                     calculate = F)
  
  mod
}



simulate_dcar_normal <- function(adj, num, precision) {
  n <- length(num)  # Number of locations
  
  # Create a sparse precision matrix Q
  Q <- matrix(0, n, n)
  
  # Fill the precision matrix based on adjacency information
  index <- 1
  for (i in 1:n) {
    neighbors <- adj[index:(index + num[i] - 1)]
    Q[i, neighbors] <- -1
    Q[neighbors, i] <- -1
    Q[i, i] <- num[i]
    index <- index + num[i]
  }
  
  # Multiply by precision to get the precision matrix for the CAR model
  Q <- precision * Q
  
  E       <- eigen(Q)
  E$vec   <- E$vec[,-n]
  E$val   <- E$val[-n]
  Y       <- E$vec%*%rnorm(n-1,0,1/sqrt(E$val))
  
  return(Y)
}
