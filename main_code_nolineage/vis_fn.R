###############################################################################
#'
#' main_code_nolineage/vis_fn.R
#' Author: Ben R. Goldstein
#' Description: Helper code for visualization script. Used for both
#'    types of models.
#'
###############################################################################


summarize_conf_by_region <- function(df) {
  df %>% 
    group_by(StateProv) %>% 
    summarize(
      mean_occ_mean = mean(cloglog_occ_mean, na.rm = TRUE),
      mean_occ_sd = mean(cloglog_occ_sd, na.rm = TRUE)
    )
}

generate_spatpred_iSDM <- function(res, species) {
  
  scaling_factors <- res$scaling_factors
  #get_scaling_factors(species, seed = 0159876, subset = 1)
  
  gridcol <- paste0("GRID_ID_", species)
  spec_range_cells <- grid_translator_wspec[!is.na(grid_translator_wspec[, gridcol]), ]
  
  # Prepare covariate grid
  covar_values_scaled <- covar_values %>% 
    filter(grid_cell %in% spec_range_cells$scale2_grid_ID &
             !is.na(grid_cell)) %>% 
    left_join(grid_translator_wspec, by = c("grid_cell" = "scale2_grid_ID")) %>% 
    mutate(Pop_den_sqrt = sqrt(Pop_den))
  covar_values_scaled$Agriculture <- covar_values_scaled$Pastureland + 
    covar_values_scaled$Arable
  for (i in 1:nrow(scaling_factors)) {
    covar_values_scaled[, paste0(scaling_factors$covar[i], "_scaled")] <- 
      (covar_values_scaled[, scaling_factors$covar[i]] - scaling_factors$mean[i]) /
      scaling_factors$sd[i]
  }
  covar_values_scaled$Intercept <- 1
  
  covar_values_scaled <- covar_values_scaled[, c(gridcol, "scale4_grid_ID", "grid_cell", occ_covars)] %>% 
    na.omit()
  
  intercept_column <- "lambda_intercept"
  target_str <- "lambda_beta\\["
  # Prepare samples. We need:
  # - lambda intercept
  # - intercept spatial ranef
  # - spatial lambda for all covars.
  samples2_names <- colnames(res$samples_list[[1]]$samples2[[1]])
  
  target_cols_spatranef <- paste0("spat_ranef[", 
                                  unique(covar_values_scaled$scale4_grid_ID), "]")
  target_cols_spatranef <- target_cols_spatranef[
    target_cols_spatranef %in% samples2_names
  ]
  samples_spat_ranef <- res$samples_list %>% 
    lapply(function(x) {
      as.data.frame(as.matrix(x[[2]]))
    }) %>% 
    bind_rows() %>% 
    select(all_of(target_cols_spatranef)) %>% 
    as.matrix()
  index_df_spatranef <- data.frame(
    index = 1:length(target_cols_spatranef),
    param = target_cols_spatranef,
    scale4_grid_ID = parse_number(target_cols_spatranef)
  )
  
  target_cols <- samples2_names[grepl("lambda_beta", samples2_names)]
  
  samples_betas <- res$samples_list %>% 
    lapply(function(x) {
      as.data.frame(as.matrix(x[[2]]))
    }) %>% 
    bind_rows() %>% 
    select(all_of(target_cols)) %>% 
    as.matrix()
  index_df_beta <- data.frame(
    index = 1:length(target_cols),
    param = target_cols
  ) %>% 
    mutate(scale4_grid_ID = parse_number(param),
           covar_ind = parse_number(str_sub(param, start = nchar(param) - 3,
                                            end = nchar(param))))
  
  
  samples_intercept <- res$samples_list %>% 
    lapply(function(x) {
      as.data.frame(as.matrix(x[[1]]))
    }) %>% 
    bind_rows() %>% 
    select(all_of(intercept_column)) %>% 
    as.matrix()
  samples_spatmag <- res$samples_list %>% 
    lapply(function(x) {
      as.data.frame(as.matrix(x[[1]]))
    }) %>% 
    bind_rows() %>% 
    select("spat_magnitude_intercept") %>% 
    as.matrix()
  
  
  # Thin the intercept to match monitors2
  samples_intercept <- samples_intercept[
    1:nrow(samples_betas) * (nrow(samples_intercept) / nrow(samples_betas))
  ]
  samples_spatmag <- samples_spatmag[
    1:nrow(samples_betas) * (nrow(samples_spatmag) / nrow(samples_betas))
  ]
  
  stopifnot(nrow(samples_intercept) == nrow(samples_betas))
  
  cov_mtx <- as.matrix(covar_values_scaled[, occ_covars])
  cloglog_occ_mtx <- matrix(NA, 
                            nrow = nrow(covar_values_scaled),
                            ncol = nrow(samples_betas))
  pb <- progress::progress_bar$new(total = nrow(covar_values_scaled))
  for (i in 1:nrow(covar_values_scaled)) {
    pb$tick()
    lambda_indexes <- index_df_beta$index[
      index_df_beta$scale4_grid_ID == covar_values_scaled$scale4_grid_ID[i]
    ]
    ranef_ind <- index_df_spatranef$index[
      index_df_spatranef$scale4_grid_ID == covar_values_scaled$scale4_grid_ID[i]
    ]
    
    cloglog_occ_mtx[i, ] <- as.numeric(
      log(samples_intercept) + 
        (samples_betas[, lambda_indexes] %*% cov_mtx[i, ]) +
        samples_spat_ranef[, ranef_ind] * samples_spatmag
    )
    
  }
  
  covar_values_scaled$cloglog_occ_mean <- rowMeans(cloglog_occ_mtx)
  covar_values_scaled$cloglog_occ_sd <- apply(cloglog_occ_mtx, 1, sd)
  
  return(covar_values_scaled)
}

generate_spatpred <- function(result, species, inat_betas = F, museum_betas = F,
                              thin, CIs = TRUE, always_use_lambda_int = F) {
  if (!exists("grid_translator_wspec")) source("main_code_lineage/main_data_prep.R")
  
  
  scaling_factors <- result$scaling_factors
  #get_scaling_factors(species, seed = 0159876, subset = 1)
  
  gridcol <- paste0("GRID_ID_", species)
  spec_range_cells <- grid_translator_wspec[!is.na(grid_translator_wspec[, gridcol]), ]
  
  covar_values_scaled <- covar_values %>% 
    filter(grid_cell %in% spec_range_cells$scale2_grid_ID) %>% 
    left_join(grid_translator_wspec, by = c("grid_cell" = "scale2_grid_ID")) %>% 
    mutate(Pop_den_sqrt = sqrt(Pop_den))
  covar_values_scaled$Agriculture <- covar_values_scaled$Pastureland + 
    covar_values_scaled$Arable
  for (i in 1:nrow(scaling_factors)) {
    covar_values_scaled[, paste0(scaling_factors$covar[i], "_scaled")] <- 
      (covar_values_scaled[, scaling_factors$covar[i]] - scaling_factors$mean[i]) /
      scaling_factors$sd[i]
  }
  covar_values_scaled$Intercept <- 1
  
  if (inat_betas) {
    target_str <- "inat_beta["
    target_str_sp <- "inat_spat_ranef["
    target_str_spatmag <- "inat_spat_magnitude"
    intercept_column <- "log_inat_intercept"
    int_link <- exp
  } else if (museum_betas) {
    target_str <- "museum_beta["
    target_str_sp <- "museum_spat_ranef["
    target_str_spatmag <- "museum_spat_magnitude"
    intercept_column <- "log_museum_intercept"
    int_link <- exp
  } else {
    target_str <- "lambda_beta["
    target_str_sp <- "spat_ranef["
    target_str_spatmag <- "spat_magnitude"
    intercept_column <- "lambda_intercept"
    int_link <- function(x) x
  }
  
  if (always_use_lambda_int) {
    intercept_column <- "lambda_intercept"
    int_link <- function(x) x
  }
  
  # Figure out if we're dealing with an explicitly spatial model
  spatial_model_fit <- any(grepl("spat_ranef", result$summary_list_spatial[[1]]$param))
  
  if (spatial_model_fit) {
    
    thin_ratio <- result$thin2 / result$thin
    
    # Predict occupancy
    samples_intercept <- result$samples_list %>% 
      lapply(function(x) {
        as.data.frame(as.matrix(x[[1]]))
      }) %>% 
      bind_rows() %>% 
      select(all_of(intercept_column)) %>% 
      as.matrix()
    
    samples <- result$samples_list %>% 
      lapply(function(x) {
        as.data.frame(as.matrix(x[[1]]))
      }) %>% 
      bind_rows() %>% 
      select(all_of(paste0(target_str, 1:length(result$occ_covars), "]"))) %>% 
      as.matrix()
    
    samples_spatial <- result$samples_list %>% 
      lapply(function(x) {
        as.data.frame(as.matrix(x[[2]]))
      }) %>% 
      bind_rows() %>% 
      select(all_of(paste0(target_str_sp, 1:max(spec_range_cells[, gridcol]), "]"))) %>% 
      as.matrix()
    
    samples_spatial_mag <- result$samples_list %>% 
      lapply(function(x) {
        as.data.frame(as.matrix(x[[1]]))
      }) %>% 
      bind_rows() %>% 
      select(all_of(target_str_spatmag)) %>% 
      as.matrix()
    
    
    ### Start the loop
    occ_covars <- c(result$occ_covars)
    cloglog_occ_mtx <- matrix(nrow = nrow(samples)/thin_ratio, ncol = nrow(covar_values_scaled))
    covar_mtx <- as.matrix(covar_values_scaled[, occ_covars])
    spatcells <- unname(unlist(covar_values_scaled[, gridcol]))
    
    pb <- progress::progress_bar$new(total = nrow(samples)/thin_ratio)
    for (i in 1:(nrow(samples)/thin_ratio)) {
      
      pb$tick()
      cloglog_occ_mtx[i, ] <- 
        int_link(samples_intercept[i*thin_ratio, ]) +
        covar_mtx %*% samples[i*thin_ratio, ] +
        samples_spatial[i, spatcells] * samples_spatial_mag[i, 1]
      
    }
    
    covar_values_scaled$cloglog_occ_mean <- colMeans(cloglog_occ_mtx)
    covar_values_scaled$cloglog_occ_sd <- apply(cloglog_occ_mtx, 2, sd)
    
  } else {
    
    # stop("I broke non-spatial when I changed intercept to log; go fix")
    # Predict occupancy
    samples_intercept <- result$samples_list %>% 
      lapply(function(x) {
        as.data.frame(as.matrix(x[[1]]))
      }) %>% 
      bind_rows() %>% 
      select(all_of(intercept_column)) %>% 
      as.matrix()
    
    samples <- result$samples_list %>% 
      lapply(function(x) {
        as.data.frame(as.matrix(x[[1]]))
      }) %>% 
      bind_rows() %>% 
      select(all_of(paste0(target_str, 1:length(result$occ_covars), "]"))) %>% 
      as.matrix()
    
    ### Start the loop
    occ_covars <- c(result$occ_covars)
    cloglog_occ_mtx <- matrix(nrow = nrow(samples)/thin, ncol = nrow(covar_values_scaled))
    covar_mtx <- as.matrix(covar_values_scaled[, occ_covars])
    
    pb <- progress::progress_bar$new(total = nrow(samples)/thin)
    for (i in 1:(nrow(samples)/thin)) {
      
      pb$tick()
      cloglog_occ_mtx[i,] <- int_link(samples_intercept[i*thin, ]) + 
        covar_mtx %*% samples[i*thin, ]
      
    }
    
    covar_values_scaled$cloglog_occ_mean <- colMeans(cloglog_occ_mtx)
    covar_values_scaled$cloglog_occ_sd <- apply(cloglog_occ_mtx, 2, sd)
    if (CIs) {
      covar_values_scaled$cloglog_occ_LB <- apply(cloglog_occ_mtx, 2, quantile, prob = 0.025, na.rm = T)
      covar_values_scaled$cloglog_occ_UB <- apply(cloglog_occ_mtx, 2, quantile, prob = 0.975, na.rm = T)
    }
    
  }
  return(covar_values_scaled)
}

# plot_spatpred <- function(result, species, suffix, write.out = TRUE, inat_betas = FALSE,
#                           thin = 50) {
#   
#   covar_values_scaled <- generate_spatpred(result, species, suffix, inat_betas, thin)
#   cloglog_occ_mean <- covar_values_scaled$cloglog_occ_mean
#   cloglog_occ_sd <- covar_values_scaled$cloglog_occ_sd
#   #### Plot of occupancy
#   
#   occ_mean <- icloglog(cloglog_occ_mean)
#   
#   occ_mean_raster <- continental_grid_scale2
#   
#   range_cells <- which(terra::values(continental_grid_scale2) %in% covar_values_scaled$grid_cell)
#   other_cells <- which(!terra::values(continental_grid_scale2) %in% covar_values_scaled$grid_cell &
#                          !is.na(terra::values(occ_mean_raster)))
#   
#   terra::values(occ_mean_raster)[!is.na(terra::values(occ_mean_raster))] <- 0
#   terra::values(occ_mean_raster)[range_cells] <- occ_mean
#   
#   
#   p1 <- ggplot() +
#     geom_spatraster(data = occ_mean_raster, maxcell = 1e8,
#                     na.rm = T) +
#     scale_fill_viridis_c("Occupancy", na.value = NA)
#   
#   
#   ### Plot of intensity
#   
#   int_mean <- cloglog_occ_mean
#   int_mean_raster <- continental_grid_scale2
#   
#   oor_color <- min(na.rm = TRUE, int_mean) - 0.1
#   
#   terra::values(int_mean_raster)[!is.na(terra::values(continental_grid_scale2))] <- oor_color
#   terra::values(int_mean_raster)[range_cells] <- int_mean
#   
#   p2 <- ggplot() +
#     geom_spatraster(data = int_mean_raster, maxcell = 1e8) +
#     scale_fill_viridis_c("Log intensity", na.value = NA)
#   
#   
#   ### Plot of intensity
#   
#   sd_raster <- continental_grid_scale2
#   
#   terra::values(sd_raster)[!is.na(terra::values(continental_grid_scale2))] <- 0
#   terra::values(sd_raster)[range_cells] <- cloglog_occ_sd
#   
#   p3 <- ggplot() +
#     geom_spatraster(data = sd_raster, maxcell = 1e8,
#                     na.rm = T) +
#     scale_fill_viridis_c("Est. uncertainty\nin log intensity", na.value = NA)
#   
#   
#   bigplot <- gridExtra::arrangeGrob(p1, #p2, 
#                                     p3, nrow = 1, 
#                                     top = paste0(species, suffix))
#   if (write.out) {
#     ggsave(plot = bigplot, paste0("plots/mod_compare_", species, suffix, ".jpg"),
#            width = 12, height = 5)
#   }
#   return(list(
#     p1 = p1,
#     p2 = p2,
#     p3 = p3,
#     arranged = bigplot
#   ))
# }

# plot_spatpred_difference <- function(result3, result4, species, suffix, 
#                                      write.out = TRUE, thin = 50) {
#   
#   covar_values_scaled_Joint <- generate_spatpred(result3, species, suffix, inat_betas = F, thin = thin)
#   covar_values_scaled_CT    <- generate_spatpred(result4, species, suffix, inat_betas = F, thin = thin)
#   covar_values_scaled_iNat  <- generate_spatpred(result4, species, suffix, inat_betas = T, thin = thin)
#   
#   #### Diff. in occu estimate
#   
#   occ_mean_CTmJoint <- icloglog(covar_values_scaled_CT$cloglog_occ_mean) - 
#     icloglog(covar_values_scaled_Joint$cloglog_occ_mean)
#   occ_mean_iNatmJoint <- icloglog(covar_values_scaled_iNat$cloglog_occ_mean) - 
#     icloglog(covar_values_scaled_Joint$cloglog_occ_mean)
#   
#   
#   occ_raster_CTmJoint <- continental_grid_scale2
#   occ_raster_iNatmJoint <- continental_grid_scale2
#   
#   range_cells <- which(terra::values(continental_grid_scale2) %in% covar_values_scaled_Joint$grid_cell)
#   other_cells <- which(!terra::values(continental_grid_scale2) %in% covar_values_scaled_Joint$grid_cell &
#                          !is.na(terra::values(continental_grid_scale2)))
#   
#   terra::values(occ_raster_CTmJoint)[!is.na(terra::values(occ_raster_CTmJoint))] <- 0
#   terra::values(occ_raster_CTmJoint)[range_cells] <- occ_mean_CTmJoint
#   terra::values(occ_raster_iNatmJoint)[!is.na(terra::values(occ_raster_iNatmJoint))] <- 0
#   terra::values(occ_raster_iNatmJoint)[range_cells] <- occ_mean_iNatmJoint
#   
#   
#   p1 <- ggplot() +
#     geom_spatraster(data = occ_raster_iNatmJoint, maxcell = 1e8,
#                     na.rm = T) +
#     ggtitle("Diff. btw. iNat and joint") +
#     scale_fill_gradient2("dOccu", mid = "lightgray", na.value = NA) +
#     theme_minimal()
#   p2 <- ggplot() +
#     geom_spatraster(data = occ_raster_CTmJoint, maxcell = 1e8,
#                     na.rm = T) +
#     ggtitle("Diff. btw. CT and joint") +
#     scale_fill_gradient2("dOccu", mid = "lightgray", na.value = NA) +
#     theme_minimal()
#   
#   
#   
#   #### Plots of uncertainty diffs
#   
#   sd_CTmJoint <- covar_values_scaled_CT$cloglog_occ_sd - 
#     covar_values_scaled_Joint$cloglog_occ_sd
#   sd_iNatmJoint <- covar_values_scaled_iNat$cloglog_occ_sd - 
#     covar_values_scaled_Joint$cloglog_occ_sd
#   
#   sd_raster_CTmJoint <- continental_grid_scale2
#   sd_raster_iNatmJoint <- continental_grid_scale2
#   terra::values(sd_raster_CTmJoint)[!is.na(terra::values(sd_raster_CTmJoint))] <- 0
#   terra::values(sd_raster_CTmJoint)[range_cells] <- sd_CTmJoint
#   terra::values(sd_raster_iNatmJoint)[!is.na(terra::values(sd_raster_iNatmJoint))] <- 0
#   terra::values(sd_raster_iNatmJoint)[range_cells] <- sd_iNatmJoint
#   
#   
#   
#   p3 <- ggplot() +
#     geom_spatraster(data = sd_raster_iNatmJoint, maxcell = 1e8,
#                     na.rm = T) +
#     ggtitle("Diff. in cloglog-scale SD (iNat-only minus joint)") +
#     scale_fill_gradient2("d SD", mid = "lightgray", na.value = NA) +
#     theme_minimal()
#   p4 <- ggplot() +
#     geom_spatraster(data = sd_raster_CTmJoint, maxcell = 1e8,
#                     na.rm = T) +
#     ggtitle("Diff. in cloglog-scale SD (CT minus joint)") +
#     scale_fill_gradient2("d SD", mid = "lightgray", na.value = NA) +
#     theme_minimal()
#   
#   
#   
#   
#   bigplot <- arrangeGrob(p1, p2, p3, p4, nrow = 2)
#   if (write.out) {
#     ggsave(plot = bigplot, paste0("plots/diffplot_", species, suffix, ".jpg"),
#            width = 12, height = 10)
#   }
#   return(list(
#     p1 = p1,
#     p2 = p2,
#     p3 = p3,
#     p4 = p4,
#     arranged = bigplot
#   ))
# }


plot_one_joint <- function(result3, species, thin = 50) {
  covar_values_scaled_Joint <- generate_spatpred(result3, species, inat_betas = F, thin = thin, CIs = F,
                                                 always_use_lambda_int = always_use_lambda_int) %>% 
    filter(!is.na(cloglog_occ_mean))
  
  range_cells <- which(terra::values(continental_grid_scale2) %in% covar_values_scaled_Joint$grid_cell)
  other_cells <- which(!terra::values(continental_grid_scale2) %in% covar_values_scaled_Joint$grid_cell &
                         !is.na(terra::values(occ_mean_raster)))
  
  occ_mean_raster <- continental_grid_scale2
  terra::values(occ_mean_raster)[!is.na(terra::values(occ_mean_raster))] <- 0
  terra::values(occ_mean_raster)[range_cells] <- icloglog(covar_values_scaled_Joint$cloglog_occ_mean)
  
  occ_intensity_raster <- continental_grid_scale2
  terra::values(occ_intensity_raster) <- NA
  terra::values(occ_intensity_raster)[range_cells] <- covar_values_scaled_Joint$cloglog_occ_mean
  # Blank out parts of the map where occupancy prob. is less than 5%:
  terra::values(occ_intensity_raster)[terra::values(occ_intensity_raster) < cloglog(0.05)] <- NA
  
  occ_sd_raster <- continental_grid_scale2
  terra::values(occ_sd_raster)[!is.na(terra::values(occ_mean_raster))] <- 0
  terra::values(occ_sd_raster)[range_cells] <- covar_values_scaled_Joint$cloglog_occ_sd
  
  p3 <- ggplot() +
    geom_spatraster(data = occ_intensity_raster, maxcell = 1e8,
                    na.rm = T) +
    scale_fill_viridis_c("Intensity", na.value = NA) +
    ggtitle(paste0("Log 'Abundance'")) +
    theme_minimal()
  p1 <- ggplot() +
    geom_spatraster(data = occ_mean_raster, maxcell = 1e8,
                    na.rm = T) +
    scale_fill_viridis_c("Occu.", na.value = NA, limits = c(0, 1)) +
    ggtitle(paste0("Occupancy")) +
    theme_minimal()
  p2 <- ggplot() +
    geom_spatraster(data = occ_sd_raster, maxcell = 1e8,
                    na.rm = T) +
    scale_fill_viridis_c("SD", na.value = NA) +
    ggtitle(paste0("Cloglog-scale SD")) +
    theme_minimal()
  
  p_arr <- gridExtra::arrangeGrob(p3, p1, p2, nrow = 1, 
                                  top = paste0(species, " - joint model results"))
  ggsave(paste0("plots/maps/", species, "-CAR.jpg"), p_arr,
         width=14, height = 6)
}


plot_spatpred_iSDM <- function(res, species) {
  if (!exists("grid_translator_wspec")) source("main_code_nolineage/main_data_prep.R")
  covar_values_scaled_Joint <- generate_spatpred_iSDM(res, species) %>% 
    filter(!is.na(cloglog_occ_mean))
  
  occ_mean_raster <- continental_grid_scale2
  range_cells <- which(terra::values(continental_grid_scale2) %in% covar_values_scaled_Joint$grid_cell)
  other_cells <- which(!terra::values(continental_grid_scale2) %in% covar_values_scaled_Joint$grid_cell &
                         !is.na(terra::values(occ_mean_raster)))
  
  
  terra::values(occ_mean_raster)[!is.na(terra::values(occ_mean_raster))] <- 0
  terra::values(occ_mean_raster)[range_cells] <- icloglog(covar_values_scaled_Joint$cloglog_occ_mean)
  
  p1 <- ggplot() +
    geom_spatraster(data = occ_mean_raster, maxcell = 1e8,
                    na.rm = T) +
    scale_fill_viridis_c("Occu.", na.value = NA, limits = c(0, 1)) +
    ggtitle(paste0("Occupancy"))
  
  sd_raster <- occ_mean_raster
  terra::values(sd_raster)[!is.na(terra::values(continental_grid_scale2))] <- 0
  terra::values(sd_raster)[range_cells] <- covar_values_scaled_Joint$cloglog_occ_sd
  
  # terra::values(sd_raster)[terra::values(sd_raster) > sd_UB] <- sd_UB
  
  p2 <- ggplot() +
    geom_spatraster(data = sd_raster, maxcell = 1e8,
                    na.rm = T) +
    scale_fill_viridis_c("SD", na.value = NA) +
    ggtitle(paste0("Uncertainty in log intensity"))
  
  p_arr <- gridExtra::arrangeGrob(p1, p2, nrow = 1, top = species)
  
  ggsave(paste0("plots/isdm_maps_nolineage/", species, ".jpg"), p_arr, width = 12, height = 6)
}

plot_spatpred_all <- function(result3, result4, species, suffix, 
                              write.out = TRUE, thin = 50, always_use_lambda_int = F) {
  
  covar_values_scaled_Joint <- generate_spatpred(result3, species, inat_betas = F, thin = thin, CIs = F,
                                                 always_use_lambda_int = always_use_lambda_int) %>% 
    filter(!is.na(cloglog_occ_mean))
  covar_values_scaled_CT    <- generate_spatpred(result4, species, thin = thin, CIs = F,
                                                 always_use_lambda_int = always_use_lambda_int) %>% 
    filter(!is.na(cloglog_occ_mean))
  covar_values_scaled_iNat  <- generate_spatpred(result4, species, inat_betas = T, thin = thin, CIs = F,
                                                 always_use_lambda_int = always_use_lambda_int) %>% 
    filter(!is.na(cloglog_occ_mean))
  
  stopifnot(
    nrow(covar_values_scaled_Joint) == nrow(covar_values_scaled_CT)
  )
  stopifnot(
    nrow(covar_values_scaled_Joint) == nrow(covar_values_scaled_iNat)
  )
  
  occ_mean_raster <- continental_grid_scale2
  range_cells <- which(terra::values(continental_grid_scale2) %in% covar_values_scaled_Joint$grid_cell)
  other_cells <- which(!terra::values(continental_grid_scale2) %in% covar_values_scaled_Joint$grid_cell &
                         !is.na(terra::values(occ_mean_raster)))
  
  cv_list <- list(covar_values_scaled_Joint, covar_values_scaled_CT, 
                  covar_values_scaled_iNat)
  sd_min <- cv_list %>% lapply(function(x) x$cloglog_occ_sd) %>% 
    unlist() %>% as.numeric() %>% min()
  sd_UB <- cv_list %>% lapply(function(x) x$cloglog_occ_sd) %>% 
    unlist() %>% as.numeric() %>% quantile(prob = 0.975)
  
  names <- c("Joint", "CT only", "iNat only")
  p1_list <- list()
  p2_list <- list()
  # p3_list <- list()
  for (i in 1:3) {
    
    #### Main plots
    terra::values(occ_mean_raster)[!is.na(terra::values(occ_mean_raster))] <- 0
    terra::values(occ_mean_raster)[range_cells] <- icloglog(cv_list[[i]]$cloglog_occ_mean)
    
    p1_list[[i]] <- ggplot() +
      geom_spatraster(data = occ_mean_raster, maxcell = 1e8,
                      na.rm = T) +
      scale_fill_viridis_c("Occu.", na.value = NA, limits = c(0, 1)) +
      ggtitle(paste0("Occupancy (", names[i], ")"))
    
    ### Plot of uncertainty
    
    sd_raster <- continental_grid_scale2
    
    terra::values(sd_raster)[!is.na(terra::values(continental_grid_scale2))] <- 0
    terra::values(sd_raster)[range_cells] <- cv_list[[i]]$cloglog_occ_sd
    
    terra::values(sd_raster)[terra::values(sd_raster) > sd_UB] <- sd_UB
    
    p2_list[[i]] <- ggplot() +
      geom_spatraster(data = sd_raster, maxcell = 1e8,
                      na.rm = T) +
      scale_fill_viridis_c("SD", na.value = NA, limits = c(sd_min, sd_UB)) +
      ggtitle(paste0("Uncertainty in log intensity (", names[i], ")"))
    
    # 
    # CI_size_raster <- continental_grid_scale2
    # 
    # terra::values(CI_size_raster)[!is.na(terra::values(continental_grid_scale2))] <- 0
    # terra::values(CI_size_raster)[range_cells] <- 
    #   icloglog(cv_list[[i]]$cloglog_occ_UB) - icloglog(cv_list[[i]]$cloglog_occ_LB)
    # 
    # p3_list[[i]] <- ggplot() +
    #   geom_spatraster(data = CI_size_raster, maxcell = 1e8,
    #                   na.rm = T) +
    #   scale_fill_viridis_c("SD", na.value = NA) +
    #   ggtitle(paste0("Width of occupancy 95%CI (", names[i], ")"))
    
  }
  
  
  #### Diff. plots
  
  occ_mean_CTmJoint <- icloglog(covar_values_scaled_CT$cloglog_occ_mean) - 
    icloglog(covar_values_scaled_Joint$cloglog_occ_mean)
  occ_mean_iNatmJoint <- icloglog(covar_values_scaled_iNat$cloglog_occ_mean) - 
    icloglog(covar_values_scaled_Joint$cloglog_occ_mean)
  
  occ_raster_CTmJoint <- continental_grid_scale2
  occ_raster_iNatmJoint <- continental_grid_scale2
  
  terra::values(occ_raster_CTmJoint)[!is.na(terra::values(occ_raster_CTmJoint))] <- 0
  terra::values(occ_raster_CTmJoint)[range_cells] <- occ_mean_CTmJoint
  terra::values(occ_raster_iNatmJoint)[!is.na(terra::values(occ_raster_iNatmJoint))] <- 0
  terra::values(occ_raster_iNatmJoint)[range_cells] <- occ_mean_iNatmJoint
  
  
  diff_p1 <- ggplot() +
    geom_spatraster(data = occ_raster_iNatmJoint, maxcell = 1e8,
                    na.rm = T) +
    ggtitle("Diff. btw. iNat and joint") +
    scale_fill_gradient2("dOccu", mid = "lightgray", na.value = NA) +
    theme_minimal()
  diff_p2 <- ggplot() +
    geom_spatraster(data = occ_raster_CTmJoint, maxcell = 1e8,
                    na.rm = T) +
    ggtitle("Diff. btw. CT and joint") +
    scale_fill_gradient2("dOccu", mid = "lightgray", na.value = NA) +
    theme_minimal()
  
  #### Plots of uncertainty diffs
  
  sd_CTmJoint <- covar_values_scaled_CT$cloglog_occ_sd - 
    covar_values_scaled_Joint$cloglog_occ_sd
  sd_iNatmJoint <- covar_values_scaled_iNat$cloglog_occ_sd - 
    covar_values_scaled_Joint$cloglog_occ_sd
  
  sd_raster_CTmJoint <- continental_grid_scale2
  sd_raster_iNatmJoint <- continental_grid_scale2
  terra::values(sd_raster_CTmJoint)[!is.na(terra::values(sd_raster_CTmJoint))] <- 0
  terra::values(sd_raster_CTmJoint)[range_cells] <- sd_CTmJoint
  terra::values(sd_raster_iNatmJoint)[!is.na(terra::values(sd_raster_iNatmJoint))] <- 0
  terra::values(sd_raster_iNatmJoint)[range_cells] <- sd_iNatmJoint
  
  
  
  diff_p3 <- ggplot() +
    geom_spatraster(data = sd_raster_iNatmJoint, maxcell = 1e8,
                    na.rm = T) +
    ggtitle("Diff. in cloglog-scale SD (iNat-only minus joint)") +
    scale_fill_gradient2("d SD", mid = "lightgray", na.value = NA) +
    theme_minimal()
  diff_p4 <- ggplot() +
    geom_spatraster(data = sd_raster_CTmJoint, maxcell = 1e8,
                    na.rm = T) +
    ggtitle("Diff. in cloglog-scale SD (CT minus joint)") +
    scale_fill_gradient2("d SD", mid = "lightgray", na.value = NA) +
    theme_minimal()
  
  # #### Plots of CI size diffs
  # 
  # CI_mag_CTmJoint <- 
  #   (icloglog(covar_values_scaled_CT$cloglog_occ_UB) - 
  #      icloglog(covar_values_scaled_CT$cloglog_occ_LB)) - 
  #   (icloglog(covar_values_scaled_Joint$cloglog_occ_UB) - 
  #      icloglog(covar_values_scaled_Joint$cloglog_occ_LB))
  # CI_mag_iNatmJoint <- 
  #   (icloglog(covar_values_scaled_iNat$cloglog_occ_UB) - 
  #      icloglog(covar_values_scaled_iNat$cloglog_occ_LB)) - 
  #   (icloglog(covar_values_scaled_Joint$cloglog_occ_UB) - 
  #      icloglog(covar_values_scaled_Joint$cloglog_occ_LB))
  # 
  # CImag_raster_CTmJoint <- continental_grid_scale2
  # CImag_raster_iNatmJoint <- continental_grid_scale2
  # terra::values(CImag_raster_CTmJoint)[!is.na(terra::values(sd_raster_CTmJoint))] <- 0
  # terra::values(CImag_raster_CTmJoint)[range_cells] <- CI_mag_CTmJoint
  # terra::values(CImag_raster_iNatmJoint)[!is.na(terra::values(sd_raster_iNatmJoint))] <- 0
  # terra::values(CImag_raster_iNatmJoint)[range_cells] <- CI_mag_iNatmJoint
  
  # diff_p5 <- ggplot() +
  #   geom_spatraster(data = CImag_raster_iNatmJoint, maxcell = 1e8,
  #                   na.rm = T) +
  #   ggtitle("Diff. in size of 95%CI of occu. (iNat-only minus joint)\nif >0, joint is more confident") +
  #   scale_fill_gradient2("d width", mid = "lightgray", na.value = NA) +
  #   theme_minimal()
  # diff_p6 <- ggplot() +
  #   geom_spatraster(data = CImag_raster_CTmJoint, maxcell = 1e8,
  #                   na.rm = T) +
  #   ggtitle("Diff. in size of 95%CI of occu. (CT minus joint)\nif >0, joint is more confident") +
  #   scale_fill_gradient2("d width", mid = "lightgray", na.value = NA) +
  #   theme_minimal()
  
  
  compareplot_mean <- arrangeGrob(grobs = p1_list, nrow = 2, top = species)
  compareplot_sd   <- arrangeGrob(grobs = p2_list, nrow = 2, top = species)
  # compareplot_CI   <- arrangeGrob(grobs = p3_list, nrow = 1, top = species)
  
  diffplot <- arrangeGrob(diff_p1, diff_p3,# diff_p5,
                          diff_p2, diff_p4,# diff_p6, 
                          nrow = 2, top = species)
  if (write.out) {
    # browser()
    ggsave(plot = compareplot_mean, paste0("plots/multiplot_mean_", species, suffix, ".jpg"),
           width = 12, height = 10)
    ggsave(plot = compareplot_sd, paste0("plots/multiplot_SD_", species, suffix, ".jpg"),
           width = 12, height = 10)
    # ggsave(plot = compareplot_CI, paste0("plots/multiplot_CIsize_", species, suffix, ".jpg"),
    #        width = 12, height = 5)
    ggsave(plot = diffplot, paste0("plots/diffplot_", species, suffix, ".jpg"),
           width = 15, height = 12)
  }
  return(list(
    compareplot_mean = compareplot_mean,
    compareplot_sd = compareplot_sd,
    # compareplot_CI = compareplot_CI,
    diffplot = diffplot
  ))
}


# Function to extract the second number from a string
extract_second_number <- function(input_string) {
  # Find all numbers in the string
  numbers <- str_extract_all(input_string, "\\d+\\.?\\d*")[[1]]
  
  # Check if there are at least two numbers
  if (length(numbers) >= 2) {
    # Parse and return the second number
    return(parse_number(numbers[2]))
  } else {
    return(NA)  # Return NA if there are not enough numbers
  }
}

plot_SVC_relimportance <- function(result, species, suffix, 
                                   write.out = TRUE, thin = 10) {
  
  occ_covars <- gsub("_scaled", "", result$occ_covars)
  
  if (!exists("grid_translator_wspec")) source("main_code/main_data_prep.R")
  
  samples_names <- colnames(result$samples_list[[1]]$samples$chain1)
  samples2_names <- colnames(result$samples_list[[1]]$samples2$chain1)
  
  cell_inds <- which(grepl(paste0("^lambda_beta"), 
                           samples2_names))
  spatcell <- parse_number(samples2_names[cell_inds])
  param_ind <- rep(1:length(occ_covars), each = max(spatcell))
  
  samples2 <- result$samples_list %>% 
    lapply(function(x) {
      
      mtx <- as.data.frame(as.matrix(x$samples2))
      rows <- thin * 1:(nrow(mtx)/thin)
      mtx[rows, cell_inds]
    }) %>% 
    bind_rows()
  
  # For each cell, calculate which variable is most important the most
  importance_df <- samples2 %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(scale4_grid_ID = spatcell, 
           parnum = param_ind) %>% 
    pivot_longer(cols = paste0("V", 1:nrow(samples2))) %>% 
    group_by(scale4_grid_ID, name) %>% 
    summarize(maxvar = parnum[abs(value) == max(abs(value))],
              .groups = "drop") %>% 
    count(scale4_grid_ID, maxvar) %>% 
    group_by(scale4_grid_ID) %>% 
    mutate(pct = n / sum(n)) %>% 
    arrange(-pct) %>% 
    mutate(rank = row_number()) %>% 
    filter(rank == 1)
  
  
  # calculate means and sds of effects
  colnames(grid_translator_wspec)[
    colnames(grid_translator_wspec) == paste0("GRID_ID_", species)
  ] <- "spec_range_gridID"
  
  cells_to_keep <- grid_translator_wspec %>% 
    distinct(scale4_grid_ID, spec_range_gridID) %>% 
    filter(!is.na(spec_range_gridID))
  
  # get the spatial locations
  parname_df <- data.frame(
    parname = occ_covars,
    parnum = 1:length(occ_covars)
  )
  points <- as.data.frame(continental_grid_scale4, xy = T) %>% 
    rename(scale4_grid_ID = lyr.1) %>% 
    left_join(importance_df, by = "scale4_grid_ID") %>% 
    rename(parnum = maxvar) %>% 
    left_join(parname_df, by = "parnum") %>% 
    mutate(parname = ifelse(scale4_grid_ID %in% cells_to_keep$scale4_grid_ID,
                            parname, "Out of range"))
  
  
  mainplot <- ggplot(points) + 
    geom_tile(aes(x, y, fill = parname)) +
    scale_fill_manual("Covariate", values = var_colors) +
    theme_void() + coord_fixed(1)
  
  points2 <- points %>% 
    mutate(parname = ifelse(parname != "Out of range" & pct < 0.5,
                            "Competitive", parname))
  
  mainplot2 <- ggplot(points2) + 
    geom_tile(aes(x, y, fill = parname)) +
    scale_fill_manual("Covariate", values = var_colors) +
    theme_void() + coord_fixed(1)
  
  arr <- gridExtra::arrangeGrob(mainplot, mainplot2, nrow = 1,
                                top = paste0("Which variables are most important for ", species))
  
  if (write.out) {
    ggsave(plot = arr, paste0("plots/relative_importance/topvars_noLineage_", species, suffix, ".jpg"),
           width = 10, height = 4, dpi = 600)
  }
  
  
  overall_imp_df <- samples2 %>%
    t() %>%
    as.data.frame() %>%
    mutate(scale4_grid_ID = spatcell,
           parnum = param_ind) %>%
    pivot_longer(cols = paste0("V", 1:nrow(samples2))) %>%
    group_by(parnum, scale4_grid_ID) %>%
    summarize(value = mean(value)) %>% 
    group_by(parnum) %>% 
    summarize(mean_abs_beta = mean(abs(value)),
              abs_beta_025 = quantile(abs(value), probs = 0.025),
              abs_beta_975 = quantile(abs(value), probs = 0.975),
              mean_beta = mean(value)) %>% 
    left_join(parname_df, by = "parnum")
  
  overallplot <- overall_imp_df %>% 
    ggplot(aes(parname, mean_abs_beta, fill = mean_beta)) +
    geom_col(col = "black") +
    geom_pointrange(aes(parname, mean_abs_beta,
                        ymin = abs_beta_025, 
                        ymax = abs_beta_975)) +
    ylab("Overall absolute effect (relative importance)\nand variability across space") + 
    xlab("") +
    scale_fill_gradient2("Mean effect\n(directional)") +
    coord_flip() + theme_minimal()
  
  if (write.out) {
    ggsave(plot = overallplot, paste0("plots/relative_importance/overall_noLineage_", species, suffix, ".jpg"),
           width = 10, height = 4, dpi = 600)
  }
  
}

# use sig_only = TRUE to gray out covar effects that overlap zero
plot_SVC_effect <- function(result, param, species, suffix, 
                            write.out = TRUE, thin = 10) {
  
  occ_covars <- gsub("_scaled", "", result$occ_covars)
  
  # What's the index of the param we asked for?
  stopifnot(param %in% occ_covars)
  param_ind <- which(occ_covars == param)
  
  if (!exists("grid_translator_wspec")) source("main_code/main_data_prep.R")
  
  samples_names <- colnames(result$samples_list[[1]]$samples$chain1)
  samples2_names <- colnames(result$samples_list[[1]]$samples2$chain1)
  
  cell_inds <- which(grepl(paste0("^lambda_beta\\[\\d*, ", param_ind, "\\]$"), 
                           samples2_names))
  spatcell <- parse_number(samples2_names[cell_inds])
  
  samples2 <- result$samples_list %>% 
    lapply(function(x) {
      
      mtx <- as.data.frame(as.matrix(x$samples2))
      rows <- thin * 1:(nrow(mtx)/thin)
      mtx[rows, cell_inds]
      
    }) %>% 
    bind_rows()
  
  # calculate means and sds of effects
  effs <- data.frame(
    scale4_grid_ID = spatcell,
    mean = colMeans(samples2),
    sd = apply(samples2, 2, sd)
  ) %>% 
    mutate(
      nonzero  = sign(mean + 1.96*sd) == sign(mean - 1.96*sd)
    )
  
  colnames(grid_translator_wspec)[
    colnames(grid_translator_wspec) == paste0("GRID_ID_", species)
  ] <- "spec_range_gridID"
  
  cells_to_keep <- grid_translator_wspec %>% 
    distinct(scale4_grid_ID, spec_range_gridID) %>% 
    filter(!is.na(spec_range_gridID))
  
  # Get the scale4 cell IDs that we need here
  effs_wcell <- effs %>% 
    # select(scale4_grid_ID, spec_range_gridID) %>% 
    # filter(!is.na(spec_range_gridID)) %>% 
    # distinct() %>% 
    # left_join(effs, by = "spec_range_gridID") %>% 
    filter(scale4_grid_ID %in% cells_to_keep$scale4_grid_ID) %>% 
    arrange(scale4_grid_ID)
  
  # browser()
  
  # We expect to have an effect for each row of this_grid_transl
  
  # Make plot with significance filtering
  thisras <- continental_grid_scale4
  values(thisras)[!is.na(values(continental_grid_scale4))] <- 0
  
  values(thisras)[values(continental_grid_scale4) %in% effs_wcell$scale4_grid_ID] <-
    effs_wcell$mean * as.numeric(effs_wcell$nonzero)
  
  mainplot2 <- ggplot() + 
    geom_spatraster(data = thisras) +
    scale_fill_gradient2("Mean eff.\n(Filtered)")
  
  
  # Make plot WITHOUT significance filtering
  thisras <- continental_grid_scale4
  values(thisras)[!is.na(values(continental_grid_scale4))] <- 0
  
  values(thisras)[values(continental_grid_scale4) %in% effs_wcell$scale4_grid_ID] <-
    effs_wcell$mean
  mainplot1 <- ggplot() + 
    geom_spatraster(data = thisras) +
    scale_fill_gradient2("Mean eff.\n(All)")
  
  
  # sdras <- thisras
  # values(sdras)[values(continental_grid_scale4) %in% effs_wcell$scale4_grid_ID] <-
  #   effs_wcell$sd
  
  
  
  # mainplot <- ggplot() + 
  #   geom_spatraster(data = thisras) +
  #   scale_fill_gradient2("Mean eff.")
  # sdplot <- ggplot() + 
  #   geom_spatraster(data = sdras) +
  #   scale_fill_viridis_c("SD")
  
  arr <- gridExtra::arrangeGrob(mainplot1, mainplot2, nrow = 1,
                                top = paste0("SVC effect of ", param, " on ", species))
  
  if (write.out) {
    ggsave(plot = arr, paste0("plots/svcplots/", species, "_", param, suffix, ".jpg"),
           width = 8, height = 3.6)
  }
  
  return(arr)
}

plot_lineage_effect <- function(result, param, species, suffix, 
                                write.out = TRUE, thin = 10) {
  
  occ_covars <- gsub("_scaled", "", result$occ_covars)
  
  lineage_transl <- result$lineage_transl # Need to actually write this out
  
  # What's the index of the param we asked for?
  stopifnot(param %in% occ_covars)
  param_ind <- which(occ_covars == param)
  
  if (!exists("grid_translator_wspec")) source("main_code/main_data_prep.R")
  
  samples_names <- colnames(result$samples_list[[1]]$samples$chain1)
  samples2_names <- colnames(result$samples_list[[1]]$samples2$chain1)
  
  maineff_ind <- which(samples_names == paste0("lambda_beta_mean[", param_ind, "]"))
  
  nLineageEffs <- nrow(lineage_transl) - 1
  lin_int_nums <- which(rep(1:length(occ_covars), each = nLineageEffs) == param_ind)
  
  interaction_inds <- which(
    samples_names %in% paste0("lineage_interactions[", lin_int_nums, "]")
  )
  
  eff_samples <- result$samples_list %>% 
    lapply(function(x) {
      mtx <- as.data.frame(as.matrix(x$samples))
      rows <- thin * 1:(nrow(mtx)/thin)
      mtx[rows, c(maineff_ind, interaction_inds)]
    }) %>% 
    bind_rows()
  
  
  # calculate means and sds of effects
  effs <- lineage_transl %>% mutate(
    mean = NA,
    sd = NA
  )
  
  for (i in 1:nrow(effs)) {
    if (i == nrow(effs)) {
      # For lineage n, beta is b0 - sum(b_i)
      beta_vec <- eff_samples[,1] - 
        rowSums(cbind(rep(0, nrow(eff_samples)), eff_samples[, 2:nrow(effs)]))
    } else {
      # For lineages 1 to n-1, beta is b0 + b_i
      beta_vec <- eff_samples[,1] + eff_samples[, i+1]
    }
    
    effs$mean[i] <- mean(beta_vec)
    effs$sd[i] <- sd(beta_vec)
  }
  
  colnames(grid_translator_wspec)[
    colnames(grid_translator_wspec) == paste0("GRID_ID_", species)
  ] <- "spec_range_gridID"
  
  colnames(lineage_df)[colnames(lineage_df) == species] <- "target"
  this_lineage_df <- lineage_df %>% 
    filter(!is.na(target))
  
  # Populate the rasters 
  lineage_ras <- continental_grid_scale4
  values(lineage_ras)[!is.na(values(continental_grid_scale4))] <- 0
  mean_ras <- sd_ras <- lineage_ras
  
  values(lineage_ras)[values(continental_grid_scale4) %in% 
                        this_lineage_df$scale4_grid_ID] <- this_lineage_df$target
  
  for (i in 1:nrow(effs)) {
    values(mean_ras)[values(lineage_ras) == effs$target_lineage[i]] <-
      effs$mean[i]
    values(sd_ras)[values(lineage_ras) == effs$target_lineage[i]] <-
      effs$sd[i]
  }
  
  
  lineageplot <- ggplot() + 
    geom_spatraster(data = lineage_ras) +
    scale_fill_viridis_c(option = "A", "Lineage")
  mainplot <- ggplot() + 
    geom_spatraster(data = mean_ras) +
    scale_fill_gradient2("Mean eff.")
  sdplot <- ggplot() + 
    geom_spatraster(data = sd_ras) +
    scale_fill_viridis_c("SD")
  
  arr <- gridExtra::arrangeGrob(lineageplot, mainplot, sdplot, nrow = 1,
                                top = paste0("SVC effect of ", param, " on ", species))
  
  if (write.out) {
    ggsave(plot = arr, paste0("plots/lineageplots/", species, "_", param, suffix, ".jpg"),
           width = 11, height = 3.6)
  }
  
  return(arr)
}




#### Helper functions ####
get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

get_hypothesis_summary_expanded <- function(temp) {
  
  samples_mtx <- temp$samples_list[[1]]$samples %>% 
    lapply(as.matrix) %>% 
    do.call(what = rbind) %>% 
    as.data.frame() %>% 
    select(all_of(indic_vars))
  vars <- colnames(samples_mtx)
  
  frac_hypothesis <- matrix(NA, nrow = length(occ_covars), 
                            ncol = 8)
  
  for (i in 1:nrow(frac_hypothesis)) {
    this_indvar_CAR <- which(vars == paste0("CAR_gamma[", i, "]"))
    this_indvar_eco <- which(vars == paste0("eco_gamma[", i, "]"))
    
    frac_hypothesis[i, 1] <- mean(samples_mtx[, this_indvar_CAR] == 0 & 
                                    samples_mtx[, this_indvar_eco] == 0)
    frac_hypothesis[i, 2] <- mean(samples_mtx[, this_indvar_CAR] == 0 & 
                                    samples_mtx[, this_indvar_eco] == 1)
    frac_hypothesis[i, 3] <- mean(samples_mtx[, this_indvar_CAR] == 1 & 
                                    samples_mtx[, this_indvar_eco] == 0)
    frac_hypothesis[i, 4] <- mean(samples_mtx[, this_indvar_CAR] == 1 & 
                                    samples_mtx[, this_indvar_eco] == 1)
    
  }
  
  return(frac_hypothesis)
}

get_hypothesis_summary_concise <- function(temp) {
  
  samples_mtx <- temp$samples_list[[1]]$samples %>% 
    lapply(as.matrix) %>% 
    do.call(what = rbind) %>% 
    as.data.frame() %>% 
    select(all_of(indic_vars))
  vars <- colnames(samples_mtx)
  
  frac_hypothesis <- matrix(NA, nrow = length(occ_covars), 
                            ncol = 2)
  
  for (i in 1:nrow(frac_hypothesis)) {
    this_indvar_CAR <- which(vars == paste0("CAR_gamma[", i, "]"))
    this_indvar_eco <- which(vars == paste0("eco_gamma[", i, "]"))
    
    frac_hypothesis[i, 1] <- mean(samples_mtx[, this_indvar_eco] == 1)
    frac_hypothesis[i, 2] <- mean(samples_mtx[, this_indvar_CAR] == 1)
    
  }
  
  return(frac_hypothesis)
}



get_var_samples <- function(temp) {
  temp$samples_list[[1]]$samples %>% 
    lapply(as.matrix) %>% 
    do.call(what = rbind) %>% 
    as.data.frame() %>% 
    select(all_of(sigma_vars)) %>% 
    mutate(iter = row_number()) %>% 
    pivot_longer(all_of(sigma_vars)) %>% 
    mutate(parname = occ_covars[parse_number(name)],
           svc_type = substr(name,7,9))
}




#### Define some useful objects ####

svc_type_levels <- c("eco", "CAR")

indic_vars <- c(
  paste0("CAR_gamma[", 1:length(occ_covars), "]"),
  paste0("eco_gamma[", 1:length(occ_covars), "]")
)

sigma_vars <- c(
  paste0("sigma_CAR[", 1:length(occ_covars), "]"),
  paste0("sigma_eco[", 1:length(occ_covars), "]")
)

supported_type_colors <- c(
  "None" = "#dddddd",
  "CAR" = "#e41a1c",
  "Ecoregion" = "#377eb8",
  "Lineage" = "#4daf4a"
)

