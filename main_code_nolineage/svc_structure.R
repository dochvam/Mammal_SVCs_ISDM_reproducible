###############################################################################
#'
#' main_code_nolineage/svc_structure.R
#' Author: Ben R. Goldstein
#' Description: Code for evaluating whether SVCs are largely explained by
#'    nonlinear species-environment relationships.
#'
###############################################################################

#### First do non-lineage models ####

source('setup_code/define_fundamentals_nolineage.R')
source("main_code_nolineage/vis_fn.R")
source("main_code_nolineage/main_data_prep.R")

library(mgcv)
library(nimble)
redo <- F


occ_covars_clean <- gsub("_scaled", "", occ_covars)
parname_df <- data.frame(
  parname = occ_covars_clean,
  parnum = 1:length(occ_covars)
)


# Retrieve all the lambda_betas

lambda_beta_df_list <- list()

for (i in 1:nrow(taxon_key)) {
  
  result <- readRDS(paste0("intermediate/integration_results/test_kfold_", 
                           taxon_key$common_name_clean[i], 
                           "_joint_SVCs_main_noLineage_samples.RDS"))
  
  target_col <- paste0("GRID_ID_", taxon_key$common_name_clean[i])
  cells_to_keep <- grid_translator_wspec[!is.na(grid_translator_wspec[, target_col]), ]
  
  this_samples <- as.matrix(result$samples_list[[1]]$samples2)
  this_samples <- this_samples[, grepl("lambda_beta\\[", colnames(this_samples))]
  
  scale4_ID_vec <- parse_number(colnames(this_samples))
  parnum <- rep(1:7, each = length(unique(scale4_ID_vec)))
  
  this_samples <- this_samples[, scale4_ID_vec %in% cells_to_keep$scale4_grid_ID]
  
  parnum <- parnum[scale4_ID_vec %in% cells_to_keep$scale4_grid_ID]
  scale4_ID_vec <- scale4_ID_vec[scale4_ID_vec %in% cells_to_keep$scale4_grid_ID]
  
  lambda_beta_df_list[[i]] <- data.frame(
    parnum = parnum,
    scale4_grid_ID = scale4_ID_vec,
    mean = colMeans(this_samples),
    se = unlist(apply(this_samples, 2, sd)),
    LB = unlist(apply(this_samples, 2, quantile, probs = 0.025)),
    UB = unlist(apply(this_samples, 2, quantile, probs = 0.975))
  ) %>% 
    mutate(species = taxon_key$common_name_clean[i]) %>%
    left_join(parname_df, by = "parnum")
  
}

lambda_beta_df_nolineage <- bind_rows(lambda_beta_df_list)


#### lineage models ####
source('setup_code/define_fundamentals.R')
source("main_code_lineage/vis_fn.R")
source("main_code_lineage/main_data_prep.R")


lambda_beta_df_nolineage <- lambda_beta_df_nolineage %>% 
  filter(!species %in% taxon_key$common_name_clean)

lambda_beta_df_list_lineage <- list()

for (i in 1:nrow(taxon_key)) {
  
  result <- readRDS(paste0("intermediate/integration_results/test_kfold_", 
                           taxon_key$common_name_clean[i], 
                           "_joint_SVCs_main_samples.RDS"))
  
  target_col <- paste0("GRID_ID_", taxon_key$common_name_clean[i])
  
  this_samples <- as.matrix(result$samples_list[[1]]$samples2)
  this_samples <- this_samples[, grepl("lambda_beta\\[", colnames(this_samples))]
  
  grid_cell <- parse_number(colnames(this_samples))
  parnum <- rep(1:7, each = length(unique(grid_cell)))
  
  grid_transl_temp <- grid_translator_wspec
  colnames(grid_transl_temp)[grepl(taxon_key$common_name_clean[i], 
                                   colnames(grid_transl_temp))] <- "grid_cell"
  grid_transl_temp <- grid_transl_temp[, c("grid_cell", "scale4_grid_ID")] %>% 
    distinct() %>% 
    na.omit()
  
  lambda_beta_df_list_lineage[[i]] <- data.frame(
    parnum = parnum,
    grid_cell = grid_cell,
    mean = colMeans(this_samples),
    se = unlist(apply(this_samples, 2, sd)),
    LB = unlist(apply(this_samples, 2, quantile, probs = 0.025)),
    UB = unlist(apply(this_samples, 2, quantile, probs = 0.975))
  ) %>% 
    mutate(species = taxon_key$common_name_clean[i]) %>%
    left_join(parname_df, by = "parnum") %>% 
    left_join(grid_transl_temp, by = "grid_cell") %>% 
    select(-grid_cell)
}

lambda_beta_df_lineage <- bind_rows(lambda_beta_df_list_lineage)
write_csv(lambda_beta_df_lineage, "intermediate/lambda_beta_summary_LINEAGE.jpg")

lambda_beta_df <- bind_rows(lambda_beta_df_nolineage, lambda_beta_df_lineage)

# Coordinates of each grid cell centroid
cell_coords <- as.data.frame(continental_grid_scale4, xy = TRUE) %>%
  rename(scale4_grid_ID = lyr.1)

# Average covar. values by S2 grid cell
covar_values_byS2 <- covar_values %>%
  rename(scale2_grid_ID = grid_cell) %>%
  left_join(grid_translator[, c("scale2_grid_ID", "scale4_grid_ID")],
            by = "scale2_grid_ID") %>%
  select(scale4_grid_ID, all_of(occ_covars_clean)) %>%
  pivot_longer(cols = all_of(occ_covars_clean)) %>%
  group_by(scale4_grid_ID, parname = name) %>%
  summarize(value = mean(value, na.rm = TRUE)) %>% 
  left_join(cell_coords, by = "scale4_grid_ID")


# Test plot
covar_values_byS2 %>%
  filter(parname == "EVI_mean", !is.na(value)) %>%
  ggplot() +
  geom_tile(aes(x, y, fill = value)) +
  theme_minimal() +
  scale_fill_viridis_c()


source('setup_code/define_fundamentals_nolineage.R')

# Loop over species. Is temp. associated w. temp effect?
plist <- list()
gam_pred_df_list <- list()
gam_summary_df_list <- list()
pb <- progress::progress_bar$new(total = nrow(taxon_key) * length(occ_covars_clean))
ct <- 0
for (i in 1:nrow(taxon_key)) {
  for (j in 1:length(occ_covars_clean)) {
    pb$tick(); ct <- ct + 1
    
    this_compare <- lambda_beta_df %>%
      filter(species == taxon_key$common_name_clean[i]) %>%
      filter(parname == occ_covars_clean[j]) %>%
      left_join(covar_values_byS2, by = c("scale4_grid_ID", "parname")) %>%
      mutate(value_scaled = as.numeric(scale(value))) %>%
      filter(!is.na(value_scaled))
    
    plist[[ct]] <- ggplot(this_compare, aes(value, mean)) + 
      geom_point() +
      geom_smooth() +
      theme_minimal() + 
      ggtitle(paste0(taxon_key$common_name_clean[i], " - ",
                     occ_covars_clean[j])) #+ xlab("Max temp (10ths of degree C)") + ylab("Mean effect")
    
    ### Fit spline model
    
    gam_fit <- gam(mean ~ value_scaled + 
                     s(value_scaled, bs = "cr",
                       k = min(50, length(unique(this_compare$value_scaled)) - 1)), 
                   weights = 1 / (se^2), data = this_compare)
    s <- summary(gam_fit)
    gam_summary_df_list[[ct]] <- data.frame(
      pval_main = s$p.pv[2],
      pval_spline = s$s.pv,
      rsq = s$r.sq,
      dev.expl = s$dev.expl,
      param = occ_covars_clean[j],
      species = taxon_key$common_name_clean[i]
    )
    
    range <- quantile(this_compare$value_scaled, probs = c(0.1, 0.9))
    pred_df <- data.frame(value_scaled = seq(range[1], range[2], length.out = 20))
    pred_df$slope <- predict(gam_fit, newdata = pred_df)
    pred_df$se <- predict(gam_fit, newdata = pred_df, se.fit = TRUE)$se.fit
    pred_df$species <- taxon_key$common_name_clean[i]
    pred_df$param <- occ_covars[j]
    
    gam_pred_df_list[[ct]] <- pred_df
  }
}

### Summaries from GAMs ###

gam_summary_df <- bind_rows(gam_summary_df_list) %>% 
  filter(species != "american_bison") %>% 
  rename(covariate = param) %>% 
  left_join(covar_type_df)

# Correct via FDR
gam_summary_df$pval_spl_corr <- p.adjust(gam_summary_df$pval_spline, method = "fdr", 
                                     n = nrow(gam_summary_df))
gam_summary_df$pval_main_corr <- p.adjust(gam_summary_df$pval_main, method = "fdr", 
                                     n = nrow(gam_summary_df))
gam_summary_df$sig_main <- gam_summary_df$pval_main_corr < 0.05
gam_summary_df$sig_spl <- gam_summary_df$pval_spl_corr < 0.05

# GAM hypothesis test
table(gam_summary_df$sig_spl)
table(gam_summary_df$sig_main)

gam_summary_df %>% 
  filter(species != "american_bison") %>% 
  ggplot() +
  geom_histogram(aes(dev.expl, fill = sig_spl)) +
  scale_fill_manual("spline pval < 0.05", values = c("#e0e0e0", "#8856a7")) +
  theme_minimal() + 
  xlab("GAM Pct. deviance explained") + ylab("Count of species/covar pairs")


ggplot(gam_summary_df) +
  geom_boxplot(aes(species, dev.expl)) +
  coord_flip()
ggplot(gam_summary_df) +
  geom_boxplot(aes(covariate, dev.expl)) +
  coord_flip()
ggplot(gam_summary_df) +
  geom_boxplot(aes(type, dev.expl), outlier.shape = NA) +
  coord_flip() +
  theme_minimal() +
  xlab("Niche dimension") + ylab("Pct. deviance explained by nonlinearity")

gam_pred_df <- bind_rows(gam_pred_df_list)

gam_pred_df %>%
  filter(se < 5) %>%
  filter(param == "Temp_max_scaled") %>%
  ggplot() +
  geom_ribbon(aes(value_scaled, ymin = slope-1.96*se, ymax = slope+1.96*se, 
                  group = species), alpha = 0.5) +
  geom_line(aes(value_scaled, slope, group = species), alpha = 0.5) +
  # facet_wrap(~species, scales = "free_y") +
  geom_hline(yintercept = 0)




#### Snowshoe hare example plot

this_compare <- lambda_beta_df %>%
  filter(species == "snowshoe_hare") %>%
  filter(parname == "Temp_max") %>%
  left_join(covar_values_byS2, by = c("scale4_grid_ID", "parname")) %>%
  mutate(value_scaled = as.numeric(scale(value))) %>%
  filter(!is.na(value_scaled))

ggplot(this_compare) +
  geom_point(aes(value_scaled, mean))
