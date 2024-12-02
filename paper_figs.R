###############################################################################
#'
#' paper_figs.R
#' Author: Ben R. Goldstein
#' Description: This behemoth of a file processes all the output from 
#'    run_models_lineage.R and run_models_nolineage.R. It produces the main
#'    graphs used in the manuscript and supplement and also produces predictions
#'    of the intensity field for all species across the study area. It depends 
#'    on main_code_nolineage/vis_fn.R for a lot of the most interesting 
#'    operations.
#'
###############################################################################

library(tidyverse)
library(coda)

redo <- FALSE

species_to_exclude <- c()

#### Produce summaries of SVC effects ####

#### Process no-lineage svc betas
if (!dir.exists('summary_final')) dir.create("summary_final")

source("main_code_nolineage/model_fn.R")
source("main_code_nolineage/model_code_main.R")
source("setup_code/define_fundamentals_nolineage.R")

result_files <- list.files("intermediate/integration_results/", full.names = TRUE)
result_files <- result_files[grepl('noLineage', result_files)]
result_files <- result_files[!grepl('NOSVCS', result_files)]

target_species <- taxon_key$common_name_clean

finished <- data.frame(species = target_species, rhat = NA, n.eff = NA)

summary_list <- list()

pb <- progress::progress_bar$new(total = nrow(finished))
for (i in 1:nrow(finished)) {
  pb$tick()
  this_results <- result_files[grep(finished$species[i], result_files)]
  
  temp <- lapply(this_results, function(x) {
    thisres <- readRDS(x)
    thisres$samples_list[[1]]$samples %>% as.list()
  }) %>% 
    do.call(what = c)
  
  temp <- mcmc.list(temp)
  
  summary_list[[i]] <- temp %>% 
    MCMCvis::MCMCsummary() %>% 
    mutate(species = target_species[i])
  summary_list[[i]]$param <- rownames(summary_list[[i]])
  rownames(summary_list[[i]]) <- NULL
  
  write_csv(summary_list[[i]], paste0("summary_final/", finished$species[i], 
                                      "_final_summary1.csv"))
  
  finished$rhat[i] <- summary_list[[i]]$Rhat[summary_list[[i]]$param == "logDens"]
  finished$n.eff[i] <- summary_list[[i]]$n.eff[summary_list[[i]]$param == "logDens"]
  
  this_samples2 <- lapply(this_results, function(x) {
    thisres <- readRDS(x)
    thisres$samples_list[[1]]$samples2 %>% as.list()
  }) %>%
    do.call(what = c)
  
  this_samples2 <- mcmc.list(this_samples2)
  this_summary2 <- this_samples2 %>%
    MCMCvis::MCMCsummary(
      Rhat = F, n.eff = F
    ) %>%
    mutate(species = target_species[i])
  
  this_summary2$param <- rownames(this_summary2)
  rownames(this_summary2) <- NULL
  
  this_summary2 %>% 
    filter(grepl("lambda_beta", param)) %>% 
    write_csv(paste0("lambda_beta_results/", finished$species[i], "_lambda_beta_summary.csv"))
}
write_csv(finished, "nolineage_finished.csv")

#### Process lineage svc betas

source("main_code_lineage/model_fn.R")
source("main_code_lineage/model_code_main.R")
source("setup_code/define_fundamentals.R")

result_files <- list.files("intermediate/integration_results/", full.names = TRUE)
result_files <- result_files[!grepl('noLineage', result_files)]
result_files <- result_files[!grepl('NOSVCS', result_files)]

target_species <- taxon_key$common_name_clean
target_species <- target_species[!target_species %in% c("mexican_flying_squirrel")]

finished_lineage <- data.frame(species = target_species, rhat = NA, n.eff = NA)

summary_list <- list()

pb <- progress::progress_bar$new(total = nrow(finished_lineage))
for (i in 1:nrow(finished_lineage)) {
  pb$tick()
  this_results <- result_files[grep(finished_lineage$species[i], result_files)]
  
  temp <- lapply(this_results, function(x) {
    thisres <- readRDS(x)
    thisres$samples_list[[1]]$samples %>% as.list()
  }) %>% 
    do.call(what = c)
  
  temp <- mcmc.list(temp)

  summary_list[[i]] <- temp %>% 
    MCMCvis::MCMCsummary() %>% 
    mutate(species = target_species[i])
  summary_list[[i]]$param <- rownames(summary_list[[i]])
  rownames(summary_list[[i]]) <- NULL
  
  write_csv(summary_list[[i]], paste0("summary_final/", finished_lineage$species[i], 
                                      "_final_summary1_lineage.csv"))
  
  finished_lineage$rhat[i] <- 
    summary_list[[i]]$Rhat[summary_list[[i]]$param == "logDens"]
  finished_lineage$n.eff[i] <- 
    summary_list[[i]]$n.eff[summary_list[[i]]$param == "logDens"]
  
  this_samples2 <- lapply(this_results, function(x) {
    thisres <- readRDS(x)
    thisres$samples_list[[1]]$samples2 %>% as.list()
  }) %>%
    do.call(what = c)
    
  this_samples2 <- mcmc.list(this_samples2)
  this_summary2 <- this_samples2 %>%
    MCMCvis::MCMCsummary(
      Rhat = F, n.eff = F
    ) %>%
    mutate(species = target_species[i])
  
  this_summary2$param <- rownames(this_summary2)
  rownames(this_summary2) <- NULL
  
  this_summary2 %>% 
    filter(grepl("lambda_beta", param)) %>% 
    mutate(type = "lineage") %>% 
    write_csv(paste0("lambda_beta_results/", finished_lineage$species[i], "_lambda_beta_summary_lineage.csv"))  
}



#### Get inclusion probabilities ####

source("setup_code/define_fundamentals.R")
source("main_code_nolineage/vis_fn.R")

taxon_key_lineage <- taxon_key

summary_threshold <- 0.5
strong_threshold <- 0.9

coarse_summary_lineage <- data.frame(
  species = "NA",
  covariate = "NA", 
  svc_type = "NA",
  supported = NA
)[0, ]


# Hypothesis support from lineage models
target_species <- taxon_key$common_name_clean
target_species <- target_species[!target_species %in% species_to_exclude]
suffix <- "_main"

occ_covars_clean <- occ_covars %>% 
  gsub(pattern = "_scaled", replacement = "")

for (i in 1:length(target_species)) {
  input_file <- paste0("res_combined_final/", target_species[i], "_final_res_lineage.RDS")
  if (file.exists(input_file)) {
    temp <- readRDS(input_file)
    
    frac_hypothesis_concise <- get_hypothesis_summary_concise_lineage(temp)
    colnames(frac_hypothesis_concise) <- c("Lineage", "Ecoregion", "CAR")
    frac_hypothesis_concise <- as.data.frame(frac_hypothesis_concise) %>% 
      mutate(param = occ_covars) %>% 
      mutate(param = gsub("_scaled", "", param)) %>% 
      pivot_longer(cols = all_of(c("Lineage", "Ecoregion", "CAR"))) 
    
    coarse_summary_lineage <- bind_rows(coarse_summary_lineage,
                                        data.frame(
                                          species   = target_species[i],
                                          covariate = frac_hypothesis_concise$param, 
                                          svc_type  = frac_hypothesis_concise$name,
                                          supported = frac_hypothesis_concise$value >= summary_threshold,
                                          strong    = frac_hypothesis_concise$value >= strong_threshold
                                        ))
  }
}

source("setup_code/define_fundamentals_nolineage.R")
source("main_code_nolineage/vis_fn.R")

taxon_key_nolineage <- taxon_key

coarse_summary_nolineage <- data.frame(
  species = "NA",
  covariate = "NA", 
  svc_type = "NA",
  supported = NA
)[0, ]

target_species <- taxon_key$common_name_clean
target_species <- target_species[!target_species %in% species_to_exclude]
suffix <- "_main"

occ_covars_clean <- occ_covars %>% 
  gsub(pattern = "_scaled", replacement = "")

for (i in 1:length(target_species)) {
  input_file <- paste0("integration_results/", target_species[i], "_final_res.RDS")
  
  if (file.exists(input_file)) {
    temp <- readRDS(input_file)
    
    frac_hypothesis_concise <- get_hypothesis_summary_concise(temp)
    colnames(frac_hypothesis_concise) <- c("Ecoregion", "CAR")
    frac_hypothesis_concise <- as.data.frame(frac_hypothesis_concise) %>% 
      mutate(param = occ_covars) %>% 
      mutate(param = gsub("_scaled", "", param)) %>% 
      pivot_longer(cols = all_of(c("Ecoregion", "CAR"))) 
    
    coarse_summary_nolineage <- bind_rows(coarse_summary_nolineage,
                                          data.frame(
                                            species = target_species[i],
                                            covariate = frac_hypothesis_concise$param, 
                                            svc_type = frac_hypothesis_concise$name,
                                            supported = frac_hypothesis_concise$value >= summary_threshold,
                                            strong    = frac_hypothesis_concise$value >= strong_threshold
                                          ))
  }
}

final_results <- bind_rows(
  coarse_summary_nolineage %>% filter(!species %in% unique(coarse_summary_lineage$species)),
  coarse_summary_lineage
)
write_csv(final_results, "intermediate/final_results_inclprob.csv")


#### Main incl prob figure ####
final_results <- read_csv("intermediate/final_results_inclprob.csv") %>% 
  filter(!species %in% species_to_exclude)

# Get the percentages
final_results %>% 
  group_by(svc_type) %>% 
  summarize(support_pct = mean(supported))

spec_cts <- final_results %>% 
  group_by(species) %>% 
  summarize(supported = sum(supported))
mean(spec_cts$supported)

# Combine the results we want to interpret
overall_ct_plot <- final_results %>% 
  mutate(support_type = ifelse(strong, "Supported (Strong)", 
                               ifelse(supported, "Supported", "Not Supported"))) %>% 
  mutate(svc_type = recode(svc_type, "Lineage" = "H1: Lineage",
                           "Ecoregion" = "H2: Ecoregion", "CAR" = "H3: Fine-scale structure")) %>% 
  ggplot() +
  geom_bar(aes(svc_type, fill = support_type), position="fill") +
  scale_x_discrete(limits = rev) +
  scale_fill_manual("", values = c("#e0e0e0", "#9ebcda", "#8856a7")) +
  theme_minimal() + ylab("% species-covariate effect pairs") + xlab("") +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank()) +
  coord_flip()

ggsave("paper_figs/overall_ct_plot.jpg", overall_ct_plot,
       width = 6, height = 1.5, dpi = 600)


#### Results by covariate ####

bycov_ct_plot <- final_results %>% 
  left_join(covar_type_df) %>% 
  mutate(support_type = ifelse(strong, "Supported (Strong)", 
                               ifelse(supported, "Supported", "Not Supported"))) %>% 
  mutate(svc_type = recode(svc_type, 
                           "Lineage" = "H1",
                           "Ecoregion" = "H2",
                           "CAR" = "H3")) %>% 
  ggplot() +
  scale_x_discrete(limits = rev) +  geom_bar(aes(svc_type, fill = support_type), position="fill",
                                             show.legend = FALSE) +
  scale_fill_manual("", values = c("#e0e0e0", "#9ebcda", "#8856a7")) +
  theme_minimal() + ylab("% species-covariate effect pairs") + xlab("") +
  coord_flip() +
  facet_wrap(~type)

ggsave("paper_figs/bycov_ct_plot.jpg", bycov_ct_plot,
       width = 6, height = 4, dpi = 600)

allcounts <- gridExtra::arrangeGrob(overall_ct_plot + ggtitle("(A) Hypothesis rates of support overall"), 
                                    bycov_ct_plot + ggtitle("(B) Hypothesis rates of support by niche dimension"), 
                                    ncol = 1)
ggsave("paper_figs/all_counts.jpg", allcounts,
       width = 7, height = 4, dpi = 600)


# percents
final_results %>% 
  count(svc_type, supported) %>% 
  group_by(svc_type) %>% 
  mutate(pct = n/sum(n))

# SVCs per species
temp <- final_results %>% 
  filter(supported) %>% 
  count(species)
mean(temp$n)


#### SVC covariate plots ####
# these get dumped into plots/svcplots/

source('main_code_nolineage/main_data_prep.R')
source('main_code_nolineage/vis_fn.R')

for (i in 1:length(target_species)) {
  for (j in 1:length(occ_covars_clean)) {
    capture <- plot_SVC_effect(param = occ_covars_clean[j], 
                               species = target_species[i], suffix = suffix,
                               write.out = T)
  }
}

#### Theta1 plot ####

summary_df <- lapply(list.files("summary_final", full.names = TRUE), function(x) {
  ignored <- suppressMessages(temp <- read_csv(x))
  temp$lineage <- grepl("lineage", x)
  temp
}) %>% 
  bind_rows()
theta_df <- summary_df %>% 
  filter(param == "theta1") %>% 
  arrange(-lineage) %>% 
  filter(!duplicated(species)) %>% 
  filter(!species %in% species_to_exclude)

theta_plot <- theta_df %>% 
  mutate(nonzero = sign(`2.5%`) == sign(`97.5%`)) %>% 
  mutate(`2.5%` = pmax(-0.5, `2.5%`), `97.5%` = pmin(1, `97.5%`)) %>% 
  ggplot() +
  geom_pointrange(aes(species, mean, ymin = `2.5%`, ymax = `97.5%`)) +
  coord_flip() +
  theme_minimal() +
  ylim(c(-0.5, 1)) +
  geom_hline(yintercept = 0)
ggsave(theta_plot, filename = "paper_figs/theta1_plot.jpg",
       width = 5, height = 8, dpi = 600)


#### Sign flips figure ####

lambda_beta_files <- data.frame(
  fn = list.files("lambda_beta_results/", full.names = TRUE)
) %>% 
  filter(!species %in% species_to_exclude) %>% 
  mutate(lin = grepl("_lineage", fn)) %>% 
  mutate(fn_clean = gsub("_lineage", "", fn)) %>% 
  arrange(-lin) %>% 
  filter(!duplicated(fn_clean))

lambda_beta_df_all <- lapply(
  lambda_beta_files$fn,
  read_csv, show_col_types = FALSE
) %>% 
  bind_rows() %>% 
  mutate(pos = `2.5%` > 0, neg = `97.5%` < 0,
         parnum = parse_number(substr(param, nchar(param) - 3, nchar(param))),
         parname = occ_covars_clean[parnum]
  )

grid_coords <- as.data.frame(continental_grid_scale4, xy = TRUE) %>% 
  rename(scale4_grid_ID = lyr.1)

lambda_beta_df_all %>% 
  mutate(scale4_grid_ID = parse_number(param)) %>% 
  left_join(grid_coords, by = "scale4_grid_ID") %>% 
  select(x, y, species, parname, mean, `2.5%`, `50%`, `97.5%`) %>% 
  write_csv("svc_estimates.csv")

# Avg. magnitude diffs
qdiff_df <- lambda_beta_df_all %>% 
  group_by(species, parname) %>% 
  summarize(q10 = quantile(mean, prob = 0.1), 
            q90 = quantile(mean, prob = 0.9),
            qdiff = q90 - q10, .groups = "drop")
quantile(qdiff_df$qdiff, probs = c(0.1, 0.5, 0.9))


sign_flip_df <- lambda_beta_df_all %>% 
  group_by(species, parname) %>% 
  summarize(npos = sum(pos), nneg = sum(neg)) %>% 
  mutate(is_sign_flip = npos > 0 & nneg > 0)

mean(sign_flip_df$is_sign_flip)

sign_flip_plot <- sign_flip_df %>% 
  ggplot() + 
  geom_bar(aes(parname, fill = is_sign_flip), position = "stack") +
  theme_minimal() + 
  scale_fill_viridis_d("Sign flip", end = 0.8) +
  coord_flip() + xlab("") + ylab("Num. species/covariate pairs")
ggsave(sign_flip_plot, filename = "paper_figs/sign_flips.jpg",
       width = 5, height = 3, dpi = 600)

#### Weighted mean effects ####
lambda_wmeans_plot <- lambda_beta_df_all %>% 
  mutate(scale4_grid_ID = parse_number(param)) %>% 
  mutate(weight = 1 / (`97.5%` - `2.5%`)) %>% 
  group_by(scale4_grid_ID, parname) %>% 
  summarize(wmean = sum(mean * weight) / sum(weight)) %>% 
  left_join(grid_coords) %>% 
  ggplot() +
  geom_tile(aes(x, y, fill = wmean)) +
  facet_wrap(~parname, nrow = 2) +
  theme_bw() + theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  ylab("") + xlab("") +
  scale_fill_gradient2("Weighted\nmean effect", mid = "#cccccc") +
  coord_fixed()
ggsave("paper_figs/weighted_means.jpg", lambda_wmeans_plot,
       width = 8, height = 4.5, dpi = 600)

#### Nonlinearity ####
library(mgcv)
source("main_code_lineage/main_data_prep.R")
grid_transl_lineage <- grid_translator_wspec
source("main_code_nolineage/main_data_prep.R")
grid_transl_nolineage <- grid_translator_wspec

occ_covars_clean <- gsub("_scaled", "", occ_covars)
parname_df <- data.frame(
  parname = occ_covars_clean,
  parnum = 1:length(occ_covars)
)
lambda_beta_df <- lambda_beta_df_all %>% 
  mutate(scale4_grid_ID = parse_number(param))

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

# Loop over species. Is temp. associated w. temp effect?
gam_pred_df_list <- list()
gam_summary_df_list <- list()
pb <- progress::progress_bar$new(total = nrow(taxon_key) * length(occ_covars_clean))
ct <- 0
taxon_key <- taxon_key %>% filter(!common_name_clean %in% species_to_exclude)

for (i in 1:nrow(taxon_key)) {
  for (j in 1:length(occ_covars_clean)) {
    pb$tick()
    
    ct <- ct + 1
    this_compare <- lambda_beta_df %>%
      filter(species == taxon_key$common_name_clean[i]) %>%
      filter(parname == occ_covars_clean[j]) %>%
      left_join(covar_values_byS2, by = c("scale4_grid_ID", "parname")) %>%
      mutate(value_scaled = as.numeric(scale(value))) %>%
      filter(!is.na(value_scaled), !is.na(scale4_grid_ID))
    
    ### Fit spline model
    gam_fit <- gam(mean ~ value_scaled + 
                     s(value_scaled, bs = "cr",
                       k = min(50, length(unique(this_compare$value_scaled)) - 1)), 
                   weights = 1 / (sd^2), data = this_compare)
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

gam_summary_df <- bind_rows(gam_summary_df_list) %>% 
  rename(covariate = param) %>% 
  left_join(covar_type_df)

incl_results <- final_results %>% 
  filter(supported)

nonlin_fig_C <- gam_summary_df %>% 
  filter(paste0(species, covariate) %in% paste0(incl_results$species, incl_results$covariate)) %>%
  ggplot() +
  geom_boxplot(aes(type, rsq), outlier.shape = NA) +
  coord_flip() +
  ylim(c(-0.01, 0.3)) +
  theme_minimal() +
  xlab("Niche dimension") + ylab("R-squared of nonlinearity GAM") +
  ggtitle("Nonlinearity in species-environment effects")

ggsave(nonlin_fig_C, filename = "paper_figs/nonlin_fig.jpg",
       width = 8, height = 3, dpi = "retina")

gam_summary_df %>% 
  filter(paste0(species, covariate) %in% paste0(incl_results$species, incl_results$covariate)) %>%
  .$rsq %>% 
  quantile(probs = c(0.1, 0.5, 0.9))


#### Finally: produce the predicted intensity rasters ####
library(coda)
library(nimble)
source("setup_code/define_fundamentals_nolineage.R")
source("main_code_nolineage/vis_fn.R")

target_species <- taxon_key$common_name_clean
target_species <- target_species[!target_species %in% species_to_exclude]

for (i in 1:length(target_species)) {
  writeLines(paste0(i, "..."))
  plot_spatpred_iSDM(target_species[i], thin = 20)
}

