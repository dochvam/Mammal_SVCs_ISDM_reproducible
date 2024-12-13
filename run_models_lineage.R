###############################################################################
#'
#' run_models_nolineage.R
#' Author: Ben R. Goldstein
#' Description: Running this script will trigger the full model estimation
#'    workflow for all 13 lineage species, including information about lineage.
#'
###############################################################################

library(nimbleEcology)
library(tidyverse)
library(parallel)
library(coda)

# Load helper fns
source("main_code_lineage/model_fn.R")
source("main_code_lineage/model_code_main.R")
source("setup_code/define_fundamentals.R")

# Retrieve the list of species
target_species <- taxon_key$common_name_clean

result_files <- list.files("intermediate/integration_results/",
                           pattern = "_lineage",
                           full.names = TRUE)
if (length(result_files) == 0) {
  this_round <- 1
} else {
  this_round <- max(parse_number(result_files)) + 1
}

res_list <- list()

for (i in 1:length(target_species)) {
  this_results <- result_files[grep(target_species[i], result_files)]
  
  if (length(this_results) > 0) {
    temp <- lapply(this_results, function(x) {
      thisres <- readRDS(x)
      thisres$samples_list[[1]]$samples %>% as.list()
    }) %>% 
      do.call(what = c)
    
    temp <- mcmc.list(temp)
    
    res_list[[i]] <- temp %>% 
      MCMCvis::MCMCsummary() %>% 
      .["logDens",] %>% 
      mutate(nres = length(this_results),
             species = target_species[i])
  }
}

res_df <- bind_rows(res_list) %>% 
  mutate(finished = Rhat <= 1.1 | n.eff >= 250)

target_species <- target_species[!target_species %in% res_df$species[res_df$finished]]

outfiles <- c()
ct <- 0
overwrite <- F

set.seed(2976746)
seed_vec <- floor(runif(100, 0, 1) * 100000)


### Now, we create a single script in the temp/source directory for each species
### so that all species models can be estimated in parallel

# Beginning info:
preamble <- '
library(nimbleEcology)
library(tidyverse)

source("main_code_lineage/model_fn.R")
source("main_code_lineage/model_code_main.R")
source("setup_code/define_fundamentals.R")

target_species <- taxon_key$common_name_clean

i <- '

# Loop over species
outfiles <- c()
ct <- 0

for (i in 1:length(target_species)) {
  ct <- ct + 1
  this_outfile <- paste0("temp/source/script", ct, ".R")
  outfiles <- c(outfiles, this_outfile)
  
  cat(paste0(preamble, i, '
  
  fit_integrated_model_with_CV(
    nfolds = 0               ,
    inat_disttype = "NB",
    modtype = "joint",
    species = "', target_species[i], '",
    spatial_model = "SVCs",
    ni = 10000, nb = 2000,
    nc = 3, nt = 2, nt2 = 10,
    seed = ', seed_vec[this_round], ', subset_CT = 1,
    subset_inat = 1, suffix = "_main_', this_round, '"
  )
  '), file = this_outfile)
  
}

cl <- makeCluster(14)

# Execute the model estimates
parLapply(cl, outfiles, source)

stopCluster(cl)



