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

# Load helper fns
source("main_code_lineage/model_fn.R")
source("main_code_lineage/model_code_main.R")
source("setup_code/define_fundamentals.R")

# Retrieve the list of species
target_species <- taxon_key$common_name_clean

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
    nc = 5, nt = 2, nt2 = 10,
    seed = 8659522, subset_CT = 1,
    subset_inat = 1, suffix = "_main"
  )
  '), file = this_outfile)
  
}

cl <- makeCluster(4)

# Execute the model estimates
parLapply(cl, outfiles, source)

stopCluster(cl)



