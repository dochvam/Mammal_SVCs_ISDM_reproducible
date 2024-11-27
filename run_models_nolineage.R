###############################################################################
#'
#' run_models_nolineage.R
#' Author: Ben R. Goldstein
#' Description: Running this script will trigger the full model estimation
#'    workflow for all 35 study species, not accounting for information about
#'    lineage.
#'
###############################################################################


# Load necessary packages
library(nimbleEcology)
library(tidyverse)
library(parallel)

# Load helper fns
source("main_code_nolineage/model_fn.R")
source("main_code_nolineage/model_code_main.R")
source("setup_code/define_fundamentals_nolineage.R")

# Don't re-run models if outputs exist
overwrite <- FALSE

# Retrieve the list of species
target_species <- taxon_key$common_name_clean

### Now, we create a single script in the temp/source directory for each species
### so that all species models can be estimated in parallel

# Beginning info:
preamble <- '
library(nimbleEcology)
library(tidyverse)

source("main_code_nolineage/model_fn.R")
source("main_code_nolineage/model_code_main.R")
source("setup_code/define_fundamentals_nolineage.R")

target_species <- taxon_key$common_name_clean

'

# Loop over species
outfiles <- c()
ct <- 0
for (i in 1:length(target_species)) {
  
  # Output file name:
  result_outfile <- paste0("intermediate/integration_results/final_", 
                           target_species[i], "_joint_SVCs_main_noLineage_samples.RDS")
  
  if (!file.exists(result_outfile) || overwrite) {
    
    ct <- ct + 1
    this_outfile <- paste0("temp/source/script", ct, ".R")
    outfiles <- c(outfiles, this_outfile)
    
    # Write out a file that can be run in a clean session to estimate one spec.
    cat(paste0(preamble, '
    
    write(paste0("Starting species ", "', target_species[i], '",
                 " at time ", Sys.time()), file = "log.txt", append = TRUE)
  
    fit_integrated_model_with_CV(
      nfolds = 0               ,
      inat_disttype = "NB",
      modtype = "joint",
      species = "', target_species[i], '",
      spatial_model = "SVCs",
      ni = 10000, nb = 2000,
      nc = 5, nt = 2, nt2 = 10,
      seed = 3063446, subset_CT = 1,
      subset_inat = 1, suffix = "_main",
      overwrite = ', overwrite, '
    )
    '), file = this_outfile)
  }
}

write("============== New run ==============", file = "log.txt", append = TRUE)

cl <- makeCluster(4)

# Execute the model estimates
parLapply(cl, outfiles, source)

stopCluster(cl)
