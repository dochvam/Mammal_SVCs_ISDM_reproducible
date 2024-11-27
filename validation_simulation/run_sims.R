library(tidyverse)

set.seed(7228114)
source("validation_simulation/simulate_fn.R")

example_specs <- c("eastern_chipmunk", "bobcat")

# Three species, three scenarios, three replicates per
tasks <- data.frame(expand.grid(species = example_specs,
                                scenario = 1:3,
                                iter = 1:8)) %>%
  mutate(seed = as.integer(runif(48) * 10000))

files_tosource <- c()
ct <- 0
overwrite <- F

for (i in 1:nrow(tasks)) {
  outfile <- file.path("intermediate", "validation_sim",
                       paste0("VS_", 
                              tasks$species[i], "_S", 
                              tasks$scenario[i], "_I", 
                              tasks$iter[i], ".RDS"))
  
  if (overwrite || !file.exists(outfile)) {
    ct <- ct + 1
    
    files_tosource[ct] <- paste0("temp/source/VS_source_", ct, ".R")
    
    writeLines(paste0(
      '
library(tidyverse)
source("validation_simulation/simulate_fn.R")

run_simulation(
  scenario = ', tasks$scenario[i], ',
  species = "', tasks$species[i], '",
  seed = ', tasks$seed[i], ',
  iter = ',tasks$iter[i], ',
  ni = 10000, nb = 2000,
  nc = 5, nt = 2, nt2 = 10
)
  
      '),
      con = files_tosource[ct]
    )
  }
}

library(parallel)
cl <- makeCluster(10)
parLapply(cl, sample(files_tosource), source)
stopCluster(cl)

