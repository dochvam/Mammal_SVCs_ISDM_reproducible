###############################################################################
#'
#' main_code_nolineage/model_code_main.R
#' Author: Ben R. Goldstein
#' Description: Helper code defining NIMBLE models and NIMBLE functions
#'    needed to run the main analysis. The core model defn is called
#'    joint_model_SVCs_asRanefs
#'
###############################################################################

library(nimbleEcology)

sumLogProb <- nimbleFunction(
  name = 'sumLogProb',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    stochNodes <- model$getNodeNames(stochOnly = TRUE)
    stochNodes <- stochNodes[stochNodes != target]
    stochNodes <- stochNodes[!grepl("_0", stochNodes)]
  },
  run = function() {
    model[[target]] <<- model$getLogProb(stochNodes)
  },
  methods = list( reset = function() {} )
)

#### calcBrierScore_oneSite ####
calcBrierScore_oneSite <- nimbleFunction(
  run = function(y = double(1),
                 psi = double(0),
                 p = double(1),
                 len = integer(0)){
    
    x <- numeric(len)
    for (i in 1:len) {
      x[i] <- y[i] * (1 - psi * p[i])^2 + (1 - y[i]) * (psi * p[i])^2
    }

    return(sum(x))
    returnType(double(0))
  })


#### calcIntensity ####

calcIntensity_weak <- nimbleFunction(run = function(
    lambda_intercept = double(0),
    theta0 = double(0),
    theta1 = double(0),
    lambda_beta = double(1),
    xdat = double(2),
    spatial_ranef = double(0)) {
  
  log_lambda <- log(lambda_intercept) + xdat %*% lambda_beta + spatial_ranef
  
  mu <- exp(theta0 + theta1 * log(sum(exp(log_lambda))))
  
  return(mu)
  returnType(double(0))
})

#### Joint model with SVCs and eco as random effects ####

joint_model_SVCs_asRanefs <- nimbleCode({
  
  ### iNaturalist data distribution
  for (g in 1:nINatCell) {
    if (inat_disttype == "NB") {
      inat_ct[g] ~ dnbinom(size = 1 / overdisp_inat,
                           prob = 1 / (1 + overdisp_inat * inat_effort[g] * mu[g]))
      
    } else {
      inat_ct[g] ~ dpois(inat_effort[g] * mu[g])
    }
    
    ### iNaturalist grid cell intensities
    mu[g] <- calcIntensity_weak(
      lambda_intercept = lambda_intercept,
      theta0 = theta0, theta1 = theta1,
      lambda_beta = lambda_beta[spatcell_inat[g], 1:nLamBeta],
      xdat = xdat_allS2[S2toS3_start[g]:S2toS3_end[g], 1:nLamBeta],
      spatial_ranef = spat_ranef[spatcell_inat[g]] * spat_magnitude_intercept
    )
  }
  
  ### Detection probabilities
  # Separate detection beta models
  for (o in 1:nobs) {
    cloglog(p[o]) <- 
      cloglog(p_intercept) + 
      inprod(p_alpha[1:nLamBeta], xdat[deployment_ID[o], 1:nLamBeta]) +
      inprod(p_gamma[1:nGamma], wdat[o, 1:nGamma]) +
      det_ranef[deployment_ID[o]]
  }
  
  for (o in 1:nobs_holdout) {
    cloglog(p_holdout[o]) <- 
      cloglog(p_intercept) + 
      inprod(p_alpha[1:nLamBeta], xdat[deployment_ID_holdout[o], 1:nLamBeta]) +
      inprod(p_gamma[1:nGamma], wdat_holdout[o, 1:nGamma]) +
      det_ranef[deployment_ID_holdout[o]]
  }
  
  
  # Occupancy probabilities
  for (i in 1:ncameras_indat) {
    cloglog(psi[i]) <- log(lambda_intercept) +
      inprod(lambda_beta[spatcell[deployment_ID[start[i]]], 1:nLamBeta], 
             xdat[deployment_ID[start[i]], 1:nLamBeta]) +
      spat_ranef[spatcell[deployment_ID[start[i]]]] * spat_magnitude_intercept
  }
  
  for (i in 1:ncameras_holdout) {
    cloglog(psi_holdout[i]) <- log(lambda_intercept) +
      inprod(lambda_beta[spatcell[deployment_ID_holdout[start_holdout[i]]], 1:nLamBeta], 
             xdat[deployment_ID_holdout[start_holdout[i]], 1:nLamBeta]) +
      spat_ranef[spatcell[deployment_ID_holdout[start_holdout[i]]]] * spat_magnitude_intercept
  }
  
  
  # Occupancy model distributions
  for (i in 1:ncameras_indat) {
    y[start[i]:end[i]] ~ dOcc_v(
      probOcc = psi[i],
      probDetect = p[start[i]:end[i]],
      len = end[i] - start[i] + 1
    )
  }
  for (i in 1:ncameras_holdout) {
    ## Can add this in to make the holdout data contribute to the likelihood (i.e. no holdout data):
    if (!holdout_data) {
      y_holdout[start_holdout[i]:end_holdout[i]] ~ dOcc_v(
        probOcc = psi_holdout[i],
        probDetect = p_holdout[start_holdout[i]:end_holdout[i]],
        len = end_holdout[i] - start_holdout[i] + 1
      )
    }
    
    pp_brier_bysite[i] <- calcBrierScore_oneSite(
      y = y_holdout[start_holdout[i]:end_holdout[i]],
      psi = psi_holdout[i],
      p = p_holdout[start_holdout[i]:end_holdout[i]],
      len = end_holdout[i] - start_holdout[i] + 1
    )
  }
  
  ### Brier score  
  brier_score <- sum(pp_brier_bysite[1:ncameras_holdout]) / nobs_holdout
  
  ### Priors we always need
  p_intercept      ~ dunif(0, 1)
  lambda_intercept ~ dgamma(shape = 1, scale = 50)
  
  # SVC Bernoulli probabilities
  p_incl_eco ~ dunif(0, 1)
  p_incl_CAR ~ dunif(0, 1)
  
  for (i in 1:nLamBeta) {
    # Prior on overall mean effect
    lambda_beta_mean[i] ~ dnorm(0, sd = 2.75)
    
    # indicator variables for SVCs
    eco_gamma[i] ~ dbern(p_incl_eco)
    CAR_gamma[i] ~ dbern(p_incl_CAR)
    
    # CAR structure in the effect of each covar
    CAR_effect[1:nSpatCell, i] ~ dcar_normal(adj = adj[1:adj_len],
                                             num = num[1:nSpatCell],
                                             tau = 1, 
                                             zero_mean = 1)
    for (j in 1:nEcoregion) {
      eco_effect[j, i] ~ dnorm(0, sd = 1)
    }
    eco_effect[nEcoregion + 1, i] <- 0
    
    # Prior on psi mean/constant effect, on spatial variance when incl.
    sigma_CAR_0[i] ~ dunif(0, 10)
    sigma_eco_0[i] ~ dunif(0, 10)
    
    # Actual spatial sigma
    sigma_CAR[i] <- sigma_CAR_0[i] * CAR_gamma[i] + cc
    sigma_eco[i] <- sigma_eco_0[i] * eco_gamma[i] + cc
    
    # Grid-cell-level betas
    for (cell in 1:nSpatCell) {
      lambda_beta[cell, i] <- 
        lambda_beta_mean[i] +
        sigma_CAR[i] * CAR_effect[cell, i] +
        sigma_eco[i] * eco_effect[ecoregion[cell], i]
    }
  }
  
  # Detection variables
  for (i in 1:nGamma) {
    p_gamma[i] ~ dnorm(0, sd = 2.75)
  }
  for (i in 1:nLamBeta) {
    p_alpha[i] ~ dnorm(0, sd = 2.75)
  }
  # Detection random effect
  for (i in 1:ndepl_total) {
    det_ranef[i] ~ dnorm(0, sd = det_ranef_sd)
  }
  det_ranef_sd ~ dinvgamma(shape = 0.5, scale = 0.5)

  # Spatial priors for intercept
  spat_magnitude_intercept ~ dinvgamma(shape = 0.5, scale = 0.5)
  spat_ranef[1:nSpatCell] ~ dcar_normal(adj = adj[1:adj_len],
                                        num = num[1:nSpatCell],
                                        tau = 1, zero_mean = 1)
  
  # Secondary data correlation params
  theta0 ~ dnorm(0, sd = 10)
  theta1 ~ dnorm(1, sd = 2.75)
  # theta0_museum ~ dnorm(0, sd = 10)
  # theta1_museum ~ dnorm(1, sd = 2.75)
  
  # Overdispersion parameter if we're using negative binomial
  if (inat_disttype == "NB") {
    overdisp_inat ~ dgamma(shape = 1, scale = 5)
    # overdisp_museum ~ dgamma(shape = 1, scale = 5)
  }
  
  # Create a node which will track model log density
  logDens ~ dflat()
})
