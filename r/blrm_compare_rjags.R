#===============================================================
# BLRM R/JAGS re-check:
# 8 independent chains, each with 20,000 burn-in + 50,000 draws
#===============================================================

suppressPackageStartupMessages({
  library(rjags)
  library(coda)
})

set.seed(20260626)

#-----------------------------
# 1. Matched comparison setting
#-----------------------------
dose_obs <- c(2, 4)
n_obs    <- c(6, 3)
dlt_obs  <- c(1, 0)

dose_grid <- c(2, 4, 6, 8, 12)

ref_dose <- 12

mu <- c(-0.9, 0)

Sigma <- matrix(
  c(
    4, 0,
    0, 1
  ),
  nrow = 2,
  byrow = TRUE
)

Omega <- solve(Sigma)

ud   <- 0.15
od   <- 0.35
ewoc <- 0.25

n_chains <- 8
n_adapt  <- 2000
n_burnin <- 20000
n_draw   <- 50000

#-----------------------------
# 2. JAGS model
#-----------------------------
jags_model_text <- "
model {

  theta[1:2] ~ dmnorm(mu[], Omega[,])

  for (i in 1:Nobs) {
    dlt[i] ~ dbin(p[i], n[i])

    logit(p[i]) <-
      theta[1] +
      exp(theta[2]) * log(dose[i] / ref_dose)
  }

}
"

jags_data <- list(
  Nobs     = length(dose_obs),
  dose     = dose_obs,
  n        = n_obs,
  dlt      = dlt_obs,
  ref_dose = ref_dose,
  mu       = mu,
  Omega    = Omega
)

make_inits <- function(chain_id) {
  list(
    theta = c(
      rnorm(1, mean = mu[1], sd = 0.50),
      rnorm(1, mean = mu[2], sd = 0.25)
    ),
    .RNG.name = "base::Wichmann-Hill",
    .RNG.seed = 30000 + chain_id
  )
}

inits <- lapply(seq_len(n_chains), make_inits)

#-----------------------------
# 3. Run MCMC
#-----------------------------
jm <- jags.model(
  file     = textConnection(jags_model_text),
  data     = jags_data,
  inits    = inits,
  n.chains = n_chains,
  n.adapt  = n_adapt,
  quiet    = TRUE
)

update(
  object       = jm,
  n.iter       = n_burnin,
  progress.bar = "none"
)

samples <- coda.samples(
  model          = jm,
  variable.names = c("theta"),
  n.iter         = n_draw,
  thin           = 1,
  progress.bar   = "none"
)

#-----------------------------
# 4. Posterior region probabilities
#-----------------------------
calc_region_prob <- function(draw_mat) {
  
  log_alpha <- draw_mat[, "theta[1]"]
  log_beta  <- draw_mat[, "theta[2]"]
  
  p_mat <- sapply(
    dose_grid,
    function(d) {
      plogis(
        log_alpha +
          exp(log_beta) * log(d / ref_dose)
      )
    }
  )
  
  if (is.null(dim(p_mat))) {
    p_mat <- matrix(p_mat, ncol = 1)
  }
  
  out <- data.frame(
    dose    = dose_grid,
    UD_Mean = colMeans(p_mat < ud),
    TT_Mean = colMeans(p_mat >= ud & p_mat < od),
    OD_Mean = colMeans(p_mat >= od)
  )
  
  out$EWOC   <- out$OD_Mean <= ewoc
  out$check1 <- rowSums(out[, c("UD_Mean", "TT_Mean", "OD_Mean")])
  
  out
}

# Pooled result across all chains
r_region_pooled <- calc_region_prob(as.matrix(samples))

# Chain-specific region probabilities
r_region_by_chain <- do.call(
  rbind,
  lapply(seq_along(samples), function(chain_id) {
    
    tmp <- calc_region_prob(as.matrix(samples[[chain_id]]))
    tmp$chain <- chain_id
    
    tmp[, c(
      "chain", "dose",
      "UD_Mean", "TT_Mean", "OD_Mean",
      "EWOC", "check1"
    )]
  })
)

#-----------------------------
# 5. Chain stability summary
#-----------------------------
r_chain_stability <- do.call(
  rbind,
  lapply(split(r_region_by_chain, r_region_by_chain$dose), function(x) {
    
    data.frame(
      dose = x$dose[1],
      
      UD_Mean_AcrossChain = mean(x$UD_Mean),
      UD_SD_AcrossChain   = sd(x$UD_Mean),
      UD_Min_AcrossChain  = min(x$UD_Mean),
      UD_Max_AcrossChain  = max(x$UD_Mean),
      
      TT_Mean_AcrossChain = mean(x$TT_Mean),
      TT_SD_AcrossChain   = sd(x$TT_Mean),
      TT_Min_AcrossChain  = min(x$TT_Mean),
      TT_Max_AcrossChain  = max(x$TT_Mean),
      
      OD_Mean_AcrossChain = mean(x$OD_Mean),
      OD_SD_AcrossChain   = sd(x$OD_Mean),
      OD_Min_AcrossChain  = min(x$OD_Mean),
      OD_Max_AcrossChain  = max(x$OD_Mean)
    )
  })
)

row.names(r_chain_stability) <- NULL

#-----------------------------
# 6. Convergence diagnostics
#-----------------------------
psrf <- gelman.diag(
  samples,
  autoburnin   = FALSE,
  multivariate = FALSE
)$psrf

ess <- effectiveSize(samples)

r_mcmc_diagnostic <- data.frame(
  parameter  = rownames(psrf),
  Rhat       = psrf[, "Point est."],
  Rhat_upper = psrf[, "Upper C.I."],
  ESS        = as.numeric(ess[rownames(psrf)])
)

#-----------------------------
# 7. Display and save outputs
#-----------------------------
options(digits = 7)

cat("\n================ R pooled region probabilities ================\n")
print(r_region_pooled)

cat("\n================ R chain stability ================\n")
print(r_chain_stability)

cat("\n================ R MCMC diagnostics ================\n")
print(r_mcmc_diagnostic)

cat("\n================ R chain-specific region probabilities ================\n")
print(r_region_by_chain)

write.csv(
  r_region_pooled,
  "r_blrm_region_pooled.csv",
  row.names = FALSE
)

write.csv(
  r_chain_stability,
  "r_blrm_chain_stability.csv",
  row.names = FALSE
)

write.csv(
  r_mcmc_diagnostic,
  "r_blrm_mcmc_diagnostic.csv",
  row.names = FALSE
)

write.csv(
  r_region_by_chain,
  "r_blrm_region_by_chain.csv",
  row.names = FALSE
)

