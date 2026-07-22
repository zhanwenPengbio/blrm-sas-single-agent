# Independent R/JAGS check for the single-agent BLRM SAS implementation
#
# This script reproduces the matched comparison scenario used in the
# PharmaSUG China 2026 presentation. It intentionally implements the model
# directly in JAGS rather than calling a second BLRM-specific package.

required_packages <- c("rjags", "coda")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install JAGS 4.x and the following R packages before running: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

library(rjags)
library(coda)

settings <- list(
  observed_dose = c(2, 4),
  n = c(6L, 3L),
  dlt = c(1L, 0L),
  candidate_dose = c(2, 4, 6, 8, 12),
  ref_dose = 12,
  prior_mean = c(-0.9, 0),
  prior_cov = matrix(c(4, 0, 0, 1), nrow = 2, byrow = TRUE),
  underdose_bound = 0.15,
  overdose_bound = 0.35,
  ewoc_bound = 0.25,
  n_chains = 4L,
  n_adapt = 5000L,
  n_burn = 20000L,
  n_iter = 40000L,
  thin = 1L,
  chain_seeds = c(9527L, 9528L, 9529L, 9530L),
  rhat_cutoff = 1.01
)

validate_settings <- function(x) {
  if (length(x$observed_dose) != length(x$n) ||
      length(x$n) != length(x$dlt)) {
    stop("observed_dose, n, and dlt must have equal lengths.", call. = FALSE)
  }
  if (any(!is.finite(x$observed_dose)) || any(x$observed_dose <= 0)) {
    stop("All observed doses must be finite and positive.", call. = FALSE)
  }
  if (any(!is.finite(x$candidate_dose)) || any(x$candidate_dose <= 0)) {
    stop("All candidate doses must be finite and positive.", call. = FALSE)
  }
  if (!is.finite(x$ref_dose) || x$ref_dose <= 0) {
    stop("ref_dose must be finite and positive.", call. = FALSE)
  }
  if (any(x$n < 0) || any(x$dlt < 0) || any(x$dlt > x$n)) {
    stop("DLT counts must satisfy 0 <= dlt <= n.", call. = FALSE)
  }
  if (!identical(dim(x$prior_cov), c(2L, 2L)) ||
      any(!is.finite(x$prior_cov)) ||
      !isTRUE(all.equal(x$prior_cov, t(x$prior_cov)))) {
    stop("prior_cov must be a finite symmetric 2 x 2 matrix.", call. = FALSE)
  }
  if (min(eigen(x$prior_cov, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop("prior_cov must be positive definite.", call. = FALSE)
  }
  if (!(0 < x$underdose_bound &&
        x$underdose_bound < x$overdose_bound &&
        x$overdose_bound < 1)) {
    stop("Toxicity bounds must satisfy 0 < underdose < overdose < 1.", call. = FALSE)
  }
  if (!(0 < x$ewoc_bound && x$ewoc_bound < 1)) {
    stop("ewoc_bound must lie strictly between 0 and 1.", call. = FALSE)
  }
  if (length(x$chain_seeds) != x$n_chains) {
    stop("Provide one random-number seed per chain.", call. = FALSE)
  }
  invisible(x)
}

validate_settings(settings)

model_text <- "
model {
  theta[1:2] ~ dmnorm(prior_mean[], prior_precision[,])

  log_alpha <- theta[1]
  log_beta <- theta[2]
  beta <- exp(log_beta)

  for (i in 1:n_obs) {
    dlt[i] ~ dbin(p_obs[i], n[i])
    logit(p_obs[i]) <- log_alpha + beta * log(observed_log_ratio[i])
  }
}
"

jags_data <- list(
  n_obs = length(settings$observed_dose),
  dlt = settings$dlt,
  n = settings$n,
  observed_log_ratio = log(settings$observed_dose / settings$ref_dose),
  prior_mean = settings$prior_mean,
  prior_precision = solve(settings$prior_cov)
)

make_initial_values <- function(chain_id) {
  set.seed(settings$chain_seeds[chain_id])
  initial_theta <- as.numeric(
    settings$prior_mean +
      t(chol(settings$prior_cov)) %*% rnorm(2)
  )
  list(
    theta = initial_theta,
    .RNG.name = "base::Wichmann-Hill",
    .RNG.seed = settings$chain_seeds[chain_id]
  )
}

initial_values <- lapply(seq_len(settings$n_chains), make_initial_values)

model_connection <- textConnection(model_text)

jags_model <- jags.model(
  file = model_connection,
  data = jags_data,
  inits = initial_values,
  n.chains = settings$n_chains,
  n.adapt = settings$n_adapt,
  quiet = TRUE
)
close(model_connection)

update(jags_model, n.iter = settings$n_burn, progress.bar = "none")

posterior <- coda.samples(
  model = jags_model,
  variable.names = c("log_alpha", "log_beta"),
  n.iter = settings$n_iter,
  thin = settings$thin,
  progress.bar = "none"
)

posterior_matrix <- as.matrix(posterior)
log_alpha <- posterior_matrix[, "log_alpha"]
log_beta <- posterior_matrix[, "log_beta"]

toxicity_draws <- vapply(
  settings$candidate_dose,
  function(dose) {
    plogis(log_alpha + exp(log_beta) * log(dose / settings$ref_dose))
  },
  numeric(length(log_alpha))
)

ewoc_summary <- data.frame(
  dose = settings$candidate_dose,
  under = colMeans(toxicity_draws < settings$underdose_bound),
  target = colMeans(
    toxicity_draws >= settings$underdose_bound &
      toxicity_draws < settings$overdose_bound
  ),
  over = colMeans(toxicity_draws >= settings$overdose_bound),
  stringsAsFactors = FALSE
)

ewoc_summary$ewoc_acceptable <- ewoc_summary$over <= settings$ewoc_bound

eligible <- ewoc_summary[ewoc_summary$ewoc_acceptable, , drop = FALSE]
if (nrow(eligible) == 0L) {
  selected_dose <- NA_real_
} else {
  selected_dose <- eligible$dose[
    order(eligible$target, eligible$dose, decreasing = TRUE)[1L]
  ]
}

# Values reported for the matched SAS run in the project paper. These values
# are comparison references, not hard pass/fail thresholds for stochastic MCMC.
sas_reference <- data.frame(
  dose = c(2, 4, 6, 8, 12),
  sas_under = c(0.79960, 0.59673, 0.45635, 0.37768, 0.29348),
  sas_target = c(0.18490, 0.33718, 0.40083, 0.39900, 0.35905),
  sas_over = c(0.01550, 0.06610, 0.14283, 0.22333, 0.34748)
)

comparison <- merge(sas_reference, ewoc_summary, by = "dose", sort = FALSE)
comparison$abs_diff_under <- abs(comparison$sas_under - comparison$under)
comparison$abs_diff_target <- abs(comparison$sas_target - comparison$target)
comparison$abs_diff_over <- abs(comparison$sas_over - comparison$over)

gelman <- gelman.diag(
  posterior,
  confidence = 0.95,
  transform = FALSE,
  autoburnin = FALSE,
  multivariate = FALSE
)$psrf

diagnostics <- data.frame(
  parameter = rownames(gelman),
  rhat = gelman[, "Point est."],
  rhat_upper_95 = gelman[, "Upper C.I."],
  effective_sample_size = as.numeric(effectiveSize(posterior)),
  rhat_cutoff = settings$rhat_cutoff,
  converged = gelman[, "Point est."] <= settings$rhat_cutoff,
  row.names = NULL
)

output_directory <- file.path(getwd(), "output")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

write.csv(
  ewoc_summary,
  file.path(output_directory, "rjags_ewoc_summary.csv"),
  row.names = FALSE
)
write.csv(
  comparison,
  file.path(output_directory, "sas_rjags_comparison.csv"),
  row.names = FALSE
)
write.csv(
  diagnostics,
  file.path(output_directory, "rjags_diagnostics.csv"),
  row.names = FALSE
)

cat("\nR/JAGS posterior region probabilities\n")
print(ewoc_summary, row.names = FALSE, digits = 5)

cat("\nSelected dose under the static EWOC rule:", selected_dose, "mg\n")
cat(
  "Largest absolute SAS-versus-R/JAGS difference:",
  format(
    max(
      comparison$abs_diff_under,
      comparison$abs_diff_target,
      comparison$abs_diff_over
    ),
    digits = 5
  ),
  "\n"
)

cat("\nMCMC diagnostics\n")
print(diagnostics, row.names = FALSE, digits = 5)

cat("\nOutput files written to:", normalizePath(output_directory), "\n")
