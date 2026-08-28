# Minimal cmdstanr example: prior sensitivity, then interpret the result by
# refitting the Stan model at the worst-case candidate `lambda_max` and
# comparing it against the reference fit with a quantile table, a KDE plot,
# and an ECDF plot.

library(cmdstanr)
library(fdsens)
source("interpretation/plots.R")

set.seed(123)
y <- rnorm(30, mean = 1)
prior_mean_ref <- 0
prior_sd_ref <- 2

stan_file <- system.file("stan", "gaussian_location.stan", package = "fdsens")
model <- cmdstan_model(stan_file, force_recompile = TRUE)

fit_reference <- model$sample(
  data = list(
    N = length(y),
    y = y,
    prior_mean = prior_mean_ref,
    prior_sd = prior_sd_ref
  ),
  seed = 123,
  chains = 2,
  parallel_chains = 2,
  iter_warmup = 500,
  iter_sampling = 1000,
  refresh = 0
)

# The Stan statement `theta ~ normal(prior_mean, prior_sd)` is detected as an
# exponential-family prior. Bounds are therefore supplied in natural
# parameters lambda1 = mu / sigma^2 and lambda2 = -1 / (2 sigma^2).
prior_result <- fd_prior_global_sensitivity(
  fit = fit_reference,
  variables = "theta",
  lambda_lower = c(lambda1 = -1, lambda2 = -2),
  lambda_upper = c(lambda1 = 1, lambda2 = -0.05),
  stan_file = stan_file,
  prior_variable = "theta",
  stan_data = list(prior_mean = prior_mean_ref, prior_sd = prior_sd_ref)
)
print(prior_result)

# `lambda_max` is in the normal prior's natural parameterisation
# (lambda1 = mu / sigma^2, lambda2 = -1 / (2 sigma^2)); invert it to refit
# the candidate model in its original (mu, sigma) parameterisation.
lambda_max <- prior_result$lambda_max
sigma_candidate <- sqrt(-1 / (2 * lambda_max["lambda2"]))
mu_candidate <- lambda_max["lambda1"] * sigma_candidate^2

fit_candidate <- model$sample(
  data = list(
    N = length(y),
    y = y,
    prior_mean = as.numeric(mu_candidate),
    prior_sd = as.numeric(sigma_candidate)
  ),
  seed = 123,
  chains = 2,
  parallel_chains = 2,
  iter_warmup = 500,
  iter_sampling = 1000,
  refresh = 0
)

fits <- list(reference = fit_reference, candidate = fit_candidate)
output_dir <- "interpretation/output/gaussian_location_prior"

save_sensitivity_result(prior_result, output_dir = output_dir)
quantile_table <- plot_quantiles(fits, variables = "theta", output_dir = output_dir)
print(quantile_table)
plot_kde(fits, variables = "theta", output_dir = output_dir)
plot_ecdf(fits, variables = "theta", output_dir = output_dir)

cat("Wrote quantile table and plots to", output_dir, "\n")
