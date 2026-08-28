# Kilpisjarvi AR(5): sensitivity decomposition under independent priors.
# This example uses rstan to mirror the original analysis.

library(rstan)
library(fdsens)
source("interpretation/plots.R")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

K_val <- 5
y_kilpisjarvi <- c(
  8.3, 10.9, 9.4, 8.1, 8.1, 7.7, 8.6, 9.1, 11.0, 10.1,
  7.6, 8.8, 8.3, 7.2, 9.3, 8.8, 7.6, 10.5, 11.0, 8.9,
  11.3, 10.0, 10.1, 6.4, 8.2, 8.4, 9.5, 9.9, 10.6, 7.6,
  7.7, 8.1, 8.4, 9.7, 9.5, 7.3, 10.3, 9.6, 10.3, 9.8,
  9.0, 9.1, 9.5, 8.7, 9.9, 10.5, 9.4, 9.0, 9.0, 9.7,
  11.4, 10.7, 10.1, 10.8, 10.4, 10.3, 8.8, 9.8, 8.8, 10.8,
  8.6, 11.1
)
y_centered <- y_kilpisjarvi - mean(y_kilpisjarvi)
stan_data <- list(K = K_val, T = length(y_centered), y = y_centered)

stan_file <- system.file("stan", "kilpisjarvi_ar5.stan", package = "fdsens")
fit <- rstan::stan(
  file = stan_file,
  data = stan_data,
  chains = 5,
  warmup = 1000,
  iter = 5000,
  seed = 123,
  refresh = 100
)

# Gaussian candidate priors use natural parameters
# eta1 = mu / sd^2 and eta2 = -1 / (2 sd^2).
normal_block <- function(variable) {
  list(
    variables = variable,
    prior_variable = variable,
    candidate_family = "normal",
    lambda_lower = c(eta1 = -8, eta2 = -8),
    lambda_upper = c(eta1 = 8, eta2 = -0.02)
  )
}

blocks <- list(
  alpha = normal_block("alpha"),
  beta1 = normal_block("beta[1]"),
  beta2 = normal_block("beta[2]"),
  beta3 = normal_block("beta[3]"),
  beta4 = normal_block("beta[4]"),
  beta5 = normal_block("beta[5]")
)

# The reference prior for sigma is half-Cauchy(0, 1), represented in Stan by
# a positive parameter with `sigma ~ cauchy(0, 1)`. Candidate priors are
# inverse-Gamma, with eta1 = -(shape + 1) and eta2 = -scale.
blocks$sigma <- list(
  variables = "sigma",
  prior_variable = "sigma",
  candidate_family = "inv_gamma",
  lambda_lower = c(eta1 = -8, eta2 = -2),
  lambda_upper = c(eta1 = -3.5, eta2 = -1 / 6),
  score_prior_ref = function(draws) {
    sigma <- draws[, "sigma"]
    matrix(-2 * sigma / (1 + sigma^2), ncol = 1)
  }
)

result <- fd_prior_global_sensitivity(
  fit = fit,
  independent = TRUE,
  blocks = blocks,
  stan_file = stan_file,
  stan_data = stan_data
)

print(result)

# Seven two-parameter blocks require 7 * 2^2 = 28 corner evaluations instead
# of 2^14 = 16,384 corners for a joint 14-dimensional box.
stopifnot(sum(result$components$corner_evaluations) == 28)

output_dir <- "interpretation/output/kilpisjarvi_ar5"
save_sensitivity_result(result, output_dir = output_dir)
plot_component_shares(result, output_dir = output_dir)

# Refit at the worst-case candidate: invert each block's `lambda_max` back to
# its original Stan parameterisation. Normal blocks (alpha, beta[1:K]) use
# eta1 = mu / sd^2, eta2 = -1 / (2 sd^2); the sigma block's candidate family
# is inverse-Gamma, with eta1 = -(shape + 1), eta2 = -scale.
normal_to_original <- function(lambda_max) {
  sd_component <- sqrt(-1 / (2 * lambda_max["eta2"]))
  mean_component <- lambda_max["eta1"] * sd_component^2
  c(mean = as.numeric(mean_component), sd = as.numeric(sd_component))
}

alpha_candidate <- normal_to_original(result$lambda_max$alpha)
beta_candidate <- vapply(
  paste0("beta", seq_len(K_val)), function(name) normal_to_original(result$lambda_max[[name]]),
  numeric(2)
)
sigma_lambda_max <- result$lambda_max$sigma
sigma_shape_candidate <- as.numeric(-sigma_lambda_max["eta1"] - 1)
sigma_scale_candidate <- as.numeric(-sigma_lambda_max["eta2"])

candidate_stan_file <- system.file(
  "stan", "kilpisjarvi_ar5_candidate.stan", package = "fdsens"
)
fit_candidate <- rstan::stan(
  file = candidate_stan_file,
  data = list(
    K = K_val,
    T = length(y_centered),
    y = y_centered,
    prior_mean_alpha = alpha_candidate["mean"],
    prior_sd_alpha = alpha_candidate["sd"],
    prior_mean_beta = beta_candidate["mean", ],
    prior_sd_beta = beta_candidate["sd", ],
    sigma_shape = sigma_shape_candidate,
    sigma_scale = sigma_scale_candidate
  ),
  chains = 5,
  warmup = 1000,
  iter = 5000,
  seed = 123,
  refresh = 100
)

fits <- list(reference = fit, candidate = fit_candidate)
variables <- c("alpha", paste0("beta[", seq_len(K_val), "]"), "sigma")

quantile_table <- plot_quantiles(fits, variables = variables, output_dir = output_dir)
print(quantile_table)
plot_kde(fits, variables = variables, output_dir = output_dir)
plot_ecdf(fits, variables = variables, output_dir = output_dir)

cat("Wrote sensitivity result and plots to", output_dir, "\n")
