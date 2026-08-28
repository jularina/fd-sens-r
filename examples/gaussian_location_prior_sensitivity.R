# Minimal cmdstanr example: prior sensitivity.

library(cmdstanr)
library(fdsens)

set.seed(123)
y <- rnorm(30, mean = 1)
prior_mean_ref <- 0
prior_sd_ref <- 2

stan_file <- system.file("stan", "gaussian_location.stan", package = "fdsens")
model <- cmdstan_model(stan_file, force_recompile = TRUE)

fit <- model$sample(
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
  fit = fit,
  variables = "theta",
  lambda_lower = c(lambda1 = -1, lambda2 = -2),
  lambda_upper = c(lambda1 = 1, lambda2 = -0.05),
  stan_file = stan_file,
  prior_variable = "theta",
  stan_data = list(prior_mean = prior_mean_ref, prior_sd = prior_sd_ref)
)
print(prior_result)
