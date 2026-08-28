# Minimal cmdstanr example: learning-rate sensitivity.

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

# For learning-rate sensitivity the prior is fixed. The unscaled negative
# log-likelihood is 0.5 * sum((y - theta)^2), whose gradient is below.
score_loss <- function(draws) {
  theta <- draws[, "theta"]
  matrix(length(y) * theta - sum(y), ncol = 1)
}

lr_result <- fd_lr_global_sensitivity(
  fit = fit,
  variables = "theta",
  lambda_ref = 1,
  lower = 0.5,
  upper = 1.5,
  score_loss = score_loss
)
print(lr_result)
