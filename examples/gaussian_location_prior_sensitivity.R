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

# For prior sensitivity the likelihood is fixed, so only prior scores are
# needed; the common likelihood score cancels.
score_prior_ref <- function(draws) {
  theta <- draws[, "theta"]
  matrix((prior_mean_ref - theta) / prior_sd_ref^2, ncol = 1)
}

score_prior_candidate <- function(draws, lambda) {
  theta <- draws[, "theta"]
  candidate_mean <- lambda["prior_mean"]
  candidate_sd <- lambda["prior_sd"]
  matrix((candidate_mean - theta) / candidate_sd^2, ncol = 1)
}

prior_result <- fd_prior_global_sensitivity(
  fit = fit,
  variables = "theta",
  lower = c(prior_mean = -2, prior_sd = 0.5),
  upper = c(prior_mean = 2, prior_sd = 4),
  score_prior_ref = score_prior_ref,
  score_prior_candidate = score_prior_candidate
)
print(prior_result)
