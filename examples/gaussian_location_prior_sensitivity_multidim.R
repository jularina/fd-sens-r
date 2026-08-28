# Minimal cmdstanr example: prior sensitivity with a 2-D mean and a 2x2
# covariance. fd_prior_global_sensitivity() box-constrains each hyperparameter
# independently, so it has no way to keep a searched covariance positive
# definite on its own -- that is on the caller. Here the box on the 3 unique
# covariance entries (sigma11, sigma22, sigma12) is chosen so every candidate
# covariance in the box is guaranteed PD: with sigma11, sigma22 >= 3 and
# |sigma12| <= 1, det = sigma11 * sigma22 - sigma12^2 >= 3 * 3 - 1^2 = 8 > 0
# everywhere in the box.

library(cmdstanr)
library(fdsens)

set.seed(123)
N <- 30
true_theta <- c(1, -1)
y <- cbind(
  rnorm(N, true_theta[1], 1),
  rnorm(N, true_theta[2], 1)
)

prior_mean_ref <- c(0, 0)
prior_cov_ref <- matrix(c(4, 1, 1, 4), nrow = 2)

# Check the reference covariance is positive semi-definite before fitting,
# rather than relying on Stan's runtime check alone.
stopifnot(all(eigen(prior_cov_ref, symmetric = TRUE, only.values = TRUE)$values >= 0))

stan_file <- system.file("stan", "gaussian_location_multidim.stan", package = "fdsens")
model <- cmdstan_model(stan_file, force_recompile = TRUE)

fit <- model$sample(
  data = list(
    N = N,
    y = y,
    prior_mean = prior_mean_ref,
    prior_cov = prior_cov_ref
  ),
  seed = 123,
  chains = 2,
  parallel_chains = 2,
  iter_warmup = 500,
  iter_sampling = 1000,
  refresh = 0
)

prior_cov_ref_inv <- solve(prior_cov_ref)

# Gradient of log N(theta; prior_mean_ref, prior_cov_ref), one row per draw.
score_prior_ref <- function(draws) {
  theta <- draws[, c("theta[1]", "theta[2]")]
  diff <- matrix(prior_mean_ref, nrow = nrow(theta), ncol = 2, byrow = TRUE) - theta
  diff %*% prior_cov_ref_inv
}

# Candidate hyperparameters: lambda = (mu1, mu2, sigma11, sigma22, sigma12),
# i.e. the 2-D mean and the 3 unique entries of the 2x2 covariance.
score_prior_candidate <- function(draws, lambda) {
  theta <- draws[, c("theta[1]", "theta[2]")]
  mu <- c(lambda["mu1"], lambda["mu2"])
  sigma <- matrix(
    c(lambda["sigma11"], lambda["sigma12"], lambda["sigma12"], lambda["sigma22"]),
    nrow = 2
  )
  diff <- matrix(mu, nrow = nrow(theta), ncol = 2, byrow = TRUE) - theta
  diff %*% solve(sigma)
}

lower <- c(mu1 = -2, mu2 = -2, sigma11 = 3, sigma22 = 3, sigma12 = -1)
upper <- c(mu1 = 2, mu2 = 2, sigma11 = 5, sigma22 = 5, sigma12 = 1)

# Confirm the box's worst case (smallest variances, largest |covariance|) is
# still positive definite before handing it to the search.
worst_case_det <- lower["sigma11"] * lower["sigma22"] -
  max(abs(lower["sigma12"]), abs(upper["sigma12"]))^2
stopifnot(worst_case_det > 0)

prior_result <- fd_prior_global_sensitivity(
  fit = fit,
  variables = "theta",
  lambda_lower = lower,
  lambda_upper = upper,
  method = "black_box",
  score_prior_ref = score_prior_ref,
  score_prior_candidate = score_prior_candidate
)
print(prior_result)
