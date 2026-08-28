data {
  int<lower=1> K;
  int<lower=1> T;
  array[T] real y;
  real prior_mean_alpha;
  real<lower=0> prior_sd_alpha;
  array[K] real prior_mean_beta;
  array[K] real<lower=0> prior_sd_beta;
  real<lower=0> sigma_shape;
  real<lower=0> sigma_scale;
}

parameters {
  real alpha;
  array[K] real beta;
  real<lower=0> sigma;
}

model {
  // Candidate priors: same normal family as the reference for alpha/beta,
  // with candidate hyperparameters; inverse-Gamma for sigma, whose reference
  // prior is half-Cauchy(0, 1) (a different family, see the independent
  // decomposition example).
  alpha ~ normal(prior_mean_alpha, prior_sd_alpha);
  for (k in 1:K) beta[k] ~ normal(prior_mean_beta[k], prior_sd_beta[k]);
  sigma ~ inv_gamma(sigma_shape, sigma_scale);

  for (t in (K + 1):T) {
    real mu = alpha;
    for (k in 1:K) mu += beta[k] * y[t - k];
    y[t] ~ normal(mu, sigma);
  }
}
