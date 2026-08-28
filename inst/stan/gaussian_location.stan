data {
  int<lower=1> N;
  vector[N] y;
  real prior_mean;
  real<lower=0> prior_sd;
}

parameters {
  real theta;
}

model {
  theta ~ normal(prior_mean, prior_sd);
  y ~ normal(theta, 1);
}

