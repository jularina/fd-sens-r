data {
  int<lower=1> K;
  int<lower=1> T;
  array[T] real y;
}

parameters {
  real alpha;
  array[K] real beta;
  real<lower=0> sigma;
}

model {
  // Independent reference priors from the Kilpisjarvi AR(5) analysis.
  alpha ~ normal(0, 5);
  beta[1] ~ normal(0, 5);
  beta[2] ~ normal(0, 5);
  beta[3] ~ normal(0, 5);
  beta[4] ~ normal(0, 5);
  beta[5] ~ normal(0, 5);
  sigma ~ cauchy(0, 1);

  for (t in (K + 1):T) {
    real mu = alpha;
    for (k in 1:K) mu += beta[k] * y[t - k];
    y[t] ~ normal(mu, sigma);
  }
}

