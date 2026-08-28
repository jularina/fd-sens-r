data {
  int<lower=1> N;
  array[N] vector[2] y;
  vector[2] prior_mean;
  cov_matrix[2] prior_cov;
}

parameters {
  vector[2] theta;
}

model {
  theta ~ multi_normal(prior_mean, prior_cov);
  for (n in 1:N) {
    y[n] ~ normal(theta, 1);
  }
}
