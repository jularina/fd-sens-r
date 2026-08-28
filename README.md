# fdsens

This is a minimal R prototype for Fisher-divergence global sensitivity analysis using a reference posterior fitted with Stan.

The package separates the two parametric sensitivity settings from the paper:

- `fd_lr_global_sensitivity()` varies a scalar learning rate while keeping the prior fixed;
- `fd_prior_global_sensitivity()` varies prior hyperparameters while keeping the likelihood or loss fixed.

Both functions reuse one set of reference-posterior draws. Candidate posteriors are not sampled separately.

## Installation

```r
install.packages(".", repos = NULL, type = "source")
```

Install `cmdstanr` and CmdStan separately by following the [`cmdstanr` installation guide](https://mc-stan.org/cmdstanr/articles/cmdstanr.html).

## Learning-rate sensitivity

Suppose the posterior has loss $\lambda l(\theta;x)$ and reference learning rate $\lambda_{\mathrm{ref}}$. With the prior fixed, the score difference is determined entirely by $\nabla_\theta l$.

```r
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
```

The estimated FD at a candidate $\lambda$ is

$$
(\lambda-\lambda_{\mathrm{ref}})^2
\frac{1}{m}\sum_{i=1}^m
\|\nabla_\theta l(\theta_i;x)\|^2.
$$

The minimum and maximum over the interval are computed analytically; no numerical optimiser is needed.

## Prior sensitivity

When only the prior changes, the common likelihood or loss score cancels. The user supplies the reference and candidate **prior** scores.

```r
score_prior_ref <- function(draws) {
  theta <- draws[, "theta"]
  matrix((prior_mean_ref - theta) / prior_sd_ref^2, ncol = 1)
}

score_prior_candidate <- function(draws, lambda) {
  theta <- draws[, "theta"]
  matrix(
    (lambda["prior_mean"] - theta) / lambda["prior_sd"]^2,
    ncol = 1
  )
}

prior_result <- fd_prior_global_sensitivity(
  fit = fit,
  variables = "theta",
  lower = c(prior_mean = -2, prior_sd = 0.5),
  upper = c(prior_mean = 2, prior_sd = 4),
  score_prior_ref = score_prior_ref,
  score_prior_candidate = score_prior_candidate
)
```

`fd_prior_global_sensitivity()` searches the hyperparameter box using a factorial grid followed by bounded local refinement. Its output contains the estimated global sensitivity, minimum and maximum FD, and the corresponding prior hyperparameters.

## Complete examples

See [`examples/gaussian_location_prior_sensitivity.R`](examples/gaussian_location_prior_sensitivity.R) and [`examples/gaussian_location_lr_sensitivity.R`](examples/gaussian_location_lr_sensitivity.R), each a standalone script fitting the same Stan reference model before running its analysis.

## Current scope

The package currently supports `CmdStanMCMC` fits and explicit R functions for loss gradients or prior scores. Keeping these functions explicit ensures that FD is computed in the model parameterisation intended by the user rather than automatically in Stan's internal unconstrained parameterisation.
