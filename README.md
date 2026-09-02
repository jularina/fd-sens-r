# fd-sens (R/Stan-based)

This is a minimal R package for Fisher-divergence global sensitivity analysis using a reference posterior fitted with Stan.

Complete methodology is described in [*A computationally-tractable measure of global sensitivity for sampling-based Bayesian inference*](https://arxiv.org/abs/2605.28099).

See [`GETTING_STARTED.md`](GETTING_STARTED.md) for the underlying methodology.

## Installation

This is not a package published on CRAN, so you first need to clone the repository and then install it from that local copy.

1. Clone the repository and move into it:

```sh
git clone https://github.com/jularina/fd-sens-r.git
cd fd-sens-r
```

2. From inside the `fd-sens-r` directory, install the package from source. In R, your working directory must be set to `fd-sens-r` (the folder containing `DESCRIPTION`) — e.g. open `fd-sens-r` as an RStudio project, or run `setwd("path/to/fd-sens-r")` — then run:

```r
install.packages(".", repos = NULL, type = "source")
```

Install `cmdstanr` and CmdStan separately by following the [`cmdstanr` installation guide](https://mc-stan.org/cmdstanr/articles/cmdstanr.html).

## Quickstart

A Gaussian location model, fit once with `cmdstanr`:

```stan
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
```

Measuring how sensitive its posterior is to the prior — over a range of plausible prior means/scales — is one function call:

```r
library(cmdstanr)
library(fdsens)

set.seed(123)
stan_file <- system.file("stan", "gaussian_location.stan", package = "fdsens")
fit <- cmdstan_model(stan_file)$sample(
  data = list(N = 30, y = rnorm(30, mean = 1), prior_mean = 0, prior_sd = 2),
  seed = 123, chains = 2, parallel_chains = 2,
  iter_warmup = 500, iter_sampling = 1000, refresh = 0
)

result <- fd_prior_global_sensitivity(
  fit, variables = "theta",
  lambda_lower = c(eta1 = -1, eta2 = -2), lambda_upper = c(eta1 = 1, eta2 = -0.05),
  stan_file = stan_file, prior_variable = "theta",
  stan_data = list(prior_mean = 0, prior_sd = 2)
)
print(result)
#> FD prior sensitivity
#>   optimisation: quadratic_corner
#>   sensitivity: 20.976
#>   minimum FD:  0 at lambda = (eta1 = 0, eta2 = -0.125)
#>   maximum FD:  20.976 at lambda = (eta1 = -1, eta2 = -2)
```

(`eta1`/`eta2` are the normal prior's natural parameters; see [Exponential-family route](#exponential-family-route) below.)

That single call searched every candidate prior in the declared box and returned `sensitivity`: the largest possible change in the posterior, together with `lambda_max`, the worst-case prior that produces it. Refitting at `lambda_max` and comparing the two posteriors confirms it — the candidate visibly shifts and tightens `theta`:

<p align="center">
  <img src="notebooks/gaussian_location_prior_sensitivity_files/figure-gfm/compare-kde-1.png" width="45%" alt="Kernel density comparison of reference vs. worst-case posterior for theta">
  <img src="notebooks/gaussian_location_prior_sensitivity_files/figure-gfm/compare-ecdf-1.png" width="45%" alt="Empirical CDF comparison of reference vs. worst-case posterior for theta">
</p>

The full runnable script (including the refit and the plotting code) is at [`examples/gaussian_location_prior_sensitivity.R`](examples/gaussian_location_prior_sensitivity.R), and [`notebooks/gaussian_location_prior_sensitivity.md`](notebooks/gaussian_location_prior_sensitivity.md) walks through it step by step.

## Repository Contents

| Path | What's there |
| --- | --- |
| [`R/`](R) | The package's core functions: `fd_prior_global_sensitivity()`, `fd_lr_global_sensitivity()`, the underlying FD estimator in `fd_sensitivity.R`, the exponential-family prior registry in `stan_exponential_family_priors.R`, and supporting helpers. |
| [`inst/stan/`](inst/stan) | The `.stan` reference models used by the examples and tests (e.g. `gaussian_location.stan`). |
| [`examples/`](examples) | Runnable, reproducible end-to-end scripts (`Rscript examples/<name>.R`) covering prior sensitivity, learning-rate sensitivity, a multidimensional prior, and sensitivity decomposition over independent prior blocks. |
| [`notebooks/`](notebooks) | A rendered walkthrough of the simplest example, step by step with output and plots (`.Rmd` source and its knitted `.md`). |
| [`interpretation/`](interpretation) | `plots.R`: helpers to turn an `fd_sensitivity_result` into tables and plots (posterior quantiles, KDE, ECDF, component-share bar chart, a JSON summary); `output/` is where the examples write their generated files. |
| [`tests/testthat/`](tests/testthat) | Unit tests for the core functions. |

## Functionality

### Prior Sensitivity

`fd_prior_global_sensitivity()` needs a few arguments to be specified:

- `fit`: the reference-posterior fit — a `CmdStanMCMC` or `stanfit` object obtained by sampling your Stan model (e.g. `cmdstan_model(stan_file)$sample(...)`).
- `variables`: the name(s) of the parameter(s), as they appear in `fit`'s draws, whose sensitivity to the prior you're measuring (e.g. `"theta"`).
- `lambda_lower` / `lambda_upper`: the box of candidate hyperparameters to search over. 

`fd_prior_global_sensitivity()` has two computational routes: exponential-family (check if the prior belongs to an [**exponential family**](https://en.wikipedia.org/wiki/Exponential_family)) or black-box.

The exponential-family route additionally needs:

- `stan_file`: path to the reference Stan program.
- `prior_variable`: the left-hand side of that prior statement, e.g. `"theta"` in `theta ~ normal(prior_mean, prior_sd)`. Defaults to `variables` when there's a single variable.
- `stan_data`: a named list resolving any fixed scalar arguments of that statement (e.g. `list(prior_mean = 0, prior_sd = 2)`), so the package can read off the reference prior's own hyperparameters.

The black-box route needs `score_prior_ref(draws)` and `score_prior_candidate(draws, lambda)`: functions you supply that return the gradient of the log-prior density with respect to the parameters, evaluated at the reference and at a candidate `lambda` respectively.

#### Exponential-family route

With `method = "auto"` (the default), the package inspects a direct prior statement in the supplied Stan program
and if it belongs to supported exponential family distributions listed below:

- `normal(mu, sigma)`,
- `gamma(alpha, beta)`,
- `beta(alpha, beta)`.

In this case, optimisation is performed using `optimization = "quadratic_corner"` route. 
If you are sure that the distribution is exponential family and is in the supported distributions - 
directly pass `method = "quadratic"`.

#### Black-box route

If the prior cannot be classified as exponential family safely, `method = "auto"` uses 
the black-box optimisation algorithm `optimization = "black_box"`. If you are sure that the distribution isn't an exponential family or is not in the supported distributions - 
directly pass `method = "black_box"`.

### Sensitivity to Independent Prior Components

If the reference and candidate priors both factorise into disjoint parameter blocks and the candidate region is a Cartesian product across those blocks, set `independent = TRUE` and pass a named `blocks` list instead of `variables`/`lambda_lower`/`lambda_upper`. 
The total minimum, maximum, and sensitivity are just the sums of each block's own minimum, maximum, and sensitivity, so each block is optimised on its own.

Each block is resolved exactly like a standalone call to `fd_prior_global_sensitivity()`.
The result is an `fd_sensitivity_decomposition` (which also inherits `fd_sensitivity_result`): `sensitivity`, `fd_min`, `fd_max`, `lambda_min`, and `lambda_max` are the totals (`lambda_min`/`lambda_max` are named lists, one entry per block),
and `components` is a data frame with one row per block giving its own `sensitivity`, `fd_min`, `fd_max`, and `sensitivity_share`.

### Learning-Rate Sensitivity

`fd_lr_global_sensitivity()` computes the learning-rate sensitivity.

### Algorithm

1. **Prepare a reference Bayesian model as a `.stan` file** (prior + likelihood/loss). Put it e.g. under `inst/stan/` and fit it with `cmdstanr::cmdstan_model()$sample()` to obtain the reference posterior draws (`CmdStanMCMC`).
2. **Choose the sensitivity analysis route**: prior sensitivity (`fd_prior_global_sensitivity()`) or learning-rate sensitivity (`fd_lr_global_sensitivity()`).
3. **If prior sensitivity:**
   - Choose whether the prior is exponential family. With `method = "auto"` (default) this is detected automatically. Choose `method = "quadratic"` when you are sure that it is exponential family, `method = "black_box"` - when not.
   - When `method = "auto"` or `method = "black_box"`, supply score functions `score_prior_ref(draws)` and `score_prior_candidate(draws, lambda)`, each returning a numeric matrix of the same shape as `draws` (one row per posterior draw, one column per variable) holding the gradient of the log-prior density with respect to the parameters. 
   - Choose candidate prior hyperparameters $\lambda$ range - the box `lambda_lower` / `lambda_upper`.
   - Run `fd_prior_global_sensitivity(...)`.
4. **If learning-rate sensitivity:** 
   - Choose `lambda_ref` and the `lower` / `upper` interval.
   - Supply `score_loss`.
   - Run `fd_lr_global_sensitivity(...)`.
5. **Receive results**: an `fd_sensitivity_result` object (print it, or access fields directly):
   - `fd_sensitivity_result$sensitivity` — the global sensitivity value $\widehat{S}_m^{\mathrm{FD}}(\Gamma)$, i.e. `fd_max - fd_min`;
   - `fd_sensitivity_result$lambda_max` — the worst-case hyperparameters $\lambda_{\sup}$, leading to largest posterior change;
   - `fd_sensitivity_result$lambda_min` — the least-sensitive hyperparameters $\lambda_{\inf}$, leading to smallest posterior change;
   - `fd_sensitivity_result$fd_max` / `fd_sensitivity_result$fd_min` — the FD estimates at `lambda_max` / `lambda_min`;
   - `fd_sensitivity_result$interval` — the searched box (`lower`, `upper`) that was passed in;
   - `fd_sensitivity_result$analysis` — `"prior"` or `"learning_rate"`;
   - `fd_sensitivity_result$draws` — the reference-posterior draws used for the estimate.

### Examples

- [`examples/gaussian_location_prior_sensitivity.R`](examples/gaussian_location_prior_sensitivity.R): prior sensitivity for Gaussian location model (simplest example);
- [`notebooks/gaussian_location_prior_sensitivity.md`](notebooks/gaussian_location_prior_sensitivity.md) walks through the same Gaussian-location as above.
- [`examples/gaussian_location_prior_sensitivity_multidim.R`](examples/gaussian_location_prior_sensitivity_multidim.R): prior sensitivity for multidimensional Gaussian location model;
- [`examples/gaussian_location_lr_sensitivity.R`](examples/gaussian_location_lr_sensitivity.R): learning-rate sensitivity;
- [`examples/kilpisjarvi_ar5_independent_prior_sensitivity.R`](examples/kilpisjarvi_ar5_independent_prior_sensitivity.R): `rstan` sensitivity to independent prior components.

Run any example file from the command line with `Rscript`, e.g.:

```sh
Rscript examples/gaussian_location_prior_sensitivity.R
```

### Interpreting Results

[`interpretation/plots.R`](interpretation/plots.R) provides reusable helpers to interpret the results:

- `save_sensitivity_result(result, output_dir)` — writes an `fd_sensitivity_result`'s (or `fd_sensitivity_decomposition`'s) values (`sensitivity`, `lambda_min`/`lambda_max`, `fd_min`/`fd_max`, `interval`, `analysis`, ...) as a JSON dict, leaving out draw-level fields at every nesting level;
- `plot_quantiles(fits, variables, output_dir)` — writes a posterior-quantile table (`.csv`) and a median/90%-interval plot (`.png`) comparing a named list of `CmdStanMCMC` or `stanfit` fits (e.g. `list(reference = ref_fit, candidate = cand_fit)`);
- `plot_kde(fits, variables, output_dir)` — writes a kernel density estimate comparison (`.png`) across those fits;
- `plot_ecdf(fits, variables, output_dir)` — writes an empirical CDF comparison (`.png`) across those fits;
- `plot_component_shares(result, output_dir)` — for an `fd_sensitivity_decomposition`, writes a 100%-stacked bar (`.png`) showing each block's percentage share of the total sensitivity.

## Contributing

1. For the exponential-family route contributions adding more families are welcome: extend `stan_prior_registry()` in [`R/stan_exponential_family_priors.R`](R/stan_exponential_family_priors.R) with a new entry keyed by the Stan distribution name, providing:

- `n_arguments`: number of arguments in the Stan sampling statement;
- `original_names`, `natural_names`: labels for the original and natural parameters;
- `to_natural(x)`: maps the original parameters to natural parameters `lambda`;
- `gradient(theta)`: the Jacobian of the sufficient statistics with respect to `theta`, one row per draw (used to build the quadratic form's Gram matrix);
- `support(theta)`: a per-draw validity check (e.g. `theta > 0` for a gamma prior);
- `valid_box(lower, upper)`: rejects a natural-parameter box that would leave the family's valid parameter space.

## Citing FD-Sens

If you use this package, please cite:

Odnoblyudova A., Dellaporta C., and Briol F.-X. (2026). A computationally-tractable measure of global sensitivity for sampling-based Bayesian inference. [arXiv:2605.28099](https://arxiv.org/abs/2605.28099).