# fd-sens (R extension)

This is a minimal R package for Fisher-divergence global sensitivity analysis using a reference posterior fitted with Stan.

Complete methodology is described in [*A computationally-tractable measure of global sensitivity for sampling-based Bayesian inference*](https://arxiv.org/abs/2605.28099). Please, cite before usage.

## Installation

```r
install.packages(".", repos = NULL, type = "source")
```

Install `cmdstanr` and CmdStan separately by following the [`cmdstanr` installation guide](https://mc-stan.org/cmdstanr/articles/cmdstanr.html).

## Prior sensitivity

`fd_prior_global_sensitivity()` has two computational routes.

### Exponential-family route

With `method = "auto"` (the default), the package inspects a direct prior statement in the supplied Stan program
and if it belongs to exponential family distributions listed below:

- `normal(mu, sigma)`,
- `gamma(alpha, beta)`,
- `beta(alpha, beta)`.

Use `stan_exponential_family_priors()` to see their natural parameters and sufficient statistics. 

In this case, optimisation is performed using `optimization = "quadratic_corner"` route. See Proposition 2 in the paper.

### Black-box route

If the prior cannot be classified as exponential family safely, `method = "auto"` uses 
the black-box optimisation algorithm `optimization = "black_box"`.

## Sensitivity to independent prior components

If the reference and candidate priors both factorise into disjoint parameter blocks and the candidate region is a Cartesian product across those blocks, set `independent = TRUE` and pass a named `blocks` list instead of `variables`/`lambda_lower`/`lambda_upper`. 
The total minimum, maximum, and sensitivity are just the sums of each block's own minimum, maximum, and sensitivity, so each block is optimised on its own.

Each block is resolved exactly like a standalone call to `fd_prior_global_sensitivity()`.
The result is an `fd_sensitivity_decomposition` (which also inherits `fd_sensitivity_result`): `sensitivity`, `fd_min`, `fd_max`, `lambda_min`, and `lambda_max` are the totals (`lambda_min`/`lambda_max` are named lists, one entry per block),
and `components` is a data frame with one row per block giving its own `sensitivity`, `fd_min`, `fd_max`, and `sensitivity_share`.

## Learning-rate sensitivity

`fd_lr_global_sensitivity()` computes the learning-rate sensitivity.

## Algorithm

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

See [`GETTING_STARTED.md`](GETTING_STARTED.md) for the underlying methodology.

## Complete examples

- [`examples/gaussian_location_prior_sensitivity.R`](examples/gaussian_location_prior_sensitivity.R): automatic Gaussian quadratic route, then interprets the result (see below);
- [`examples/gaussian_location_prior_sensitivity_multidim.R`](examples/gaussian_location_prior_sensitivity_multidim.R): explicit black-box route;
- [`examples/gaussian_location_lr_sensitivity.R`](examples/gaussian_location_lr_sensitivity.R): learning-rate sensitivity;
- [`examples/kilpisjarvi_ar5_independent_prior_sensitivity.R`](examples/kilpisjarvi_ar5_independent_prior_sensitivity.R): `rstan` AR(5) analysis with independent-prior decomposition, then interprets the result (see below).
- [`notebooks/gaussian_location_prior_sensitivity.md`](notebooks/gaussian_location_prior_sensitivity.md) walks through the same Gaussian-location workflow as a narrated notebook, with the printed sensitivity result.

## Interpreting results

[`interpretation/plots.R`](interpretation/plots.R) provides reusable helpers, writing output into a directory:

- `save_sensitivity_result(result, output_dir)` — writes an `fd_sensitivity_result`'s (or `fd_sensitivity_decomposition`'s) values (`sensitivity`, `lambda_min`/`lambda_max`, `fd_min`/`fd_max`, `interval`, `analysis`, ...) as a JSON dict, leaving out draw-level fields at every nesting level;
- `plot_quantiles(fits, variables, output_dir)` — writes a posterior-quantile table (`.csv`) and a median/90%-interval plot (`.png`) comparing a named list of `CmdStanMCMC` or `stanfit` fits (e.g. `list(reference = ref_fit, candidate = cand_fit)`);
- `plot_kde(fits, variables, output_dir)` — writes a kernel density estimate comparison (`.png`) across those fits;
- `plot_ecdf(fits, variables, output_dir)` — writes an empirical CDF comparison (`.png`) across those fits;
- `plot_component_shares(result, output_dir)` — for an `fd_sensitivity_decomposition` (from `independent = TRUE`), writes a 100%-stacked bar (`.png`) showing each block's percentage share of the total sensitivity.

## Contributing

1. For the exponential-family route contributions adding more families are welcome: extend `stan_prior_registry()` in [`R/stan_exponential_family_priors.R`](R/stan_exponential_family_priors.R) with a new entry keyed by the Stan distribution name, providing:

- `n_arguments`: number of arguments in the Stan sampling statement;
- `original_names`, `natural_names`: labels for the original and natural parameters;
- `to_natural(x)`: maps the original parameters to natural parameters `lambda`;
- `gradient(theta)`: the Jacobian of the sufficient statistics with respect to `theta`, one row per draw (used to build the quadratic form's Gram matrix);
- `support(theta)`: a per-draw validity check (e.g. `theta > 0` for a gamma prior);
- `valid_box(lower, upper)`: rejects a natural-parameter box that would leave the family's valid parameter space.

## References
[*A computationally-tractable measure of global sensitivity for sampling-based Bayesian inference*](https://arxiv.org/abs/2605.28099).