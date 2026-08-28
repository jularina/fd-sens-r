# fdsens

This is a minimal R prototype for Fisher-divergence global sensitivity analysis with a reference posterior fitted in Stan.

The first version supports:

- reference draws from a `cmdstanr` `CmdStanMCMC` fit;
- a reference-posterior score derived automatically from the fit's own
  log-density gradient (for unconstrained parameters), or supplied by hand;
- a user-supplied candidate posterior score function;
- FD estimation from one set of reference-posterior draws (`fd_sensitivity()`);
- global sensitivity for one or more hyperparameters over a closed
  hyperrectangle (`fd_global_sensitivity()`).

## Installation

Install the local prototype with:

```r
install.packages(".", repos = NULL, type = "source")
```

For the Stan example, install `cmdstanr` and CmdStan separately following the [`cmdstanr` installation guide](https://mc-stan.org/cmdstanr/articles/cmdstanr.html).

## Minimal workflow

```r
result <- fd_global_sensitivity(
  fit = fit,
  variables = "theta",
  lower = -2,
  upper = 2,
  score_candidate = score_candidate
)

print(result)
```

`score_candidate` receives the reference draw matrix and the proposed hyperparameter vector `lambda` (same length as `lower`/`upper`; length 1 for a single hyperparameter), and must return one posterior-score vector per draw. `score_ref` is optional: when omitted, the reference score is derived automatically from `fit`'s log-density gradient, which requires `variables` to have unconstrained support (no lower/upper bounds or other transforms). Supply `score_ref` explicitly for constrained parameters.

Automatic `score_ref` also requires `cmdstan_model(..., force_recompile = TRUE)`: cmdstanr only exposes a fit's log-density gradient when the model was compiled in the current R session, not when a cached executable is reused.

`lower`/`upper` may have any length `k`: `fd_global_sensitivity()` searches a `k`-dimensional hyperrectangle (a factorial grid, then local refinement via `stats::optim()`'s `L-BFGS-B`). The grid grows combinatorially with `k`, so `grid_size` is treated as an approximate total, not a per-dimension count.

See [`examples/cmdstanr_gaussian_location.R`](examples/cmdstanr_gaussian_location.R) for a complete Stan example searching jointly over a prior mean and sd.
