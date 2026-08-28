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

## Learning-rate sensitivity

`fd_lr_global_sensitivity()` computes the learning-rate sensitivity.

## Complete examples

- [`examples/gaussian_location_prior_sensitivity.R`](examples/gaussian_location_prior_sensitivity.R): automatic Gaussian quadratic route;
- [`examples/gaussian_location_prior_sensitivity_multidim.R`](examples/gaussian_location_prior_sensitivity_multidim.R): explicit black-box route;
- [`examples/gaussian_location_lr_sensitivity.R`](examples/gaussian_location_lr_sensitivity.R): learning-rate sensitivity.

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