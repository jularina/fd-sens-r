# fd-sens (R extension)

This is a minimal R package for Fisher-divergence global sensitivity analysis using a reference posterior fitted with Stan.

## Installation

```r
install.packages(".", repos = NULL, type = "source") # install fdsens
```

Install `cmdstanr` and CmdStan separately by following the [`cmdstanr` installation guide](https://mc-stan.org/cmdstanr/articles/cmdstanr.html).

## Complete examples

See [`examples/gaussian_location_prior_sensitivity.R`](examples/gaussian_location_prior_sensitivity.R), [`examples/gaussian_location_prior_multidim.R`](examples/gaussian_location_prior_multidim.R) for prior sensitivity.
See [`examples/gaussian_location_lr_sensitivity.R`](examples/gaussian_location_lr_sensitivity.R) for learning rate sensitivity. 

