# A minimal stand-in for a CmdStanMCMC fit, exposing just the `$`-accessible
# methods fdsens calls (`draws()`, `init_model_methods()`,
# `unconstrain_draws()`, `variable_skeleton()`, `grad_log_prob()`). `draws`
# must have unconstrained support so `unconstrain_draws()` can be the
# identity.
mock_cmdstan_fit <- function(draws, grad_fun) {
  skeleton <- stats::setNames(
    replicate(ncol(draws), 0, simplify = FALSE),
    colnames(draws)
  )

  structure(
    list(
      draws = function(variables = NULL, format = "matrix", ...) {
        if (!is.null(variables)) draws[, variables, drop = FALSE] else draws
      },
      init_model_methods = function(...) invisible(NULL),
      unconstrain_draws = function(draws, format = "draws_matrix", ...) draws,
      variable_skeleton = function(...) skeleton,
      grad_log_prob = function(unconstrained_variables, ...) grad_fun(unconstrained_variables)
    ),
    class = "CmdStanMCMC"
  )
}
