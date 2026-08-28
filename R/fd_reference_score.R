# Derive the reference-posterior score (the gradient of the log density)
# directly from a CmdStanMCMC fit's own log-density machinery, instead of
# requiring the user to supply `score_ref`.
#
# This only works for parameters with unconstrained (identity-transform)
# support: cmdstanr's `grad_log_prob()` differentiates with respect to the
# *unconstrained* parameterization, which equals the constrained one only
# when there is no bounding transform. We verify that before trusting the
# gradient.
fd_reference_score <- function(fit, draws, variables) {
  tryCatch(
    fit$init_model_methods(),
    error = function(e) {
      hint <- if (grepl("pre-compiled", conditionMessage(e), fixed = TRUE)) {
        paste(
          "cmdstanr only exposes a fit's log-density gradient when the model",
          "was compiled in the current R session, not when a cached",
          "executable is reused. Recompile with",
          "`cmdstan_model(..., force_recompile = TRUE)`, or supply",
          "`score_ref` explicitly instead."
        )
      } else {
        paste0(
          "Failed to initialize Stan model methods for the automatic ",
          "`score_ref`: ", conditionMessage(e),
          ". Supply `score_ref` explicitly instead."
        )
      }
      stop(hint, call. = FALSE)
    }
  )

  unconstrained <- as.matrix(
    fit$unconstrain_draws(draws = fit$draws(), format = "draws_matrix")
  )
  # `unconstrain_draws()` does not propagate variable names (columns come
  # back as placeholders like "...1"); recover the flattened parameter order
  # from `variable_skeleton()`, which grad_log_prob() also follows.
  colnames(unconstrained) <- unconstrained_variable_names(fit)

  missing_variables <- setdiff(variables, colnames(unconstrained))
  if (length(missing_variables) > 0L) {
    stop(
      "Could not locate unconstrained draws for: ",
      paste(missing_variables, collapse = ", "),
      call. = FALSE
    )
  }

  if (!isTRUE(all.equal(
    unconstrained[, variables, drop = FALSE],
    draws,
    check.attributes = FALSE,
    tolerance = 1e-6
  ))) {
    stop(
      "Automatic `score_ref` only supports parameters with unconstrained ",
      "support (no lower/upper bounds or other transforms). Supply ",
      "`score_ref` explicitly for: ", paste(variables, collapse = ", "),
      call. = FALSE
    )
  }

  gradients <- matrix(
    NA_real_,
    nrow = nrow(unconstrained),
    ncol = ncol(unconstrained),
    dimnames = list(NULL, colnames(unconstrained))
  )
  for (i in seq_len(nrow(unconstrained))) {
    gradients[i, ] <- as.numeric(
      fit$grad_log_prob(unconstrained_variables = as.numeric(unconstrained[i, ]))
    )
  }

  gradients[, variables, drop = FALSE]
}

# The flattened, in-order parameter names matching `unconstrain_draws()`'s
# column order and `grad_log_prob()`'s input/output order.
unconstrained_variable_names <- function(fit) {
  skeleton <- fit$variable_skeleton()
  unlist(lapply(names(skeleton), function(name) {
    n <- length(skeleton[[name]])
    if (n <= 1L) name else paste0(name, "[", seq_len(n), "]")
  }), use.names = FALSE)
}
