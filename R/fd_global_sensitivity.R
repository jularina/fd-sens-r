#' FD global sensitivity for one or more hyperparameters
#'
#' Estimate the minimum and maximum Fisher divergence over a hyperrectangle of
#' candidate hyperparameter values. The same reference-posterior draws are
#' reused for every candidate.
#'
#' @param fit A `CmdStanMCMC` fit.
#' @param variables Character vector of parameters to extract from `fit`.
#' @param lower,upper Finite numeric vectors of the same length, one entry per
#'   hyperparameter dimension, with `lower < upper` elementwise. A length-1
#'   `lower`/`upper` recovers the original single-hyperparameter search. Name
#'   the entries (e.g. `c(prior_mean = -2, prior_sd = 0.5)`) to address
#'   hyperparameters by name instead of position; `upper` must then use the
#'   same names as `lower` (any order). `score_candidate`'s `lambda` argument
#'   carries the same names, e.g. `lambda["prior_mean"]`.
#' @param score_candidate Function of the draw matrix and the hyperparameter
#'   vector `lambda` (same length, and names if any, as `lower`/`upper`)
#'   returning candidate scores.
#' @param score_ref Optional function of the draw matrix returning
#'   reference-posterior scores. When omitted, the reference score is derived
#'   automatically from `fit`'s log-density gradient, which requires
#'   `variables` to have unconstrained support.
#' @param grid_size Approximate total number of grid points used to locate the
#'   global extrema before local refinement, spread evenly across
#'   hyperparameter dimensions (at least 3 per dimension). Grows
#'   combinatorially with the number of hyperparameters.
#' @param tol Gradient-projection tolerance (`pgtol`) passed to
#'   [stats::optim()]'s `L-BFGS-B` method.
#'
#' @return An object of class `fd_sensitivity_result` containing the estimated
#'   global sensitivity, extrema, hyperparameters, and reference draws.
#' @export
fd_global_sensitivity <- function(
    fit,
    variables,
    lower,
    upper,
    score_candidate,
    score_ref = NULL,
    grid_size = 201L,
    tol = .Machine$double.eps^0.25) {
  bounds <- validate_bounds(lower, upper)
  lower <- bounds$lower
  upper <- bounds$upper
  param_names <- names(lower)

  if (length(grid_size) != 1L || !is.finite(grid_size) || grid_size < 3L) {
    stop("`grid_size` must be an integer of at least 3.", call. = FALSE)
  }
  k <- length(lower)
  points_per_dim <- max(3L, floor(grid_size^(1 / k)))

  draws <- as_reference_draws(fit, variables)
  ref_scores <- if (is.null(score_ref)) {
    fd_reference_score(fit, draws, variables)
  } else {
    validate_scores(score_ref(draws), draws, "score_ref")
  }

  objective <- function(lambda) {
    candidate_scores <- validate_scores(
      score_candidate(draws, lambda),
      draws,
      "score_candidate"
    )
    mean(rowSums((ref_scores - candidate_scores)^2))
  }

  # A factorial grid over each dimension's endpoint-inclusive sequence always
  # contains every corner of the hyperrectangle. For k = 1 this is exactly
  # the original evenly spaced 1-D grid.
  dims <- lapply(seq_len(k), function(i) seq(lower[i], upper[i], length.out = points_per_dim))
  grid <- as.matrix(expand.grid(dims))
  # `expand.grid()` names unnamed columns "Var1", "Var2", ...; only keep
  # column names when the caller actually named `lower`/`upper`.
  colnames(grid) <- param_names
  values <- apply(grid, 1L, objective)

  minimum <- refine_extreme(grid, values, objective, lower, upper, maximum = FALSE, tol = tol)
  maximum <- refine_extreme(grid, values, objective, lower, upper, maximum = TRUE, tol = tol)

  result <- list(
    sensitivity = maximum$objective - minimum$objective,
    fd_min = minimum$objective,
    fd_max = maximum$objective,
    lambda_min = minimum$lambda,
    lambda_max = maximum$lambda,
    interval = list(lower = lower, upper = upper),
    draws = draws
  )
  class(result) <- "fd_sensitivity_result"
  result
}

validate_bounds <- function(lower, upper) {
  if (!is.numeric(lower) || !is.numeric(upper)) {
    stop("`lower` and `upper` must be numeric vectors.", call. = FALSE)
  }
  if (length(lower) == 0L || length(lower) != length(upper)) {
    stop("`lower` and `upper` must be the same, non-zero, length.", call. = FALSE)
  }

  lower_names <- names(lower)
  upper_names <- names(upper)
  named <- !is.null(lower_names) || !is.null(upper_names)
  if (named) {
    if (is.null(lower_names) || is.null(upper_names) ||
        any(!nzchar(lower_names)) || any(!nzchar(upper_names)) ||
        anyDuplicated(lower_names) || anyDuplicated(upper_names)) {
      stop(
        "`lower` and `upper` must both have complete, unique names, or ",
        "both be unnamed.",
        call. = FALSE
      )
    }
    if (!setequal(lower_names, upper_names)) {
      stop("`lower` and `upper` must use the same hyperparameter names.", call. = FALSE)
    }
    upper <- upper[lower_names]
  }

  if (any(!is.finite(lower)) || any(!is.finite(upper)) || any(lower >= upper)) {
    stop(
      "`lower` and `upper` must be finite, with `lower < upper` elementwise.",
      call. = FALSE
    )
  }

  list(lower = lower, upper = upper)
}

refine_extreme <- function(grid, values, objective, lower, upper, maximum, tol) {
  index <- if (maximum) which.max(values) else which.min(values)
  start <- grid[index, ]

  objective_signed <- if (maximum) function(lambda) -objective(lambda) else objective
  refined <- stats::optim(
    par = start,
    fn = objective_signed,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(pgtol = tol)
  )
  refined_value <- if (maximum) -refined$value else refined$value

  # optim() can fail to improve on the coarse grid point (e.g. when the true
  # extremum sits exactly at a corner already in the grid); keep whichever
  # candidate is actually better.
  best_is_refined <- if (maximum) {
    refined_value > values[index]
  } else {
    refined_value < values[index]
  }

  if (best_is_refined) {
    list(lambda = refined$par, objective = refined_value)
  } else {
    list(lambda = start, objective = values[index])
  }
}

#' @export
print.fd_sensitivity_result <- function(x, ...) {
  format_lambda <- function(lambda) {
    values <- vapply(lambda, format, character(1L), digits = 6L)
    if (!is.null(names(lambda))) {
      values <- paste0(names(lambda), " = ", values)
    }
    if (length(values) == 1L) values else paste0("(", paste(values, collapse = ", "), ")")
  }
  cat("FD global sensitivity\n")
  cat("  sensitivity:", format(x$sensitivity, digits = 6L), "\n")
  cat("  minimum FD: ", format(x$fd_min, digits = 6L),
      "at lambda =", format_lambda(x$lambda_min), "\n")
  cat("  maximum FD: ", format(x$fd_max, digits = 6L),
      "at lambda =", format_lambda(x$lambda_max), "\n")
  invisible(x)
}
