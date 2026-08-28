#' FD global sensitivity to prior hyperparameters
#'
#' Use exact corner maximisation when a supported exponential-family prior is
#' detected in a Stan program. Otherwise use black-box optimisation with
#' user-supplied prior-score functions.
#'
#' @param fit A `CmdStanMCMC` fit containing reference-posterior draws.
#' @param variables Character vector of model parameters to extract from `fit`.
#' @param lambda_lower,lambda_upper Finite vectors defining the candidate box.
#'   For the quadratic route these are bounds on the prior's **natural
#'   parameters**, in the order returned by [stan_prior_info()]. For the
#'   black-box route they parameterise `score_prior_candidate` directly.
#' @param method One of `"auto"`, `"quadratic"`, or `"black_box"`. `"auto"`
#'   uses the quadratic route only when the Stan prior can be classified
#'   conservatively; otherwise it falls back to black-box optimisation.
#' @param stan_file Path to the Stan program. Required for automatic detection.
#' @param prior_variable Name on the left-hand side of the Stan prior statement,
#'   for example `"theta"` in `theta ~ normal(mu, sigma)`.
#' @param stan_data Named list used to resolve scalar prior arguments appearing
#'   in the Stan statement. Numeric literals need not be included.
#' @param score_prior_ref Function returning the reference-prior score. Required
#'   only for the black-box route.
#' @param score_prior_candidate Function of draws and `lambda` returning the
#'   candidate-prior score. Required only for the black-box route.
#' @param grid_size Approximate total number of black-box grid points.
#' @param tol Optimisation tolerance passed to [stats::optim()].
#'
#' @return An object of class `fd_sensitivity_result`. Quadratic results also
#'   contain `A`, `b`, `c`, the detected family, reference natural parameters,
#'   and all evaluated corners.
#' @export
fd_prior_global_sensitivity <- function(
    fit,
    variables,
    lambda_lower,
    lambda_upper,
    method = c("auto", "quadratic", "black_box"),
    stan_file = NULL,
    prior_variable = NULL,
    stan_data = list(),
    score_prior_ref = NULL,
    score_prior_candidate = NULL,
    grid_size = 201L,
    tol = .Machine$double.eps^0.25) {
  method <- match.arg(method)
  prior_info <- NULL

  if (method != "black_box") {
    if (!is.null(stan_file) && !is.null(prior_variable)) {
      prior_info <- stan_prior_info(stan_file, prior_variable, stan_data)
    } else {
      prior_info <- list(
        supported = FALSE,
        reason = "`stan_file` and `prior_variable` were not both supplied."
      )
    }

    if (isTRUE(prior_info$supported)) {
      return(fd_prior_quadratic(
        fit = fit,
        variables = variables,
        lambda_lower = lambda_lower,
        lambda_upper = lambda_upper,
        prior_info = prior_info,
        tol = tol
      ))
    }
    if (method == "quadratic") {
      stop("Quadratic prior route unavailable: ", prior_info$reason, call. = FALSE)
    }
  }

  if (!is.function(score_prior_ref) || !is.function(score_prior_candidate)) {
    reason <- if (!is.null(prior_info$reason)) paste0(" Detection result: ", prior_info$reason) else ""
    stop(
      "The black-box route requires `score_prior_ref` and ",
      "`score_prior_candidate`.", reason,
      call. = FALSE
    )
  }

  result <- fd_prior_black_box(
    fit, variables, lambda_lower, lambda_upper, score_prior_ref,
    score_prior_candidate, grid_size, tol
  )
  result$detection <- prior_info
  result
}

fd_prior_quadratic <- function(
    fit, variables, lambda_lower, lambda_upper, prior_info, tol) {
  registry <- stan_prior_registry()[[prior_info$family]]
  bounds <- validate_natural_bounds(
    lambda_lower, lambda_upper, registry, prior_info$family
  )
  lambda_lower <- bounds$lower
  lambda_upper <- bounds$upper
  draws <- as_reference_draws(fit, variables)
  validate_detected_variable(draws, prior_info$variable)

  A <- sufficient_statistic_gram(draws, registry)
  lambda_ref <- prior_info$natural_parameters
  b <- as.numeric(-2 * A %*% lambda_ref)
  names(b) <- names(lambda_ref)
  c_value <- as.numeric(crossprod(lambda_ref, A %*% lambda_ref))
  objective <- function(lambda) {
    as.numeric(crossprod(lambda, A %*% lambda) + crossprod(b, lambda) + c_value)
  }

  corner_grid <- as.matrix(expand.grid(lapply(seq_along(lambda_lower), function(i) {
    c(lambda_lower[i], lambda_upper[i])
  })))
  colnames(corner_grid) <- names(lambda_lower)
  corner_values <- apply(corner_grid, 1L, objective)
  maximum_index <- which.max(corner_values)
  lambda_max <- corner_grid[maximum_index, ]
  fd_max <- corner_values[maximum_index]

  if (all(lambda_ref >= lambda_lower) && all(lambda_ref <= lambda_upper)) {
    lambda_min <- lambda_ref
    fd_min <- 0
  } else {
    start <- pmin(pmax(lambda_ref, lambda_lower), lambda_upper)
    minimum <- stats::optim(
      par = start,
      fn = objective,
      method = "L-BFGS-B",
      lower = lambda_lower,
      upper = lambda_upper,
      control = list(pgtol = tol)
    )
    lambda_min <- minimum$par
    fd_min <- max(0, minimum$value)
  }

  result <- new_fd_sensitivity_result(
    fd_max - fd_min,
    fd_min,
    fd_max,
    lambda_min,
    lambda_max,
    list(lower = lambda_lower, upper = lambda_upper),
    draws,
    "prior"
  )
  result$optimization <- "quadratic_corner"
  result$prior_family <- prior_info$family
  result$reference_natural_parameters <- lambda_ref
  result$A <- A
  result$b <- b
  result$c <- c_value
  result$corners <- corner_grid
  result$corner_fd <- corner_values
  result$detection <- prior_info
  result
}

fd_prior_black_box <- function(
    fit, variables, lambda_lower, lambda_upper, score_prior_ref,
    score_prior_candidate, grid_size, tol) {
  bounds <- validate_bounds(lambda_lower, lambda_upper)
  lower <- bounds$lower
  upper <- bounds$upper
  param_names <- names(lower)
  if (length(grid_size) != 1L || !is.finite(grid_size) || grid_size < 3L) {
    stop("`grid_size` must be an integer of at least 3.", call. = FALSE)
  }

  k <- length(lower)
  points_per_dim <- max(3L, floor(grid_size^(1 / k)))
  draws <- as_reference_draws(fit, variables)
  ref_scores <- validate_scores(score_prior_ref(draws), draws, "score_prior_ref")
  objective <- function(lambda) {
    candidate_scores <- validate_scores(
      score_prior_candidate(draws, lambda), draws, "score_prior_candidate"
    )
    mean(rowSums((ref_scores - candidate_scores)^2))
  }

  dims <- lapply(seq_len(k), function(i) {
    seq(lower[i], upper[i], length.out = points_per_dim)
  })
  grid <- as.matrix(expand.grid(dims))
  colnames(grid) <- param_names
  values <- apply(grid, 1L, objective)
  minimum <- refine_extreme(
    grid, values, objective, lower, upper, maximum = FALSE, tol = tol
  )
  maximum <- refine_extreme(
    grid, values, objective, lower, upper, maximum = TRUE, tol = tol
  )

  result <- new_fd_sensitivity_result(
    maximum$objective - minimum$objective,
    minimum$objective,
    maximum$objective,
    minimum$lambda,
    maximum$lambda,
    list(lower = lower, upper = upper),
    draws,
    "prior"
  )
  result$optimization <- "black_box"
  result
}
