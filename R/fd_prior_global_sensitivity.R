#' FD global sensitivity to prior hyperparameters
#'
#' Analyse one prior block, or decompose the analysis across explicitly
#' declared independent prior blocks. Each block uses exact corner
#' maximisation when a supported exponential-family candidate is available and
#' otherwise uses black-box optimisation.
#'
#' @param fit A `CmdStanMCMC` or `stanfit` reference-posterior fit.
#' @param variables Parameters belonging to a single prior block.
#' @param lambda_lower,lambda_upper Candidate-box bounds. On the quadratic
#'   route these must be natural-parameter bounds.
#' @param method One of `"auto"`, `"quadratic"`, or `"black_box"`.
#' @param stan_file Path to the reference Stan program.
#' @param prior_variable Left-hand side of the reference prior statement.
#' @param candidate_family Optional supported Stan family for the candidate
#'   prior. This is needed when the candidate family differs from the reference
#'   family, such as a half-Cauchy reference and inverse-Gamma candidates.
#' @param stan_data Named list resolving scalar arguments in the Stan prior.
#' @param score_prior_ref Optional reference-prior score function. Required for
#'   a black-box route and for a quadratic route whose reference family differs
#'   from the candidate exponential family.
#' @param score_prior_candidate Candidate-prior score function for black-box
#'   optimisation.
#' @param independent If `TRUE`, use the additive FD decomposition across the
#'   explicitly specified `blocks`.
#' @param blocks Named list of independent block specifications. Each block may
#'   contain `variables`, `prior_variable`, `candidate_family`,
#'   `lambda_lower`, `lambda_upper`, `method`, and score callbacks.
#' @param grid_size Approximate total number of black-box grid points per block.
#' @param tol Optimisation tolerance passed to [stats::optim()].
#'
#' @return An `fd_sensitivity_result`, or an
#'   `fd_sensitivity_decomposition` when `independent = TRUE`.
#' @export
fd_prior_global_sensitivity <- function(
    fit,
    variables = NULL,
    lambda_lower = NULL,
    lambda_upper = NULL,
    method = c("auto", "quadratic", "black_box"),
    stan_file = NULL,
    prior_variable = NULL,
    candidate_family = NULL,
    stan_data = list(),
    score_prior_ref = NULL,
    score_prior_candidate = NULL,
    independent = FALSE,
    blocks = NULL,
    grid_size = 201L,
    tol = .Machine$double.eps^0.25) {
  method <- match.arg(method)
  if (isTRUE(independent)) {
    return(fd_prior_independent(
      fit, blocks, stan_file, stan_data, grid_size, tol
    ))
  }
  fd_prior_single(
    fit, variables, lambda_lower, lambda_upper, method, stan_file,
    prior_variable, candidate_family, stan_data, score_prior_ref,
    score_prior_candidate, grid_size, tol
  )
}

fd_prior_single <- function(
    fit, variables, lambda_lower, lambda_upper, method, stan_file,
    prior_variable, candidate_family, stan_data, score_prior_ref,
    score_prior_candidate, grid_size, tol) {
  method <- match.arg(method, c("auto", "quadratic", "black_box"))
  if (is.null(prior_variable) && length(variables) == 1L) {
    prior_variable <- variables
  }
  prior_info <- NULL
  if (method != "black_box") {
    prior_info <- candidate_prior_info(
      stan_file, prior_variable, candidate_family, stan_data
    )
    if (isTRUE(prior_info$supported)) {
      return(fd_prior_quadratic(
        fit, variables, lambda_lower, lambda_upper, prior_info,
        score_prior_ref, tol
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

candidate_prior_info <- function(
    stan_file, prior_variable, candidate_family, stan_data) {
  reference_info <- NULL
  if (!is.null(stan_file) && !is.null(prior_variable)) {
    reference_info <- stan_prior_info(stan_file, prior_variable, stan_data)
  }

  if (is.null(candidate_family)) {
    if (is.null(reference_info)) {
      return(list(
        supported = FALSE,
        reason = "`stan_file` and `prior_variable` were not both supplied."
      ))
    }
    return(reference_info)
  }

  registry <- stan_prior_registry()
  if (length(candidate_family) != 1L || !candidate_family %in% names(registry)) {
    return(list(
      supported = FALSE,
      reason = paste0(
        "Candidate family `", candidate_family,
        "` is not in the supported exponential-family registry."
      )
    ))
  }
  reference_natural <- NULL
  if (!is.null(reference_info) && isTRUE(reference_info$supported) &&
      identical(reference_info$family, candidate_family)) {
    reference_natural <- reference_info$natural_parameters
  }
  list(
    supported = TRUE,
    variable = prior_variable,
    family = candidate_family,
    natural_parameters = reference_natural,
    reference = reference_info,
    reason = "Supported candidate exponential-family prior specified."
  )
}

fd_prior_quadratic <- function(
    fit, variables, lambda_lower, lambda_upper, prior_info,
    score_prior_ref, tol) {
  registry <- stan_prior_registry()[[prior_info$family]]
  bounds <- validate_natural_bounds(
    lambda_lower, lambda_upper, registry, prior_info$family
  )
  lambda_lower <- bounds$lower
  lambda_upper <- bounds$upper
  draws <- as_reference_draws(fit, variables)
  validate_detected_variable(draws, prior_info$variable)

  lambda_ref <- prior_info$natural_parameters
  if (is.function(score_prior_ref)) {
    reference_scores <- validate_scores(
      score_prior_ref(draws), draws, "score_prior_ref"
    )
  } else if (!is.null(lambda_ref)) {
    reference_scores <- sufficient_statistic_scores(draws, registry, lambda_ref)
  } else {
    stop(
      "`score_prior_ref` is required because the detected reference prior ",
      "does not match candidate family `", prior_info$family, "`.",
      call. = FALSE
    )
  }

  terms <- sufficient_statistic_quadratic(draws, registry, reference_scores)
  A <- terms$A
  b <- terms$b
  c_value <- terms$c
  objective <- function(lambda) {
    as.numeric(crossprod(lambda, A %*% lambda) + crossprod(b, lambda) + c_value)
  }

  corners <- as.matrix(expand.grid(lapply(seq_along(lambda_lower), function(i) {
    c(lambda_lower[i], lambda_upper[i])
  })))
  colnames(corners) <- names(lambda_lower)
  corner_fd <- apply(corners, 1L, objective)
  maximum_index <- which.max(corner_fd)
  lambda_max <- corners[maximum_index, ]
  fd_max <- corner_fd[maximum_index]

  unconstrained <- tryCatch(
    as.numeric(-0.5 * qr.solve(A, b)),
    error = function(e) NULL
  )
  if (!is.null(unconstrained) &&
      all(unconstrained >= lambda_lower) && all(unconstrained <= lambda_upper)) {
    lambda_min <- stats::setNames(unconstrained, names(lambda_lower))
    fd_min <- max(0, objective(lambda_min))
  } else {
    start <- if (!is.null(unconstrained)) {
      pmin(pmax(unconstrained, lambda_lower), lambda_upper)
    } else if (!is.null(lambda_ref)) {
      pmin(pmax(lambda_ref, lambda_lower), lambda_upper)
    } else {
      (lambda_lower + lambda_upper) / 2
    }
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
    fd_max - fd_min, fd_min, fd_max, lambda_min, lambda_max,
    list(lower = lambda_lower, upper = lambda_upper), draws, "prior"
  )
  result$optimization <- "quadratic_corner"
  result$prior_family <- prior_info$family
  result$reference_natural_parameters <- lambda_ref
  result$A <- A
  result$b <- b
  result$c <- c_value
  result$corners <- corners
  result$corner_fd <- corner_fd
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
    minimum$objective, maximum$objective,
    minimum$lambda, maximum$lambda,
    list(lower = lower, upper = upper), draws, "prior"
  )
  result$optimization <- "black_box"
  result
}

fd_prior_independent <- function(
    fit, blocks, stan_file, stan_data, grid_size, tol) {
  if (!is.list(blocks) || length(blocks) == 0L) {
    stop("`blocks` must be a non-empty named list.", call. = FALSE)
  }
  if (is.null(names(blocks)) || any(!nzchar(names(blocks))) || anyDuplicated(names(blocks))) {
    stop("`blocks` must have complete, unique names.", call. = FALSE)
  }

  results <- lapply(names(blocks), function(block_name) {
    block <- blocks[[block_name]]
    if (!is.list(block)) stop("Each entry of `blocks` must be a list.", call. = FALSE)
    fd_prior_single(
      fit = fit,
      variables = block$variables,
      lambda_lower = block$lambda_lower,
      lambda_upper = block$lambda_upper,
      method = match.arg(
        block$method %||% "auto", c("auto", "quadratic", "black_box")
      ),
      stan_file = block$stan_file %||% stan_file,
      prior_variable = block$prior_variable,
      candidate_family = block$candidate_family,
      stan_data = block$stan_data %||% stan_data,
      score_prior_ref = block$score_prior_ref,
      score_prior_candidate = block$score_prior_candidate,
      grid_size = block$grid_size %||% grid_size,
      tol = block$tol %||% tol
    )
  })
  names(results) <- names(blocks)

  components <- data.frame(
    block = names(results),
    sensitivity = vapply(results, `[[`, numeric(1L), "sensitivity"),
    fd_min = vapply(results, `[[`, numeric(1L), "fd_min"),
    fd_max = vapply(results, `[[`, numeric(1L), "fd_max"),
    corner_evaluations = vapply(
      results,
      function(x) if (is.null(x$corners)) NA_integer_ else nrow(x$corners),
      integer(1L)
    ),
    stringsAsFactors = FALSE
  )
  total_sensitivity <- sum(components$sensitivity)
  components$sensitivity_share <- if (total_sensitivity > 0) {
    components$sensitivity / total_sensitivity
  } else {
    rep(NA_real_, nrow(components))
  }
  result <- list(
    sensitivity = total_sensitivity,
    fd_min = sum(components$fd_min),
    fd_max = sum(components$fd_max),
    lambda_min = lapply(results, `[[`, "lambda_min"),
    lambda_max = lapply(results, `[[`, "lambda_max"),
    components = components,
    block_results = results,
    corner_evaluations = sum(components$corner_evaluations, na.rm = TRUE),
    analysis = "prior",
    optimization = "independent_decomposition"
  )
  class(result) <- c("fd_sensitivity_decomposition", "fd_sensitivity_result")
  result
}

`%||%` <- function(x, y) if (is.null(x)) y else x

#' @export
print.fd_sensitivity_decomposition <- function(x, ...) {
  cat("FD prior sensitivity: independent-block decomposition\n")
  cat("  sensitivity:", format(x$sensitivity, digits = 6L), "\n")
  cat("  minimum FD: ", format(x$fd_min, digits = 6L), "\n")
  cat("  maximum FD: ", format(x$fd_max, digits = 6L), "\n")
  print(x$components, row.names = FALSE)
  invisible(x)
}
