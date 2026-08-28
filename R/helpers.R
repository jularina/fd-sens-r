# Internal helpers shared by fd_sensitivity(), fd_prior_global_sensitivity(),
# and fd_lr_global_sensitivity(). Nothing here is exported except the
# print() method, which every fd_sensitivity_result carries regardless of
# which analysis produced it.

validate_scores <- function(scores, draws, name) {
  if (is.vector(scores) && ncol(draws) == 1L && length(scores) == nrow(draws)) {
    scores <- matrix(scores, ncol = 1L)
  }

  if (!is.matrix(scores) || !is.numeric(scores)) {
    stop("`", name, "` must return a numeric matrix.", call. = FALSE)
  }
  if (!identical(dim(scores), dim(draws))) {
    stop(
      "`", name, "` must return one score vector per draw, with dimensions ",
      nrow(draws), " by ", ncol(draws), ".",
      call. = FALSE
    )
  }
  if (any(!is.finite(scores))) {
    stop("`", name, "` returned non-finite scores.", call. = FALSE)
  }

  scores
}

validate_bounds <- function(lower, upper) {
  if (!is.numeric(lower) || !is.numeric(upper)) {
    stop("`lower` and `upper` must be numeric vectors.", call. = FALSE)
  }
  if (length(lower) == 0L || length(lower) != length(upper)) {
    stop("`lower` and `upper` must have the same non-zero length.", call. = FALSE)
  }

  lower_names <- names(lower)
  upper_names <- names(upper)
  named <- !is.null(lower_names) || !is.null(upper_names)
  if (named) {
    if (is.null(lower_names) || is.null(upper_names) ||
        any(!nzchar(lower_names)) || any(!nzchar(upper_names)) ||
        anyDuplicated(lower_names) || anyDuplicated(upper_names)) {
      stop(
        "`lower` and `upper` must both have complete, unique names, or both be unnamed.",
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
    par = start, fn = objective_signed, method = "L-BFGS-B",
    lower = lower, upper = upper, control = list(pgtol = tol)
  )
  refined_value <- if (maximum) -refined$value else refined$value
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

new_fd_sensitivity_result <- function(
    sensitivity, fd_min, fd_max, lambda_min, lambda_max, interval, draws,
    analysis) {
  result <- list(
    sensitivity = sensitivity, fd_min = fd_min, fd_max = fd_max,
    lambda_min = lambda_min, lambda_max = lambda_max, interval = interval,
    draws = draws, analysis = analysis
  )
  class(result) <- "fd_sensitivity_result"
  result
}

#' @export
print.fd_sensitivity_result <- function(x, ...) {
  format_lambda <- function(lambda) {
    values <- vapply(lambda, format, character(1L), digits = 6L)
    if (!is.null(names(lambda))) values <- paste0(names(lambda), " = ", values)
    if (length(values) == 1L) values else paste0("(", paste(values, collapse = ", "), ")")
  }
  label <- if (identical(x$analysis, "learning_rate")) {
    "FD learning-rate sensitivity"
  } else {
    "FD prior sensitivity"
  }
  cat(label, "\n")
  cat("  sensitivity:", format(x$sensitivity, digits = 6L), "\n")
  cat("  minimum FD: ", format(x$fd_min, digits = 6L),
      "at lambda =", format_lambda(x$lambda_min), "\n")
  cat("  maximum FD: ", format(x$fd_max, digits = 6L),
      "at lambda =", format_lambda(x$lambda_max), "\n")
  invisible(x)
}
