#' FD global sensitivity to a learning rate
#'
#' Compute global sensitivity when a scalar learning rate multiplies a fixed
#' loss and the prior is unchanged. The prior score cancels and the FD is a
#' quadratic function of the learning rate, so no numerical optimisation is
#' needed.
#'
#' @param fit A `CmdStanMCMC` fit containing reference-posterior draws.
#' @param variables Character vector of model parameters to extract from `fit`.
#' @param lambda_ref Reference learning rate.
#' @param lower,upper Finite scalar endpoints of the candidate learning-rate
#'   interval.
#' @param score_loss Function of the draw matrix returning one gradient
#'   vector of the unscaled loss per posterior draw.
#'
#' @return An object of class `fd_sensitivity_result`.
#' @export
fd_lr_global_sensitivity <- function(
    fit, variables, lambda_ref, lower, upper, score_loss) {
  if (length(lambda_ref) != 1L || !is.finite(lambda_ref)) {
    stop("`lambda_ref` must be a finite scalar.", call. = FALSE)
  }
  if (length(lower) != 1L || length(upper) != 1L ||
      !is.finite(lower) || !is.finite(upper) || lower >= upper) {
    stop("`lower` and `upper` must be finite scalars with lower < upper.", call. = FALSE)
  }

  draws <- as_reference_draws(fit, variables)
  gradients <- validate_scores(score_loss(draws), draws, "score_loss")
  gradient_energy <- mean(rowSums(gradients^2))
  lambda_min <- min(max(lambda_ref, lower), upper)
  endpoints <- c(lower, upper)
  lambda_max <- endpoints[which.max((endpoints - lambda_ref)^2)]
  fd_min <- (lambda_min - lambda_ref)^2 * gradient_energy
  fd_max <- (lambda_max - lambda_ref)^2 * gradient_energy

  result <- new_fd_sensitivity_result(
    fd_max - fd_min, fd_min, fd_max,
    c(lambda = lambda_min), c(lambda = lambda_max),
    list(lower = c(lambda = lower), upper = c(lambda = upper)),
    draws, "learning_rate"
  )
  result$lambda_ref <- lambda_ref
  result$gradient_energy <- gradient_energy
  result
}
