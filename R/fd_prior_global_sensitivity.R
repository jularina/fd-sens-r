#' FD global sensitivity to prior hyperparameters
#'
#' Estimate the minimum and maximum Fisher divergence when the prior varies
#' over a hyperrectangle and the likelihood or loss is kept fixed. Because the
#' common likelihood score cancels, only reference and candidate prior scores
#' are required. Searches with a factorial grid over the hyperrectangle
#' (always including every corner), then bounded quasi-Newton refinement of
#' each extremum found on that grid.
#'
#' @param fit A `CmdStanMCMC` fit containing reference-posterior draws.
#' @param variables Character vector of model parameters to extract from `fit`.
#' @param lower,upper Finite numeric vectors defining the candidate prior-
#'   hyperparameter box. Named vectors are recommended.
#' @param score_prior_ref Function of the draw matrix returning the score of
#'   the reference prior.
#' @param score_prior_candidate Function of the draw matrix and prior-
#'   hyperparameter vector `lambda` returning the candidate prior score.
#' @param grid_size Approximate total number of grid points used to locate the
#'   global extrema before local refinement, with at least 3 per dimension.
#' @param tol Gradient-projection tolerance passed to [stats::optim()].
#'
#' @return An object of class `fd_sensitivity_result`.
#' @export
fd_prior_global_sensitivity <- function(
    fit, variables, lower, upper, score_prior_ref, score_prior_candidate,
    grid_size = 201L, tol = .Machine$double.eps^0.25) {
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

  new_fd_sensitivity_result(
    maximum$objective - minimum$objective,
    minimum$objective,
    maximum$objective,
    minimum$lambda,
    maximum$lambda,
    list(lower = lower, upper = upper),
    draws,
    "prior"
  )
}
