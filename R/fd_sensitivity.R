#' Estimate the Fisher divergence between a fit's posterior and a candidate
#'
#' Both score functions must return an `m` by `d` numeric matrix, where `m` is
#' the number of posterior draws and `d` is the dimension of the posterior.
#'
#' @param fit A `CmdStanMCMC` fit.
#' @param variables Character vector naming the model parameters to extract.
#' @param score_candidate Function of the draw matrix and `lambda` returning
#'   candidate-posterior scores.
#' @param lambda Candidate hyperparameter value.
#' @param score_ref Optional function of the draw matrix returning
#'   reference-posterior scores. When omitted, the reference score is derived
#'   automatically from `fit`'s log-density gradient, which requires
#'   `variables` to have unconstrained support.
#'
#' @return A non-negative scalar estimate of the Fisher divergence.
#' @export
fd_sensitivity <- function(fit, variables, score_candidate, lambda, score_ref = NULL) {
  draws <- as_reference_draws(fit, variables)

  ref <- if (is.null(score_ref)) {
    fd_reference_score(fit, draws, variables)
  } else {
    score_ref(draws)
  }
  candidate <- score_candidate(draws, lambda)

  ref <- validate_scores(ref, draws, "score_ref")
  candidate <- validate_scores(candidate, draws, "score_candidate")

  mean(rowSums((ref - candidate)^2))
}
