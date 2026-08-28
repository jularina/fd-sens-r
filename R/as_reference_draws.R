#' Extract reference-posterior draws
#'
#' Extract a numeric matrix of posterior draws, one draw per row, from a
#' `cmdstanr` MCMC fit.
#'
#' @param fit A `CmdStanMCMC` object.
#' @param variables Character vector naming the model parameters to extract.
#'
#' @return A numeric matrix.
#' @export
as_reference_draws <- function(fit, variables) {
  if (!inherits(fit, "CmdStanMCMC")) {
    stop("`fit` must be a `CmdStanMCMC` fit.", call. = FALSE)
  }
  if (is.null(variables) || length(variables) == 0L) {
    stop("`variables` must be supplied.", call. = FALSE)
  }

  draws <- fit$draws(variables = variables, format = "matrix")
  storage.mode(draws) <- "double"

  if (nrow(draws) == 0L || ncol(draws) == 0L) {
    stop("No posterior draws were extracted.", call. = FALSE)
  }
  if (any(!is.finite(draws))) {
    stop("Reference draws must contain only finite values.", call. = FALSE)
  }

  draws
}
