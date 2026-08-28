#' Extract reference-posterior draws
#'
#' Extract a numeric matrix of posterior draws, one draw per row, from a
#' `cmdstanr` or `rstan` MCMC fit.
#'
#' @param fit A `CmdStanMCMC` or `stanfit` object.
#' @param variables Character vector naming the model parameters to extract.
#'
#' @return A numeric matrix.
#' @export
as_reference_draws <- function(fit, variables) {
  if (is.null(variables) || length(variables) == 0L) {
    stop("`variables` must be supplied.", call. = FALSE)
  }

  if (inherits(fit, "CmdStanMCMC")) {
    draws <- fit$draws(variables = variables, format = "matrix")
  } else if (inherits(fit, "stanfit")) {
    if (!requireNamespace("rstan", quietly = TRUE)) {
      stop("Package `rstan` is required to read a `stanfit` object.", call. = FALSE)
    }
    all_draws <- as.matrix(fit)
    missing_variables <- setdiff(variables, colnames(all_draws))
    if (length(missing_variables) > 0L) {
      stop(
        "Variables not found in `fit`: ",
        paste(missing_variables, collapse = ", "),
        call. = FALSE
      )
    }
    draws <- all_draws[, variables, drop = FALSE]
  } else {
    stop("`fit` must be a `CmdStanMCMC` or `stanfit` object.", call. = FALSE)
  }
  storage.mode(draws) <- "double"

  if (nrow(draws) == 0L || ncol(draws) == 0L) {
    stop("No posterior draws were extracted.", call. = FALSE)
  }
  if (any(!is.finite(draws))) {
    stop("Reference draws must contain only finite values.", call. = FALSE)
  }

  draws
}
