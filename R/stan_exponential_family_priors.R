# Supported Stan prior families and their natural-parameter representations.

#' Supported Stan exponential-family priors
#'
#' List the prior statements currently recognised by the automatic quadratic
#' route.
#'
#' @return A data frame describing original parameters, natural parameters,
#'   and sufficient statistics.
#' @export
stan_exponential_family_priors <- function() {
  data.frame(
    stan_family = c("normal", "gamma", "beta", "inv_gamma"),
    original_parameters = c(
      "mu, sigma", "alpha, beta (shape, rate)", "alpha, beta",
      "alpha, beta (shape, scale)"
    ),
    natural_parameters = c(
      "eta1 = mu / sigma^2; eta2 = -1 / (2 sigma^2)",
      "eta1 = alpha - 1; eta2 = -beta",
      "eta1 = alpha - 1; eta2 = beta - 1",
      "eta1 = -(alpha + 1); eta2 = -beta"
    ),
    sufficient_statistics = c(
      "theta, theta^2", "log(theta), theta",
      "log(theta), log(1-theta)", "log(theta), 1/theta"
    ),
    stringsAsFactors = FALSE
  )
}

#' Inspect a Stan prior statement
#'
#' Conservatively detect a supported direct sampling statement such as
#' `theta ~ normal(prior_mean, prior_sd);` and resolve its reference natural
#' parameters.
#'
#' @param stan_file Path to a Stan program.
#' @param prior_variable Variable on the left-hand side of the prior statement.
#' @param stan_data Named list used to resolve scalar statement arguments.
#'
#' @return A list with `supported`, `family`, `natural_parameters`, and a
#'   diagnostic `reason`.
#' @export
stan_prior_info <- function(stan_file, prior_variable, stan_data = list()) {
  if (length(stan_file) != 1L || !file.exists(stan_file)) {
    stop("`stan_file` must identify an existing Stan program.", call. = FALSE)
  }
  if (length(prior_variable) != 1L ||
      !grepl("^[A-Za-z][A-Za-z0-9_]*(\\[[0-9]+\\])?$", prior_variable)) {
    stop(
      "`prior_variable` must be a Stan identifier or one fixed indexed entry.",
      call. = FALSE
    )
  }
  code <- paste(readLines(stan_file, warn = FALSE), collapse = "\n")
  code <- gsub("(?s)/\\*.*?\\*/", " ", code, perl = TRUE)
  code <- gsub("//[^\n]*", " ", code, perl = TRUE)
  escaped_variable <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", prior_variable)
  pattern <- paste0(
    "\\b", escaped_variable,
    "\\s*~\\s*([A-Za-z][A-Za-z0-9_]*)\\s*\\(([^;()]*)\\)\\s*;"
  )
  locations <- gregexpr(pattern, code, perl = TRUE)
  statements <- regmatches(code, locations)[[1L]]
  if (identical(statements, character(0L)) || identical(statements, "")) {
    return(list(
      supported = FALSE,
      variable = prior_variable,
      reason = paste0("No simple direct prior statement was found for `", prior_variable, "`.")
    ))
  }
  if (length(statements) != 1L) {
    return(list(
      supported = FALSE,
      variable = prior_variable,
      reason = paste0("Multiple direct prior statements were found for `", prior_variable, "`.")
    ))
  }

  groups <- regmatches(statements, regexec(pattern, statements, perl = TRUE))[[1L]]
  family <- groups[2L]
  registry <- stan_prior_registry()
  if (!family %in% names(registry)) {
    return(list(
      supported = FALSE,
      variable = prior_variable,
      family = family,
      reason = paste0("Stan prior family `", family, "` is not in the supported registry.")
    ))
  }

  arguments <- trimws(strsplit(groups[3L], ",", fixed = TRUE)[[1L]])
  if (length(arguments) != registry[[family]]$n_arguments) {
    return(list(
      supported = FALSE,
      variable = prior_variable,
      family = family,
      reason = "The prior arguments could not be matched to the supported form."
    ))
  }
  resolved <- lapply(arguments, resolve_stan_scalar, stan_data = stan_data)
  failures <- vapply(resolved, function(x) !isTRUE(x$ok), logical(1L))
  if (any(failures)) {
    return(list(
      supported = FALSE,
      variable = prior_variable,
      family = family,
      reason = paste0(
        "Could not resolve scalar Stan argument(s): ",
        paste(arguments[failures], collapse = ", "),
        ". Supply them by name in `stan_data`."
      )
    ))
  }
  original <- vapply(resolved, `[[`, numeric(1L), "value")
  natural <- registry[[family]]$to_natural(original)

  list(
    supported = TRUE,
    variable = prior_variable,
    family = family,
    statement = statements,
    original_arguments = stats::setNames(original, registry[[family]]$original_names),
    natural_parameters = natural,
    reason = "Supported exponential-family prior detected."
  )
}

stan_prior_registry <- function() {
  list(
    normal = list(
      n_arguments = 2L,
      original_names = c("mu", "sigma"),
      natural_names = c("eta1", "eta2"),
      to_natural = function(x) {
        if (x[2L] <= 0) stop("Normal `sigma` must be positive.", call. = FALSE)
        c(eta1 = x[1L] / x[2L]^2, eta2 = -1 / (2 * x[2L]^2))
      },
      gradient = function(theta) cbind(1, 2 * theta),
      support = function(theta) rep(TRUE, length(theta)),
      valid_box = function(lower, upper) upper["eta2"] < 0
    ),
    gamma = list(
      n_arguments = 2L,
      original_names = c("alpha", "beta"),
      natural_names = c("eta1", "eta2"),
      to_natural = function(x) {
        if (any(x <= 0)) stop("Gamma shape and rate must be positive.", call. = FALSE)
        c(eta1 = x[1L] - 1, eta2 = -x[2L])
      },
      gradient = function(theta) cbind(1 / theta, 1),
      support = function(theta) theta > 0,
      valid_box = function(lower, upper) lower["eta1"] > -1 && upper["eta2"] < 0
    ),
    beta = list(
      n_arguments = 2L,
      original_names = c("alpha", "beta"),
      natural_names = c("eta1", "eta2"),
      to_natural = function(x) {
        if (any(x <= 0)) stop("Beta shape parameters must be positive.", call. = FALSE)
        c(eta1 = x[1L] - 1, eta2 = x[2L] - 1)
      },
      gradient = function(theta) cbind(1 / theta, -1 / (1 - theta)),
      support = function(theta) theta > 0 & theta < 1,
      valid_box = function(lower, upper) all(lower > -1)
    ),
    inv_gamma = list(
      n_arguments = 2L,
      original_names = c("alpha", "beta"),
      natural_names = c("eta1", "eta2"),
      to_natural = function(x) {
        if (any(x <= 0)) stop("Inverse-Gamma shape and scale must be positive.", call. = FALSE)
        c(eta1 = -(x[1L] + 1), eta2 = -x[2L])
      },
      gradient = function(theta) cbind(1 / theta, -1 / theta^2),
      support = function(theta) theta > 0,
      valid_box = function(lower, upper) upper["eta1"] < -1 && upper["eta2"] < 0
    )
  )
}

resolve_stan_scalar <- function(expression, stan_data) {
  numeric_value <- suppressWarnings(as.numeric(expression))
  if (!is.na(numeric_value)) return(list(ok = TRUE, value = numeric_value))
  if (grepl("^[A-Za-z][A-Za-z0-9_]*$", expression) &&
      !is.null(stan_data[[expression]]) &&
      is.numeric(stan_data[[expression]]) && length(stan_data[[expression]]) == 1L &&
      is.finite(stan_data[[expression]])) {
    return(list(ok = TRUE, value = as.numeric(stan_data[[expression]])))
  }
  list(ok = FALSE, value = NA_real_)
}

validate_natural_bounds <- function(lower, upper, registry, family) {
  bounds <- validate_bounds(lower, upper)
  expected <- registry$natural_names
  if (length(bounds$lower) != length(expected)) {
    stop(
      "The `", family, "` quadratic route requires ", length(expected),
      " natural-parameter bounds: ", paste(expected, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (is.null(names(bounds$lower))) {
    names(bounds$lower) <- expected
    names(bounds$upper) <- expected
  } else {
    if (!setequal(names(bounds$lower), expected)) {
      stop(
        "Natural-parameter bounds must be named: ",
        paste(expected, collapse = ", "), ".",
        call. = FALSE
      )
    }
    bounds$lower <- bounds$lower[expected]
    bounds$upper <- bounds$upper[expected]
  }
  if (!registry$valid_box(bounds$lower, bounds$upper)) {
    stop(
      "The proposed box leaves the valid natural-parameter space for the `",
      family, "` prior.",
      call. = FALSE
    )
  }
  bounds
}

validate_detected_variable <- function(draws, prior_variable) {
  if (is.null(colnames(draws))) {
    stop("Reference draws must have parameter column names.", call. = FALSE)
  }
  if (grepl("\\[", prior_variable)) {
    belongs <- colnames(draws) == prior_variable
  } else {
    pattern <- paste0("^", prior_variable, "($|\\[)")
    belongs <- grepl(pattern, colnames(draws))
  }
  if (any(!belongs)) {
    stop(
      "All entries in `variables` must belong to detected prior variable `",
      prior_variable, "`.",
      call. = FALSE
    )
  }
}

sufficient_statistic_scores <- function(draws, registry, lambda) {
  if (any(!registry$support(as.numeric(draws)))) {
    stop("Reference draws fall outside the detected prior family's support.", call. = FALSE)
  }
  scores <- matrix(NA_real_, nrow = nrow(draws), ncol = ncol(draws))
  for (i in seq_len(nrow(draws))) {
    J <- registry$gradient(as.numeric(draws[i, ]))
    scores[i, ] <- as.numeric(J %*% lambda)
  }
  colnames(scores) <- colnames(draws)
  scores
}

sufficient_statistic_quadratic <- function(draws, registry, reference_scores) {
  if (any(!registry$support(as.numeric(draws)))) {
    stop("Reference draws fall outside the candidate prior family's support.", call. = FALSE)
  }
  p <- length(registry$natural_names)
  A <- matrix(0, nrow = p, ncol = p)
  b <- numeric(p)
  c_value <- 0
  for (i in seq_len(nrow(draws))) {
    J <- registry$gradient(as.numeric(draws[i, ]))
    ref <- as.numeric(reference_scores[i, ])
    A <- A + crossprod(J)
    b <- b - 2 * as.numeric(crossprod(J, ref))
    c_value <- c_value + sum(ref^2)
  }
  A <- A / nrow(draws)
  b <- b / nrow(draws)
  c_value <- c_value / nrow(draws)
  dimnames(A) <- list(registry$natural_names, registry$natural_names)
  names(b) <- registry$natural_names
  list(A = A, b = b, c = c_value)
}
