# Helpers for interpreting an `fd_sensitivity_result`: compare a reference
# fit against a candidate fit (typically the one at `lambda_max`) via
# posterior quantiles, kernel density estimates, and empirical CDFs.
#
# All functions take `fits`, a named list of `CmdStanMCMC` fits, e.g.
#   list(reference = ref_fit, candidate = cand_fit)
# and write one plot file (and, for `plot_quantiles()`, one CSV) per call
# into `output_dir`, which is created if it does not already exist.

stopifnot(requireNamespace("posterior", quietly = TRUE))
stopifnot(requireNamespace("ggplot2", quietly = TRUE))
stopifnot(requireNamespace("jsonlite", quietly = TRUE))

# Blue for the first fit (typically the reference), purple for the second
# (typically the candidate), recycled if more fits are supplied.
FIT_PALETTE <- c("#2563EB", "#7C3AED")

fit_colors <- function(fit_names) {
  stats::setNames(rep_len(FIT_PALETTE, length(fit_names)), fit_names)
}

# Long data frame of draws: one row per (fit, variable, draw).
stack_fit_draws <- function(fits, variables) {
  if (is.null(names(fits)) || any(names(fits) == "")) {
    stop("`fits` must be a named list, e.g. list(reference = ..., candidate = ...).", call. = FALSE)
  }
  rows <- lapply(names(fits), function(fit_name) {
    draws <- fits[[fit_name]]$draws(variables = variables, format = "draws_matrix")
    long <- as.data.frame(as.matrix(draws))
    long$fit <- fit_name
    long <- stats::reshape(
      long,
      varying = variables,
      v.names = "value",
      timevar = "variable",
      times = variables,
      direction = "long"
    )
    long[, c("fit", "variable", "value")]
  })
  do.call(rbind, rows)
}

#' Posterior quantile table across fits
#'
#' @param fits Named list of `CmdStanMCMC` fits.
#' @param variables Character vector of parameters to summarise.
#' @param probs Quantile probabilities.
#' @return A data frame with one row per fit x variable (mean, sd, quantiles).
summarise_quantiles <- function(fits, variables, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) {
  if (is.null(names(fits)) || any(names(fits) == "")) {
    stop("`fits` must be a named list, e.g. list(reference = ..., candidate = ...).", call. = FALSE)
  }
  rows <- lapply(names(fits), function(fit_name) {
    draws <- fits[[fit_name]]$draws(variables = variables, format = "draws_matrix")
    summary <- posterior::summarise_draws(
      draws, mean, sd, ~ stats::quantile(.x, probs = probs)
    )
    cbind(fit = fit_name, summary)
  })
  do.call(rbind, rows)
}

#' Table and interval plot of posterior quantiles across fits
#'
#' Writes `<output_dir>/<file_stem>.csv` (the quantile table) and
#' `<output_dir>/<file_stem>.png` (a median/90%-interval comparison plot,
#' one colour per fit).
#'
#' @param fits Named list of `CmdStanMCMC` fits.
#' @param variables Character vector of parameters to summarise.
#' @param output_dir Directory to write the table and plot into.
#' @param probs Quantile probabilities.
#' @param file_stem Base file name (without extension) for the outputs.
#' @return The quantile table (data frame), invisibly.
plot_quantiles <- function(
    fits, variables, output_dir,
    probs = c(0.05, 0.25, 0.5, 0.75, 0.95), file_stem = "quantiles") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  table <- summarise_quantiles(fits, variables, probs)
  utils::write.csv(
    table, file.path(output_dir, paste0(file_stem, ".csv")), row.names = FALSE
  )

  long <- stack_fit_draws(fits, variables)
  interval_data <- do.call(rbind, lapply(split(long, list(long$fit, long$variable)), function(d) {
    q <- stats::quantile(d$value, probs = c(0.05, 0.5, 0.95))
    data.frame(fit = d$fit[1L], variable = d$variable[1L], low = q[1L], mid = q[2L], high = q[3L])
  }))

  plot <- ggplot2::ggplot(
    interval_data,
    ggplot2::aes(x = mid, y = variable, xmin = low, xmax = high, color = fit)
  ) +
    ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = 0.4)) +
    ggplot2::labs(x = "value (median, 90% interval)", y = NULL, color = "fit") +
    ggplot2::scale_color_manual(values = fit_colors(names(fits))) +
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(output_dir, paste0(file_stem, ".png")), plot, width = 6, height = 4, bg = "white")

  invisible(table)
}

#' Kernel density estimate comparison across fits
#'
#' Writes `<output_dir>/<file_stem>.png`, one density panel per variable with
#' one curve per fit.
#'
#' @param fits Named list of `CmdStanMCMC` fits.
#' @param variables Character vector of parameters to plot.
#' @param output_dir Directory to write the plot into.
#' @param file_stem Base file name (without extension) for the output.
#' @return The `ggplot` object, invisibly.
plot_kde <- function(fits, variables, output_dir, file_stem = "kde") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  long <- stack_fit_draws(fits, variables)

  plot <- ggplot2::ggplot(long, ggplot2::aes(x = value, color = fit, fill = fit)) +
    ggplot2::geom_density(alpha = 0.2) +
    ggplot2::facet_wrap(~variable, scales = "free") +
    ggplot2::labs(x = NULL, y = "density") +
    ggplot2::scale_color_manual(values = fit_colors(names(fits))) +
    ggplot2::scale_fill_manual(values = fit_colors(names(fits))) +
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(output_dir, paste0(file_stem, ".png")), plot, width = 6, height = 4, bg = "white")

  invisible(plot)
}

# Recursively convert named atomic vectors to named lists, so jsonlite emits
# a JSON object (preserving keys like "lambda1") instead of a bare array.
preserve_names <- function(x) {
  if (is.list(x)) return(lapply(x, preserve_names))
  if (!is.null(names(x))) return(as.list(x))
  x
}

#' Save an `fd_sensitivity_result`'s values to a JSON dict
#'
#' Writes `<output_dir>/<file_stem>.json` holding the result's scalar and
#' vector fields (`sensitivity`, `lambda_min`, `lambda_max`, `fd_min`,
#' `fd_max`, `interval`, `analysis`, and any other non-draw-level fields such
#' as `optimization` or `prior_family`). Draw-level fields (`draws`,
#' `corners`, `corner_fd`, `A`, `b`, `c`, `detection`) are left out to keep
#' the file a compact summary.
#'
#' @param result An `fd_sensitivity_result` object.
#' @param output_dir Directory to write the file into.
#' @param file_stem Base file name (without extension) for the output.
#' @return The list that was written, invisibly.
save_sensitivity_result <- function(result, output_dir, file_stem = "sensitivity_result") {
  stopifnot(inherits(result, "fd_sensitivity_result"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  drop <- c("draws", "corners", "corner_fd", "A", "b", "c", "detection")
  values <- unclass(result)[setdiff(names(result), drop)]
  # Preserve names (e.g. lambda1, lambda2) as JSON object keys instead of
  # jsonlite's default of dropping them and emitting a bare array.
  values <- preserve_names(values)

  jsonlite::write_json(
    values, file.path(output_dir, paste0(file_stem, ".json")),
    auto_unbox = TRUE, digits = NA, pretty = TRUE
  )

  invisible(values)
}

#' Empirical CDF comparison across fits
#'
#' Writes `<output_dir>/<file_stem>.png`, one ECDF panel per variable with one
#' curve per fit.
#'
#' @param fits Named list of `CmdStanMCMC` fits.
#' @param variables Character vector of parameters to plot.
#' @param output_dir Directory to write the plot into.
#' @param file_stem Base file name (without extension) for the output.
#' @return The `ggplot` object, invisibly.
plot_ecdf <- function(fits, variables, output_dir, file_stem = "ecdf") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  long <- stack_fit_draws(fits, variables)

  plot <- ggplot2::ggplot(long, ggplot2::aes(x = value, color = fit)) +
    ggplot2::stat_ecdf() +
    ggplot2::facet_wrap(~variable, scales = "free_x") +
    ggplot2::labs(x = NULL, y = "empirical CDF") +
    ggplot2::scale_color_manual(values = fit_colors(names(fits))) +
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(output_dir, paste0(file_stem, ".png")), plot, width = 6, height = 4, bg = "white")

  invisible(plot)
}
