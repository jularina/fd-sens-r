# Helpers for interpreting an `fd_sensitivity_result`: compare a reference
# fit against a candidate fit (typically the one at `lambda_max`) via
# posterior quantiles, kernel density estimates, and empirical CDFs.
#
# All functions take `fits`, a named list of `CmdStanMCMC` or `stanfit` fits,
# e.g. list(reference = ref_fit, candidate = cand_fit) (mixing both kinds is
# fine, since draws are read via `fdsens::as_reference_draws()`), and write
# one plot file (and, for `plot_quantiles()`, one CSV) per call into
# `output_dir`, which is created if it does not already exist.

stopifnot(requireNamespace("fdsens", quietly = TRUE))
stopifnot(requireNamespace("posterior", quietly = TRUE))
stopifnot(requireNamespace("ggplot2", quietly = TRUE))
stopifnot(requireNamespace("jsonlite", quietly = TRUE))

# Blue for the first fit (typically the reference), purple for the second
# (typically the candidate), recycled if more fits are supplied.
FIT_PALETTE <- c("#2563EB", "#7C3AED")

fit_colors <- function(fit_names) {
  stats::setNames(rep_len(FIT_PALETTE, length(fit_names)), fit_names)
}

# A sequential blue gradient, one shade per block, ordered light-to-dark.
component_colors <- function(block_names) {
  n <- length(block_names)
  colors <- if (n <= 1L) {
    "#1D4ED8"
  } else {
    grDevices::colorRampPalette(c("#93C5FD", "#1E3A8A"))(n)
  }
  stats::setNames(colors, block_names)
}

# "black" or "white" per colour, whichever contrasts more (WCAG relative
# luminance), so labels stay legible across the whole light-to-dark gradient.
readable_text_color <- function(colors) {
  relative_luminance <- function(hex) {
    channel <- grDevices::col2rgb(hex) / 255
    channel <- ifelse(channel <= 0.03928, channel / 12.92, ((channel + 0.055) / 1.055)^2.4)
    sum(c(0.2126, 0.7152, 0.0722) * channel)
  }
  contrast <- function(hex, against) {
    l1 <- relative_luminance(hex)
    l2 <- relative_luminance(against)
    (max(l1, l2) + 0.05) / (min(l1, l2) + 0.05)
  }
  stats::setNames(
    vapply(colors, function(hex) {
      if (contrast(hex, "#FFFFFF") >= contrast(hex, "#000000")) "white" else "black"
    }, character(1L)),
    names(colors)
  )
}

# Long data frame of draws: one row per (fit, variable, draw).
stack_fit_draws <- function(fits, variables) {
  if (is.null(names(fits)) || any(names(fits) == "")) {
    stop("`fits` must be a named list, e.g. list(reference = ..., candidate = ...).", call. = FALSE)
  }
  rows <- lapply(names(fits), function(fit_name) {
    draws <- fdsens::as_reference_draws(fits[[fit_name]], variables)
    long <- as.data.frame(draws)
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
#' @param fits Named list of `CmdStanMCMC` or `stanfit` fits.
#' @param variables Character vector of parameters to summarise.
#' @param probs Quantile probabilities.
#' @return A data frame with one row per fit x variable (mean, sd, quantiles).
summarise_quantiles <- function(fits, variables, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) {
  if (is.null(names(fits)) || any(names(fits) == "")) {
    stop("`fits` must be a named list, e.g. list(reference = ..., candidate = ...).", call. = FALSE)
  }
  rows <- lapply(names(fits), function(fit_name) {
    draws <- posterior::as_draws_matrix(fdsens::as_reference_draws(fits[[fit_name]], variables))
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
#' @param fits Named list of `CmdStanMCMC` or `stanfit` fits.
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
#' @param fits Named list of `CmdStanMCMC` or `stanfit` fits.
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

# Recursively (a) drop draw-level fields by name at every nesting level, so
# a decomposition's `block_results` doesn't smuggle its blocks' draws back
# in, and (b) convert named atomic vectors to named lists, so jsonlite emits
# a JSON object (preserving keys like "eta1") instead of a bare array.
sanitize_for_json <- function(x, drop) {
  if (is.data.frame(x)) return(x)
  if (is.list(x)) {
    if (!is.null(names(x))) x <- x[setdiff(names(x), drop)]
    return(lapply(x, sanitize_for_json, drop = drop))
  }
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
#' the file a compact summary; for an `fd_sensitivity_decomposition`, this
#' also strips them out of each block in `block_results`.
#'
#' @param result An `fd_sensitivity_result` (or `fd_sensitivity_decomposition`).
#' @param output_dir Directory to write the file into.
#' @param file_stem Base file name (without extension) for the output.
#' @return The list that was written, invisibly.
save_sensitivity_result <- function(result, output_dir, file_stem = "sensitivity_result") {
  stopifnot(inherits(result, "fd_sensitivity_result"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  drop <- c("draws", "corners", "corner_fd", "A", "b", "c", "detection")
  values <- sanitize_for_json(unclass(result), drop = drop)

  jsonlite::write_json(
    values, file.path(output_dir, paste0(file_stem, ".json")),
    auto_unbox = TRUE, digits = NA, pretty = TRUE
  )

  invisible(values)
}

#' 100%-stacked bar of each block's share of total sensitivity
#'
#' Writes `<output_dir>/<file_stem>.png`: a single horizontal bar spanning
#' 0-100%, divided into one segment per block of an
#' `fd_sensitivity_decomposition`'s `components$sensitivity_share`, each
#' segment labelled with its percentage.
#'
#' @param result An `fd_sensitivity_decomposition`, i.e. the result of
#'   `fd_prior_global_sensitivity(..., independent = TRUE)`.
#' @param output_dir Directory to write the plot into.
#' @param file_stem Base file name (without extension) for the output.
#' @return The `ggplot` object, invisibly.
plot_component_shares <- function(result, output_dir, file_stem = "component_shares") {
  stopifnot(inherits(result, "fd_sensitivity_decomposition"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  components <- result$components
  if (any(is.na(components$sensitivity_share))) {
    stop(
      "`sensitivity_share` is undefined (total sensitivity is zero), so ",
      "shares cannot be plotted.",
      call. = FALSE
    )
  }
  components$block <- factor(components$block, levels = components$block)
  components$percent <- 100 * components$sensitivity_share
  palette <- component_colors(levels(components$block))
  components$text_color <- readable_text_color(palette)[as.character(components$block)]

  plot <- ggplot2::ggplot(
    components,
    ggplot2::aes(x = percent, y = "sensitivity", fill = block)
  ) +
    ggplot2::geom_col(position = ggplot2::position_stack(reverse = TRUE), width = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f%%", percent), color = text_color),
      position = ggplot2::position_stack(vjust = 0.5, reverse = TRUE),
      size = 3.2, fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
    ggplot2::labs(x = "share of total sensitivity (%)", y = NULL, fill = "block") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
  ggplot2::ggsave(file.path(output_dir, paste0(file_stem, ".png")), plot, width = 7, height = 2.5, bg = "white")

  invisible(plot)
}

#' Empirical CDF comparison across fits
#'
#' Writes `<output_dir>/<file_stem>.png`, one ECDF panel per variable with one
#' curve per fit.
#'
#' @param fits Named list of `CmdStanMCMC` or `stanfit` fits.
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
