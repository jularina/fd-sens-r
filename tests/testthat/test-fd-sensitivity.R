test_that("fd_sensitivity returns zero for identical scores with an explicit score_ref", {
  draws <- matrix(c(-1, 0, 1), ncol = 1, dimnames = list(NULL, "theta"))
  fit <- mock_cmdstan_fit(draws, grad_fun = function(u) -u)
  score <- function(x) -x

  expect_equal(
    fd_sensitivity(
      fit,
      variables = "theta",
      score_candidate = function(x, lambda) score(x),
      lambda = 1,
      score_ref = score
    ),
    0
  )
})

test_that("fd_sensitivity derives score_ref automatically from the fit's gradient", {
  draws <- matrix(c(-1, 0, 1), ncol = 1, dimnames = list(NULL, "theta"))
  fit <- mock_cmdstan_fit(draws, grad_fun = function(u) -u)

  expect_equal(
    fd_sensitivity(
      fit,
      variables = "theta",
      score_candidate = function(x, lambda) -x,
      lambda = 1
    ),
    0
  )
})

test_that("fd_sensitivity rejects an automatic score_ref for constrained variables", {
  draws <- matrix(c(-1, 0, 1), ncol = 1, dimnames = list(NULL, "theta"))
  # Simulate a bounded parameter: the "unconstrained" draws differ from draws.
  fit <- mock_cmdstan_fit(draws, grad_fun = function(u) -u)
  fit$unconstrain_draws <- function(draws, format = "draws_matrix", ...) draws * 2

  expect_error(
    fd_sensitivity(
      fit,
      variables = "theta",
      score_candidate = function(x, lambda) -x,
      lambda = 1
    ),
    "unconstrained support"
  )
})

test_that("fd_prior_global_sensitivity finds the gaussian-prior endpoints", {
  set.seed(1)
  theta <- rnorm(100)
  draws <- matrix(theta, ncol = 1, dimnames = list(NULL, "theta"))
  prior_sd <- 2
  fit <- mock_cmdstan_fit(draws, grad_fun = function(u) -u / prior_sd^2)

  score_prior_ref <- function(x) -x / prior_sd^2
  score_prior_candidate <- function(x, lambda) (lambda - x) / prior_sd^2

  result <- fd_prior_global_sensitivity(
    fit,
    variables = "theta",
    lower = -2,
    upper = 2,
    score_prior_ref = score_prior_ref,
    score_prior_candidate = score_prior_candidate
  )

  expect_s3_class(result, "fd_sensitivity_result")
  expect_equal(result$fd_min, 0, tolerance = 1e-8)
  expect_equal(result$fd_max, 0.25, tolerance = 1e-6)
  expect_equal(result$sensitivity, 0.25, tolerance = 1e-6)
  expect_equal(abs(result$lambda_max), 2, tolerance = 1e-4)
})

test_that("fd_prior_global_sensitivity supports a 2-D hyperparameter search", {
  set.seed(2)
  theta <- rnorm(100)
  draws <- matrix(theta, ncol = 1, dimnames = list(NULL, "theta"))
  fit <- mock_cmdstan_fit(draws, grad_fun = function(u) -u)

  score_prior_candidate <- function(x, lambda) (lambda[1] - x) / lambda[2]^2

  lower <- c(-2, 0.5)
  upper <- c(2, 2)

  result <- fd_prior_global_sensitivity(
    fit,
    variables = "theta",
    lower = lower,
    upper = upper,
    score_prior_candidate = score_prior_candidate,
    score_prior_ref = function(x) -x
  )

  expect_s3_class(result, "fd_sensitivity_result")
  expect_length(result$lambda_min, 2)
  expect_length(result$lambda_max, 2)
  expect_true(all(result$lambda_min >= lower - 1e-8 & result$lambda_min <= upper + 1e-8))
  expect_true(all(result$lambda_max >= lower - 1e-8 & result$lambda_max <= upper + 1e-8))

  # The reference score is exactly recovered at lambda = (0, 1), inside the box.
  expect_equal(result$fd_min, 0, tolerance = 1e-6)

  # The objective is convex in lambda, so its maximum over the box is
  # attained at one of the 4 corners; compute those independently to check.
  corners <- expand.grid(lambda1 = c(lower[1], upper[1]), lambda2 = c(lower[2], upper[2]))
  corner_fd <- apply(corners, 1L, function(lambda) {
    candidate <- (lambda[1] - theta) / lambda[2]^2
    mean((-theta - candidate)^2)
  })
  expect_equal(result$fd_max, max(corner_fd), tolerance = 1e-6)
})

test_that("fd_prior_global_sensitivity supports named hyperparameter bounds", {
  set.seed(2)
  theta <- rnorm(100)
  draws <- matrix(theta, ncol = 1, dimnames = list(NULL, "theta"))
  fit <- mock_cmdstan_fit(draws, grad_fun = function(u) -u)

  seen_names <- NULL
  score_prior_candidate <- function(x, lambda) {
    seen_names <<- names(lambda)
    (lambda["mean"] - x) / lambda["sd"]^2
  }

  result <- fd_prior_global_sensitivity(
    fit,
    variables = "theta",
    lower = c(mean = -2, sd = 0.5),
    upper = c(sd = 2, mean = 2), # deliberately different order than `lower`
    score_prior_candidate = score_prior_candidate,
    score_prior_ref = function(x) -x
  )

  expect_setequal(seen_names, c("mean", "sd"))
  expect_named(result$lambda_min, c("mean", "sd"))
  expect_named(result$lambda_max, c("mean", "sd"))
})

test_that("fd_prior_global_sensitivity rejects mismatched bound names", {
  fit <- mock_cmdstan_fit(
    matrix(0, ncol = 1, dimnames = list(NULL, "theta")),
    grad_fun = function(u) -u
  )

  expect_error(
    fd_prior_global_sensitivity(
      fit,
      variables = "theta",
      lower = c(mean = -2, sd = 0.5),
      upper = c(mean = 2, scale = 2),
      score_prior_ref = function(x) -x,
      score_prior_candidate = function(x, lambda) -x
    ),
    "same hyperparameter names"
  )
})

test_that("fd_lr_global_sensitivity uses the analytic learning-rate objective", {
  theta <- c(-1, 0, 2)
  draws <- matrix(theta, ncol = 1, dimnames = list(NULL, "theta"))
  fit <- mock_cmdstan_fit(draws, grad_fun = function(u) -u)
  score_loss <- function(x) 2 * x

  result <- fd_lr_global_sensitivity(
    fit,
    variables = "theta",
    lambda_ref = 1,
    lower = 0.5,
    upper = 2,
    score_loss = score_loss
  )

  energy <- mean((2 * theta)^2)
  expect_equal(result$fd_min, 0)
  expect_equal(result$fd_max, energy)
  expect_equal(result$sensitivity, energy)
  expect_equal(unname(result$lambda_min), 1)
  expect_equal(unname(result$lambda_max), 2)
  expect_identical(result$analysis, "learning_rate")
})

test_that("fd_lr_global_sensitivity handles a reference rate outside the interval", {
  draws <- matrix(c(-1, 1), ncol = 1, dimnames = list(NULL, "theta"))
  fit <- mock_cmdstan_fit(draws, grad_fun = function(u) -u)

  result <- fd_lr_global_sensitivity(
    fit,
    variables = "theta",
    lambda_ref = 3,
    lower = 0.5,
    upper = 2,
    score_loss = function(x) x
  )

  expect_equal(unname(result$lambda_min), 2)
  expect_equal(unname(result$lambda_max), 0.5)
  expect_equal(result$fd_min, 1)
  expect_equal(result$fd_max, 6.25)
  expect_equal(result$sensitivity, 5.25)
})
