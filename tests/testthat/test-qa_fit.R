# Package functions (qa_fit) are available via devtools::load_all()
# / library(tstsqa) -- no manual source() needed once this is a real package.

# ---------------------------------------------------------------------------
# Shared toy dataset, reused across tests. Log-linear DGP with substantial
# unexplained variance, matching the dataset used during interactive
# validation (see tests/toy_example.R history).
# ---------------------------------------------------------------------------
make_toy_data <- function(seed = 123, n = 5000) {
  set.seed(seed)
  X <- rnorm(n)
  eps <- rnorm(n, sd = 1.5)
  y <- exp(1 + 0.8 * X + eps)
  list(
    donor  = data.frame(y = y[1:(n / 2)], X = X[1:(n / 2)]),
    target = data.frame(X = X[(n / 2 + 1):n]),
    true_y_target = y[(n / 2 + 1):n]
  )
}

# ---------------------------------------------------------------------------
# Core behavior: distributional recovery
# ---------------------------------------------------------------------------

test_that("quantile adjustment recovers the full target distribution (log scale)", {
  d <- make_toy_data()
  result <- suppressWarnings(qa_fit(
    d$donor, d$target, y_var = "y", x_vars = "X", outcome_scale = "log", n_grid = 200
  ))

  # Baseline imputation (no adjustment) mechanically compresses variance
  expect_lt(var(exp(result$y_hat_target)), var(d$true_y_target) / 10)

  # The adjusted distribution should closely track the true target quantiles
  probs <- c(.1, .25, .5, .75, .9)
  q_adj  <- quantile(result$y_adjusted, probs = probs)
  q_true <- quantile(d$true_y_target, probs = probs)
  expect_equal(unname(q_adj), unname(q_true), tolerance = 0.15)

  # And a formal two-sample KS test should not reject at the 1% level
  ks <- suppressWarnings(ks.test(result$y_adjusted, d$true_y_target))
  expect_gt(ks$p.value, 0.01)
})

test_that("quantile adjustment works on the level scale via a log-link GLM", {
  d <- make_toy_data()
  result <- suppressWarnings(qa_fit(
    d$donor, d$target, y_var = "y", x_vars = "X", outcome_scale = "level", n_grid = 200
  ))
  expect_s3_class(result$model, "glm")
  expect_true(is.numeric(result$y_adjusted))
  expect_length(result$y_adjusted, nrow(d$target))
})

test_that("family defaults to gaussian(link='log'), matching pre-existing behavior", {
  d <- make_toy_data()
  result_default  <- suppressWarnings(qa_fit(
    d$donor, d$target, y_var = "y", x_vars = "X", outcome_scale = "level"
  ))
  result_explicit <- suppressWarnings(qa_fit(
    d$donor, d$target, y_var = "y", x_vars = "X", outcome_scale = "level",
    family = gaussian(link = "log")
  ))
  expect_equal(unname(result_default$y_adjusted), unname(result_explicit$y_adjusted))
  expect_identical(result_default$model$family$family, "gaussian")
  expect_identical(result_default$model$family$link, "log")
})

test_that("a non-default family runs on the level scale and changes the fit", {
  d <- make_toy_data()
  result_gamma <- suppressWarnings(qa_fit(
    d$donor, d$target, y_var = "y", x_vars = "X", outcome_scale = "level",
    family = Gamma(link = "log")
  ))
  expect_s3_class(result_gamma$model, "glm")
  expect_identical(result_gamma$model$family$family, "Gamma")
  expect_length(result_gamma$y_adjusted, nrow(d$target))
})

test_that("family must be a family object", {
  d <- make_toy_data()
  expect_error(
    qa_fit(d$donor, d$target, "y", "X", outcome_scale = "level", family = "gaussian"),
    "family object"
  )
})

test_that("a non-default family on the log scale warns and is ignored", {
  d <- make_toy_data()
  # This toy data also triggers an unrelated common-support warning; isolate
  # the one this test actually checks for, same approach as the existing
  # annotate_p extrapolation test above.
  suppress_common_support <- function(expr) {
    withCallingHandlers(expr, warning = function(w) {
      if (grepl("common support", conditionMessage(w))) invokeRestart("muffleWarning")
    })
  }
  expect_warning(
    result <- suppress_common_support(
      qa_fit(d$donor, d$target, "y", "X", outcome_scale = "log", family = Gamma(link = "log"))
    ),
    "ignored when outcome_scale"
  )
  # log scale always fits lm(), never a GLM, regardless of the family argument
  expect_s3_class(result$model, "lm")
})

test_that("print() shows the family only for the level-scale branch", {
  d <- make_toy_data()
  result_log   <- suppressWarnings(qa_fit(d$donor, d$target, "y", "X", outcome_scale = "log"))
  result_level <- suppressWarnings(qa_fit(d$donor, d$target, "y", "X", outcome_scale = "level",
                                           family = Gamma(link = "log")))
  expect_false(grepl("Family:", paste(capture.output(print(result_log)), collapse = "\n")))
  expect_output(print(result_level), "Family:\\s+Gamma\\(link = \"log\"\\)")
})

test_that("the default family never spuriously triggers the log-scale ignored warning", {
  # Regression test: identical() on two family objects is FALSE even when
  # both are stats::gaussian(link = "log"), since each call to gaussian()
  # creates fresh closures with distinct environments. An earlier version
  # of this check used identical(family, stats::gaussian(link="log")) and
  # fired the "ignored" warning on EVERY outcome_scale = "log" call,
  # regardless of whether family was actually non-default.
  d <- make_toy_data()
  expect_no_ignored_warning <- function(expr) {
    fired <- FALSE
    withCallingHandlers(expr, warning = function(w) {
      if (grepl("ignored when outcome_scale", conditionMessage(w))) fired <<- TRUE
      invokeRestart("muffleWarning")
    })
    expect_false(fired)
  }
  expect_no_ignored_warning(qa_fit(d$donor, d$target, "y", "X", outcome_scale = "log"))
  expect_no_ignored_warning(qa_fit(d$donor, d$target, "y", "X", outcome_scale = "log",
                                    family = gaussian(link = "log")))
})

test_that("weighted quantiles run without error and preserve output length", {
  d <- make_toy_data()
  w <- runif(nrow(d$donor), 0.5, 1.5)
  result <- suppressWarnings(qa_fit(
    d$donor, d$target, y_var = "y", x_vars = "X", outcome_scale = "log",
    donor_weights = w
  ))
  expect_length(result$y_adjusted, nrow(d$target))
})

# ---------------------------------------------------------------------------
# Input validation: types and structure
# ---------------------------------------------------------------------------

test_that("non-data-frame donor_data/target_data are rejected", {
  d <- make_toy_data()
  expect_error(qa_fit(as.matrix(d$donor), d$target, "y", "X"),
               "must be a data frame")
  expect_error(qa_fit(d$donor, as.matrix(d$target), "y", "X"),
               "must be a data frame")
})

test_that("y_var and x_vars type/length checks fire", {
  d <- make_toy_data()
  expect_error(qa_fit(d$donor, d$target, c("y", "y"), "X"),
               "single character string")
  expect_error(qa_fit(d$donor, d$target, 5, "X"),
               "single character string")
  expect_error(qa_fit(d$donor, d$target, "y", character(0)),
               "at least one column name")
})

test_that("missing columns are reported clearly", {
  d <- make_toy_data()
  expect_error(qa_fit(d$donor, d$target, "y", "Z"),
               "donor_data is missing column")
  target_no_x <- data.frame(W = d$target$X)
  expect_error(qa_fit(d$donor, target_no_x, "y", "X"),
               "target_data is missing column")
})

test_that("non-numeric y_var is rejected", {
  d <- make_toy_data()
  d$donor$y <- as.character(d$donor$y)
  expect_error(qa_fit(d$donor, d$target, "y", "X"), "must be numeric")
})

# ---------------------------------------------------------------------------
# Input validation: missing values
# ---------------------------------------------------------------------------

test_that("missing values in donor_data are rejected", {
  d <- make_toy_data()
  d$donor$X[5] <- NA
  expect_error(qa_fit(d$donor, d$target, "y", "X"), "missing values")
})

test_that("missing values in target_data are rejected", {
  d <- make_toy_data()
  d$target$X[3] <- NA
  expect_error(qa_fit(d$donor, d$target, "y", "X"), "missing values")
})

# ---------------------------------------------------------------------------
# Input validation: outcome scale / weights / n_grid / annotate_p
# ---------------------------------------------------------------------------

test_that("non-positive y under log scale is rejected", {
  d <- make_toy_data()
  d$donor$y[1] <- -5
  expect_error(qa_fit(d$donor, d$target, "y", "X", outcome_scale = "log"),
               "strictly positive")
})

test_that("donor_weights validation fires for length mismatch, negatives, all-zero", {
  d <- make_toy_data()
  n_d <- nrow(d$donor)
  expect_error(qa_fit(d$donor, d$target, "y", "X",
                                    donor_weights = rep(1, n_d - 1)),
               "must match")
  expect_error(qa_fit(d$donor, d$target, "y", "X",
                                    donor_weights = rep(-1, n_d)),
               "nonnegative")
  expect_error(qa_fit(d$donor, d$target, "y", "X",
                                    donor_weights = rep(0, n_d)),
               "all zero")
})

test_that("n_grid validation fires for too-small, non-integer, non-scalar", {
  d <- make_toy_data()
  expect_error(qa_fit(d$donor, d$target, "y", "X", n_grid = 3),
               "at least 4")
  expect_error(qa_fit(d$donor, d$target, "y", "X", n_grid = 10.5),
               "at least 4")
  expect_error(qa_fit(d$donor, d$target, "y", "X", n_grid = c(10, 20)),
               "at least 4")
})

test_that("annotate_p validation fires in plot.qa_fit()", {
  skip_if_not_installed("ggplot2")
  d <- make_toy_data()
  result <- suppressWarnings(qa_fit(d$donor, d$target, "y", "X"))
  expect_error(plot(result, annotate_p = 1.5), "strictly between 0 and 1")
  expect_error(plot(result, annotate_p = 0), "strictly between 0 and 1")
  expect_error(plot(result, annotate_p = c(0.1, 0.9)), "strictly between 0 and 1")
})

# ---------------------------------------------------------------------------
# Warnings (not errors): common support, annotate_p extrapolation
# ---------------------------------------------------------------------------

test_that("common support violations warn but still return a result", {
  d <- make_toy_data()
  target_extreme <- d$target
  target_extreme$X[1] <- 10  # far outside donor's X range
  expect_warning(
    result <- qa_fit(d$donor, target_extreme, "y", "X"),
    "common support"
  )
  expect_length(result$y_adjusted, nrow(target_extreme))
})

test_that("annotate_p outside the estimated grid range warns", {
  skip_if_not_installed("ggplot2")
  d <- make_toy_data()
  # Suppress the (expected, incidental) common-support warning that a very
  # coarse n_grid = 4 also triggers, so this test isolates only the
  # annotate_p extrapolation warning it's actually checking for.
  suppress_common_support <- function(expr) {
    withCallingHandlers(expr, warning = function(w) {
      if (grepl("common support", conditionMessage(w))) invokeRestart("muffleWarning")
    })
  }
  result <- suppress_common_support(qa_fit(d$donor, d$target, "y", "X", n_grid = 4))
  expect_warning(plot(result, annotate_p = 0.999), "extrapolation")
})
test_that("target level absent from donor is rejected upfront", {
  d <- make_toy_data()
  d$donor$region  <- factor(sample(c("A", "B"), nrow(d$donor), replace = TRUE))
  d$target$region <- factor(sample(c("A", "B", "C"), nrow(d$target), replace = TRUE))  # "C" unseen by donor
  expect_error(
    qa_fit(d$donor, d$target, "y", c("X", "region")),
    "not present in donor_data"
  )
})

test_that("donor level absent from target warns but does not stop", {
  d <- make_toy_data()
  d$donor$region  <- factor(sample(c("A", "B", "C"), nrow(d$donor), replace = TRUE))
  d$target$region <- factor(sample(c("A", "B"), nrow(d$target), replace = TRUE))  # donor has extra "C"
  suppress_common_support <- function(expr) {
    withCallingHandlers(expr, warning = function(w) {
      if (grepl("common support", conditionMessage(w))) invokeRestart("muffleWarning")
    })
  }
  expect_warning(
    result <- suppress_common_support(qa_fit(d$donor, d$target, "y", c("X", "region"))),
    "not present in target_data"
  )
  expect_true(is.list(result))
})

# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

test_that("plot() on a qa_fit result returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  d <- make_toy_data()
  result <- suppressWarnings(qa_fit(d$donor, d$target, "y", "X"))
  p <- plot(result, annotate_p = 0.5)
  expect_s3_class(p, "ggplot")
})

test_that("plot() works without annotate_p", {
  skip_if_not_installed("ggplot2")
  d <- make_toy_data()
  result <- suppressWarnings(qa_fit(d$donor, d$target, "y", "X"))
  p <- plot(result)
  expect_s3_class(p, "ggplot")
})

test_that("qa_fit result has class qa_fit and a working print method", {
  d <- make_toy_data()
  result <- suppressWarnings(qa_fit(d$donor, d$target, "y", "X"))
  expect_s3_class(result, "qa_fit")
  expect_output(print(result), "qa_fit result")
})
