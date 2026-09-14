# Package functions (qa_diagnose, qa_fit) are
# available via devtools::load_all() / library(tstsqa) -- no manual
# source() needed once this is a real package.

# ---------------------------------------------------------------------------
# Shared toy dataset: X1, X3 genuinely predict y; X2 barely predicts y but
# strongly predicts z1 (the pathological predictor Section 5.4's S_i(z_k)
# screen is meant to catch); z2 depends only on the "good" predictors.
# Matches the dataset used throughout interactive validation.
# ---------------------------------------------------------------------------
make_toy_data <- function(seed = 7, n = 3000) {
  set.seed(seed)
  X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rnorm(n)
  y  <- exp(1 + 0.6 * X1 + 0.4 * X3 + rnorm(n, sd = 0.8))
  z1 <- 0.3 * X1 + 0.9 * X2 + rnorm(n, sd = 0.5)
  z2 <- 0.5 * X1 + 0.5 * X3 + rnorm(n, sd = 0.5)
  donor_idx  <- 1:(n / 2)
  target_idx <- (n / 2 + 1):n
  list(
    donor = data.frame(y = y[donor_idx], X1 = X1[donor_idx], X2 = X2[donor_idx], X3 = X3[donor_idx]),
    target = data.frame(X1 = X1[target_idx], X2 = X2[target_idx], X3 = X3[target_idx],
                         z1 = z1[target_idx], z2 = z2[target_idx])
  )
}

# ---------------------------------------------------------------------------
# Core behavior: the pathological predictor gets removed, real ones survive
# ---------------------------------------------------------------------------

test_that("the predictor that only explains z (not y) gets removed", {
  d <- make_toy_data()
  result <- suppressWarnings(qa_diagnose(
    d$donor, d$target, y_var = "y", z_vars = c("z1", "z2"), x_vars = c("X1", "X2", "X3"),
    outcome_scale = "log", tau = 10, B = 20, k_folds = 5, verbose = FALSE
  ))

  expect_false("X2" %in% result$selected_predictors)
  expect_true(all(c("X1", "X3") %in% result$selected_predictors))
  expect_true("X2" %in% vapply(result$removed_predictors, `[[`, character(1), "predictor"))
})

test_that("output structure is internally consistent", {
  d <- make_toy_data()
  result <- suppressWarnings(qa_diagnose(
    d$donor, d$target, y_var = "y", z_vars = c("z1", "z2"), x_vars = c("X1", "X2", "X3"),
    B = 20, k_folds = 5, verbose = FALSE
  ))

  expect_true(all(result$selected_predictors %in% result$stage1_survivors))
  expect_length(result$rho_star, 2)
  expect_named(result$rho_star, c("z1", "z2"))
  expect_equal(nrow(result$spec_search_results), 2^length(result$stage1_survivors) - 1)
  expect_true(is.numeric(result$R2_y_donor_oos_mean))
  expect_length(result$R2_y_donor_oos_ci, 2)
})

test_that("verbose = TRUE prints the documented progress messages", {
  d <- make_toy_data()
  expect_output(
    suppressWarnings(qa_diagnose(
      d$donor, d$target, y_var = "y", z_vars = c("z1", "z2"), x_vars = c("X1", "X2", "X3"),
      B = 5, k_folds = 5, verbose = TRUE
    )),
    "Removing 'X2'"
  )
})

test_that("factor x_vars are accepted (not just numeric)", {
  d <- make_toy_data()
  d$donor$region  <- factor(sample(c("A", "B", "C"), nrow(d$donor), replace = TRUE))
  d$target$region <- factor(sample(c("A", "B", "C"), nrow(d$target), replace = TRUE))
  result <- suppressWarnings(qa_diagnose(
    d$donor, d$target, y_var = "y", z_vars = c("z1", "z2"),
    x_vars = c("X1", "X2", "X3", "region"), B = 5, verbose = FALSE
  ))
  expect_true(is.list(result))
})

test_that("z_vars must still be numeric (factor z is rejected with a clear error)", {
  d <- make_toy_data()
  d$target$z1 <- factor(d$target$z1 > 0)
  expect_error(
    qa_diagnose(d$donor, d$target, y_var = "y", z_vars = c("z1", "z2"),
                 x_vars = c("X1", "X2", "X3"), B = 5, verbose = FALSE),
    "must be numeric"
  )
})

# ---------------------------------------------------------------------------
# Input validation: types and structure
# ---------------------------------------------------------------------------

test_that("non-data-frame donor_data/target_data are rejected", {
  d <- make_toy_data()
  expect_error(qa_diagnose(as.matrix(d$donor), d$target, "y", "z1", c("X1", "X2")),
               "must be a data frame")
  expect_error(qa_diagnose(d$donor, as.matrix(d$target), "y", "z1", c("X1", "X2")),
               "must be a data frame")
})

test_that("y_var/z_vars/x_vars type and length checks fire", {
  d <- make_toy_data()
  expect_error(qa_diagnose(d$donor, d$target, 5, "z1", c("X1", "X2")),
               "single character string")
  expect_error(qa_diagnose(d$donor, d$target, "y", character(0), c("X1", "X2")),
               "at least one column name")
  expect_error(qa_diagnose(d$donor, d$target, "y", "z1", character(0)),
               "at least one column name")
})

test_that("duplicate names in z_vars or x_vars are rejected", {
  d <- make_toy_data()
  expect_error(qa_diagnose(d$donor, d$target, "y", c("z1", "z1"), c("X1", "X2")),
               "duplicate")
  expect_error(qa_diagnose(d$donor, d$target, "y", "z1", c("X1", "X1")),
               "duplicate")
})

test_that("overlap between x_vars/z_vars/y_var is rejected", {
  d <- make_toy_data()
  expect_error(qa_diagnose(d$donor, d$target, "y", "X1", c("X1", "X2")),
               "must not overlap")
  expect_error(qa_diagnose(d$donor, d$target, "y", "z1", c("y", "X2")),
               "must not also appear")
})

test_that("missing columns are reported clearly", {
  d <- make_toy_data()
  expect_error(qa_diagnose(d$donor, d$target, "y", "zzz", c("X1", "X2")),
               "target_data is missing column")
  expect_error(qa_diagnose(d$donor, d$target, "y", "z1", c("X1", "XXX")),
               "donor_data is missing column")
})

test_that("non-numeric y_var is rejected", {
  d <- make_toy_data()
  d$donor$y <- as.character(d$donor$y)
  expect_error(qa_diagnose(d$donor, d$target, "y", "z1", c("X1", "X2")),
               "must be numeric")
})

# ---------------------------------------------------------------------------
# Input validation: missing values / log-scale positivity
# ---------------------------------------------------------------------------

test_that("missing values in donor_data or target_data are rejected", {
  d <- make_toy_data()
  d$donor$X1[1] <- NA
  expect_error(qa_diagnose(d$donor, d$target, "y", "z1", c("X1", "X2")),
               "missing values")

  d2 <- make_toy_data()
  d2$target$z1[1] <- NA
  expect_error(qa_diagnose(d2$donor, d2$target, "y", "z1", c("X1", "X2")),
               "missing values")
})

test_that("non-positive y under log scale is rejected", {
  d <- make_toy_data()
  d$donor$y[1] <- -1
  expect_error(qa_diagnose(d$donor, d$target, "y", "z1", c("X1", "X2"),
                                          outcome_scale = "log"),
               "strictly positive")
})

# ---------------------------------------------------------------------------
# Input validation: tau / B / k_folds / n_grid
# ---------------------------------------------------------------------------

test_that("tau, B, k_folds, n_grid validation fires", {
  d <- make_toy_data()
  expect_error(qa_diagnose(d$donor, d$target, "y", "z1", c("X1", "X2"), tau = -1),
               "positive numeric")
  expect_error(qa_diagnose(d$donor, d$target, "y", "z1", c("X1", "X2"), B = 1),
               "at least 2")
  expect_error(qa_diagnose(d$donor, d$target, "y", "z1", c("X1", "X2"), k_folds = 1),
               "at least 2")
  expect_error(qa_diagnose(d$donor, d$target, "y", "z1", c("X1", "X2"), n_grid = 2),
               "at least 4")
})

test_that("small B triggers a warning but still runs", {
  d <- make_toy_data()
  suppress_common_support <- function(expr) {
    withCallingHandlers(expr, warning = function(w) {
      if (grepl("common support", conditionMessage(w))) invokeRestart("muffleWarning")
    })
  }
  expect_warning(
    suppress_common_support(
      qa_diagnose(d$donor, d$target, "y", "z1", c("X1", "X2"),
                                 B = 5, verbose = FALSE)
    ),
    "quite small"
  )
})

# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

test_that("plotting = TRUE returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  d <- make_toy_data()
  result <- suppressWarnings(qa_diagnose(
    d$donor, d$target, y_var = "y", z_vars = c("z1", "z2"), x_vars = c("X1", "X2", "X3"),
    B = 20, k_folds = 5, verbose = FALSE, plotting = TRUE
  ))
  expect_s3_class(result$S_plot, "ggplot")
})

test_that("plotting = FALSE returns NULL for S_plot", {
  d <- make_toy_data()
  result <- suppressWarnings(qa_diagnose(
    d$donor, d$target, y_var = "y", z_vars = c("z1", "z2"), x_vars = c("X1", "X2", "X3"),
    B = 20, k_folds = 5, verbose = FALSE
  ))
  expect_null(result$S_plot)
})
