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
    outcome_scale = "log", tau = 10, B = 20, verbose = 0
  ))

  expect_false("X2" %in% result$selected_predictors)
  expect_true(all(c("X1", "X3") %in% result$selected_predictors))
  expect_true("X2" %in% vapply(result$removed_predictors, `[[`, character(1), "predictor"))
})

test_that("output structure is internally consistent", {
  d <- make_toy_data()
  result <- suppressWarnings(qa_diagnose(
    d$donor, d$target, y_var = "y", z_vars = c("z1", "z2"), x_vars = c("X1", "X2", "X3"),
    B = 20, verbose = 0
  ))

  expect_true(all(c("X1", "X3") %in% result$selected_predictors))
  expect_named(result$R2_y_donor, c("mean", "ci_lower", "ci_upper"))
  expect_true(is.numeric(result$R2_y_donor$mean))
  expect_named(result$rho_star, c("mean", "ci_lower", "ci_upper"))
  expect_named(result$rho_star$mean, c("z1", "z2"))
  expect_named(result$rho_star$ci_lower, c("z1", "z2"))
  expect_named(result$rho_star$ci_upper, c("z1", "z2"))
})

test_that("verbose = 2 prints the documented progress messages", {
  d <- make_toy_data()
  expect_output(
    suppressWarnings(qa_diagnose(
      d$donor, d$target, y_var = "y", z_vars = c("z1", "z2"), x_vars = c("X1", "X2", "X3"),
      B = 5, verbose = 2
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
    x_vars = c("X1", "X2", "X3", "region"), B = 5, verbose = 0
  ))
  expect_true(is.list(result))
})

test_that("z_vars must still be numeric (factor z is rejected with a clear error)", {
  d <- make_toy_data()
  d$target$z1 <- factor(d$target$z1 > 0)
  expect_error(
    qa_diagnose(d$donor, d$target, y_var = "y", z_vars = c("z1", "z2"),
                 x_vars = c("X1", "X2", "X3"), B = 5, verbose = 0),
    "must be numeric"
  )
})
test_that("outcome_scale is printed as the very first line", {
  d <- make_toy_data()
  out <- capture.output(
    suppressWarnings(qa_diagnose(
      d$donor, d$target, y_var = "y", z_vars = c("z1", "z2"), x_vars = c("X1", "X2", "X3"),
      B = 5, verbose = 2, outcome_scale = "log"
    ))
  )
  expect_equal(out[1], "Outcome scale: log")
})
test_that("S_table has every pair, but plot() shows no more than 5 bars", {
  skip_if_not_installed("ggplot2")
  set.seed(3)
  n <- 3000
  X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rnorm(n); X4 <- rnorm(n)
  y  <- exp(1 + 0.5 * X1 + 0.3 * X3 + 0.2 * X4 + rnorm(n, sd = 0.8))
  z1 <- 0.3 * X1 + rnorm(n, sd = 0.5)
  z2 <- 0.4 * X3 + rnorm(n, sd = 0.5)
  donor  <- data.frame(y = y[1:1500], X1 = X1[1:1500], X2 = X2[1:1500],
                        X3 = X3[1:1500], X4 = X4[1:1500])
  target <- data.frame(X1 = X1[1501:3000], X2 = X2[1501:3000], X3 = X3[1501:3000],
                        X4 = X4[1501:3000], z1 = z1[1501:3000], z2 = z2[1501:3000])
  result <- suppressWarnings(qa_diagnose(donor, target, y_var = "y", z_vars = c("z1", "z2"),
                                           x_vars = c("X1", "X2", "X3", "X4"),
                                           B = 10, verbose = 0))
  # Full table has every (predictor, z_var) pair among survivors
  expect_equal(nrow(result$S_table), length(result$selected_predictors) * 2)
  # plot() defaults to showing only the top 5
  p <- plot(result)
  expect_lte(nrow(p$data), 5)
  # explicit n overrides the default
  p3 <- plot(result, n = 3)
  expect_lte(nrow(p3$data), 3)
})

test_that("donor_weights runs end-to-end and changes the result", {
  d <- make_toy_data()
  set.seed(99)
  w <- runif(nrow(d$donor), 0.5, 2)
  result_unweighted <- suppressWarnings(qa_diagnose(
    d$donor, d$target, "y", c("z1", "z2"), c("X1", "X2", "X3"), B = 10, verbose = 0
  ))
  result_weighted <- suppressWarnings(qa_diagnose(
    d$donor, d$target, "y", c("z1", "z2"), c("X1", "X2", "X3"), B = 10, verbose = 0,
    donor_weights = w
  ))
  expect_true(is.list(result_weighted))
  expect_false(isTRUE(all.equal(result_unweighted$R2_y_donor$mean,
                                  result_weighted$R2_y_donor$mean)))
})

test_that("donor_weights validation fires for length mismatch, negatives, all-zero", {
  d <- make_toy_data()
  n_d <- nrow(d$donor)
  expect_error(
    qa_diagnose(d$donor, d$target, "y", "z1", c("X1", "X2"), donor_weights = rep(1, n_d - 1)),
    "must match"
  )
  expect_error(
    qa_diagnose(d$donor, d$target, "y", "z1", c("X1", "X2"), donor_weights = rep(-1, n_d)),
    "nonnegative"
  )
  expect_error(
    qa_diagnose(d$donor, d$target, "y", "z1", c("X1", "X2"), donor_weights = rep(0, n_d)),
    "all zero"
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
test_that("subsample_cap controls the bootstrap draw size", {
  d <- make_toy_data()
  # Default (5000) should behave exactly as before for a small toy sample
  # (n=5000, so n* = min(5000,5000) = 5000, i.e. no change vs prior tests).
  # Just confirm the parameter is accepted and validated.
  expect_error(
    qa_diagnose(d$donor, d$target, "y", "z1", c("X1", "X2"),
                 subsample_cap = -1, verbose = 0),
    "subsample_cap must be"
  )
  result <- suppressWarnings(qa_diagnose(d$donor, d$target, "y", "z1", c("X1", "X2"),
                                           B = 5, verbose = 0, subsample_cap = Inf))
  expect_true(is.list(result))
})

# ---------------------------------------------------------------------------
# Input validation: tau / B / n_grid
# ---------------------------------------------------------------------------

test_that("tau, B, n_grid validation fires", {
  d <- make_toy_data()
  expect_error(qa_diagnose(d$donor, d$target, "y", "z1", c("X1", "X2"), tau = -1),
               "positive numeric")
  expect_error(qa_diagnose(d$donor, d$target, "y", "z1", c("X1", "X2"), B = 1),
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
                                 B = 5, verbose = 0)
    ),
    "quite small"
  )
})
test_that("target level absent from donor is rejected upfront (before any bootstrap)", {
  d <- make_toy_data()
  d$donor$region  <- factor(sample(c("A", "B"), nrow(d$donor), replace = TRUE))
  d$target$region <- factor(sample(c("A", "B", "C"), nrow(d$target), replace = TRUE))  # "C" unseen by donor
  expect_error(
    qa_diagnose(d$donor, d$target, "y", "z1", c("X1", "region"), verbose = 0),
    "not present in donor_data"
  )
})
test_that("a rare factor level unlucky in one bootstrap subsample does not crash the run", {
  set.seed(99)
  n <- 3000
  X1 <- rnorm(n)
  # Rare level: only ~4 rows total, occasionally isolated entirely into
  # donor or target by chance during a bootstrap draw in Final Estimation.
  region <- factor(sample(c("A", "B", "RARE"), n, replace = TRUE, prob = c(0.60, 0.399, 0.001)))
  y <- exp(1 + 0.5 * X1 + rnorm(n, sd = 0.8))
  z1 <- 0.4 * X1 + rnorm(n, sd = 0.5)
  donor  <- data.frame(y = y[1:1500], X1 = X1[1:1500], region = region[1:1500])
  target <- data.frame(X1 = X1[1501:3000], region = region[1501:3000], z1 = z1[1501:3000])

  result <- tryCatch({
    suppressWarnings(qa_diagnose(donor, target, y_var = "y", z_vars = "z1",
                                   x_vars = c("X1", "region"), B = 10, verbose = 0,
                                   subsample_cap = Inf))
    "ok"
  }, error = function(e) conditionMessage(e))
  expect_equal(result, "ok")
})

# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

test_that("plot() on a qa_diagnose result returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  d <- make_toy_data()
  result <- suppressWarnings(qa_diagnose(
    d$donor, d$target, y_var = "y", z_vars = c("z1", "z2"), x_vars = c("X1", "X2", "X3"),
    B = 20, verbose = 0
  ))
  expect_s3_class(plot(result), "ggplot")
})

test_that("plot() has no error bars", {
  skip_if_not_installed("ggplot2")
  d <- make_toy_data()
  result <- suppressWarnings(qa_diagnose(
    d$donor, d$target, y_var = "y", z_vars = c("z1", "z2"), x_vars = c("X1", "X2", "X3"),
    B = 20, verbose = 0
  ))
  p <- plot(result)
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_false("GeomErrorbar" %in% layer_classes)
})

test_that("qa_diagnose result has class qa_diagnose and a working print method", {
  d <- make_toy_data()
  result <- suppressWarnings(qa_diagnose(
    d$donor, d$target, y_var = "y", z_vars = c("z1", "z2"), x_vars = c("X1", "X2", "X3"),
    B = 20, verbose = 0
  ))
  expect_s3_class(result, "qa_diagnose")
  expect_output(print(result), "qa_diagnose result")
  expect_output(print(result), "rho")
})
