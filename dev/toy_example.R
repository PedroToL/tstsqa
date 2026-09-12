# ============================================================================
# Toy dataset to exercise qa_fit() end-to-end:
#   - one donor/target pair simulated so both a log-scale and a level-scale
#     first stage make sense
#   - one call with outcome_scale = "log"
#   - one call with outcome_scale = "level"
#   - both with plotting = TRUE and an annotated quantile position
# ============================================================================

source("R/qa_fit.R")

set.seed(42)

n <- 4000

# Two predictors, so x_vars has more than one entry (exercises reformulate()
# with a vector rather than a single string).
X1 <- rnorm(n)
X2 <- rnorm(n, mean = 0, sd = 1)

# True log-linear DGP: log(y) = 1 + 0.6*X1 - 0.3*X2 + eps
eps <- rnorm(n, sd = 1.2)
y <- exp(1 + 0.6 * X1 - 0.3 * X2 + eps)

donor_idx  <- 1:2000
target_idx <- 2001:4000

donor  <- data.frame(y = y[donor_idx],  X1 = X1[donor_idx],  X2 = X2[donor_idx])
target <- data.frame(                   X1 = X1[target_idx], X2 = X2[target_idx])

true_y_target <- y[target_idx]

# ---------------------------------------------------------------------------
# Model 1: log scale (matches the true DGP; lm(log(y) ~ X1 + X2))
# ---------------------------------------------------------------------------
cat("=== Log-scale model ===\n")
result_log <- qa_fit(
  donor_data = donor, target_data = target,
  y_var = "y", x_vars = c("X1", "X2"),
  outcome_scale = "log", n_grid = 200,
  plotting = TRUE, annotate_p = 0.9
)

cat("class(model):", paste(class(result_log$model), collapse = "/"), "\n")
cat("var(y_adjusted):", var(result_log$y_adjusted),
    " | var(true_y_target):", var(true_y_target), "\n")

qprobs <- c(.1, .25, .5, .75, .9)
cat("\nQuantiles -- adjusted vs. true target:\n")
print(round(rbind(
  adjusted = quantile(result_log$y_adjusted, probs = qprobs),
  true     = quantile(true_y_target,         probs = qprobs)
), 3))

cat("\nKS test (adjusted vs. true target):\n")
print(ks.test(result_log$y_adjusted, true_y_target))

cat("\nFor comparison, KS test (raw exponentiated prediction vs. true target),\n")
cat("i.e. baseline imputation with no adjustment at all:\n")
print(ks.test(exp(result_log$y_hat_target), true_y_target))

cat("class(eta_plot):", paste(class(result_log$eta_plot), collapse = "/"), "\n")

ggplot2::ggsave("/home/claude/eta_plot_log.png", result_log$eta_plot,
                 width = 6, height = 4, dpi = 120)

# Sanity check: annotated eta(0.9) should match eta_spline's own prediction
eta_from_spline <- predict(result_log$eta_spline, 0.9)$y
eta_from_grid   <- approx(result_log$p_grid, result_log$Q_y_d, xout = 0.9)$y -
                    approx(result_log$p_grid, result_log$Q_yhat_d, xout = 0.9)$y
cat("eta(0.9) via spline:", eta_from_spline, " | via raw grid quantiles:", eta_from_grid, "\n")

# ---------------------------------------------------------------------------
# Model 2: level scale (glm gaussian log-link, applied to the same data,
# even though the DGP is log-linear -- this just exercises the code path;
# fit quality is not expected to be as good as the correctly-specified log model)
# ---------------------------------------------------------------------------
cat("\n=== Level-scale model ===\n")
result_level <- qa_fit(
  donor_data = donor, target_data = target,
  y_var = "y", x_vars = c("X1", "X2"),
  outcome_scale = "level", n_grid = 200,
  plotting = TRUE, annotate_p = 0.5
)

cat("class(model):", paste(class(result_level$model), collapse = "/"), "\n")
cat("var(y_adjusted):", var(result_level$y_adjusted),
    " | var(true_y_target):", var(true_y_target), "\n")

cat("\nQuantiles -- adjusted vs. true target:\n")
print(round(rbind(
  adjusted = quantile(result_level$y_adjusted, probs = qprobs),
  true     = quantile(true_y_target,           probs = qprobs)
), 3))

cat("\nKS test (adjusted vs. true target):\n")
print(ks.test(result_level$y_adjusted, true_y_target))

cat("\nFor comparison, KS test (raw prediction vs. true target),\n")
cat("i.e. baseline imputation with no adjustment at all:\n")
print(ks.test(result_level$y_hat_target, true_y_target))

ggplot2::ggsave("/home/claude/eta_plot_level.png", result_level$eta_plot,
                 width = 6, height = 4, dpi = 120)

cat("\nBoth runs completed. Plots saved to eta_plot_log.png and eta_plot_level.png\n")
