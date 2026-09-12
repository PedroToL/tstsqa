# ============================================================================
# Toy example for bootstrap_regime_diagnostics(), reusing the same DGP as
# tests/toy_select_predictors.R: X2 barely predicts y but strongly predicts
# z1, so it should be removed in the full-sample selection; z2 depends only
# on the "good" predictors X1, X3.
# ============================================================================

source("R/bootstrap_regime_diagnostics.R")

set.seed(7)
n <- 3000

X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rnorm(n)
eps_y  <- rnorm(n, sd = 0.8)
y  <- exp(1 + 0.6 * X1 + 0.4 * X3 + eps_y)
eps_z1 <- rnorm(n, sd = 0.5)
z1 <- 0.3 * X1 + 0.9 * X2 + eps_z1
eps_z2 <- rnorm(n, sd = 0.5)
z2 <- 0.5 * X1 + 0.5 * X3 + eps_z2

donor_idx  <- 1:1500
target_idx <- 1501:3000
donor  <- data.frame(y = y[donor_idx], X1 = X1[donor_idx], X2 = X2[donor_idx], X3 = X3[donor_idx])
target <- data.frame(X1 = X1[target_idx], X2 = X2[target_idx], X3 = X3[target_idx],
                      z1 = z1[target_idx], z2 = z2[target_idx])

boot_result <- bootstrap_regime_diagnostics(
  donor_data = donor, target_data = target,
  y_var = "y", z_vars = c("z1", "z2"), x_vars = c("X1", "X2", "X3"),
  outcome_scale = "log", tau = 10, B = 30, plotting = TRUE
)

cat("=== n* and sampling scheme ===\n")
cat("n_star_donor:", boot_result$n_star_donor, " replace:", boot_result$donor_replace, "\n")
cat("n_star_target:", boot_result$n_star_target, " replace:", boot_result$target_replace, "\n")

cat("\n=== Full-sample selected predictors (fixed across all subsamples) ===\n")
print(boot_result$full_result$selected_predictors)

cat("\n=== S_i(z_k): mean across", dim(boot_result$S_boot)[3], "subsamples ===\n")
print(round(boot_result$S_summary$mean, 2))

cat("\n=== S_i(z_k): IQR across subsamples ===\n")
print(round(boot_result$S_summary$iqr, 2))

cat("\n=== rho*: mean across subsamples ===\n")
print(round(boot_result$rho_summary$mean, 4))

cat("\n=== rho*: IQR across subsamples ===\n")
print(round(boot_result$rho_summary$iqr, 4))

cat("\n=== Full-sample point estimates (for comparison) ===\n")
cat("rho* (full sample):\n")
print(round(boot_result$full_result$rho_star, 4))

cat("\n=== top5 table (already printed during run via verbose=TRUE) ===\n")
print(boot_result$top5)

if (!is.null(boot_result$S_plot)) {
  ggplot2::ggsave("/home/claude/S_barplot.png", boot_result$S_plot, width = 6, height = 4, dpi = 120)
  cat("\nSaved S(z_k) bar plot to S_barplot.png\n")
}
