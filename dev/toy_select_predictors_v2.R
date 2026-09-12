# ============================================================================
# Toy example for the rewritten qa_diagnose():
#   Stage 1 (bootstrap S(z_k) screening) should still remove X2 (barely
#   predicts y, strongly predicts z1).
#   Stage 2 (subset search by OOS R^2) should then pick the best subset of
#   the Stage-1 survivors (X1, X3) -- with this DGP, that should be both of
#   them, since both genuinely help predict y.
# ============================================================================

source("R/qa_fit.R")
source("R/qa_diagnose.R")

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

result <- qa_diagnose(
  donor_data = donor, target_data = target,
  y_var = "y", z_vars = c("z1", "z2"), x_vars = c("X1", "X2", "X3"),
  outcome_scale = "log", tau = 10, B = 20, k_folds = 5, plotting = TRUE
)

cat("\n\n================ SUMMARY ================\n")
cat("Stage 1 survivors:", paste(result$stage1_survivors, collapse = ", "), "\n")
cat("Final selected predictors (post Stage 2):", paste(result$selected_predictors, collapse = ", "), "\n")

cat("\nRemoval log:\n")
for (entry in result$removed_predictors) {
  cat(sprintf("  iter %d: removed '%s' (mean S(%s) = %.2f, 95%% CI [%.2f, %.2f])\n",
              entry$iteration, entry$predictor, entry$z_var, entry$mean_S, entry$ci_lower, entry$ci_upper))
}

cat("\nStage 2 specification search results:\n")
print(result$spec_search_results)

cat("\nBest OOS R^2 (mean, 95% CI):", result$R2_y_donor_oos_mean,
    "[", result$R2_y_donor_oos_ci["lower"], ",", result$R2_y_donor_oos_ci["upper"], "]\n")
cat("\nrho*:\n")
print(result$rho_star)

if (!is.null(result$S_plot)) {
  ggplot2::ggsave("/home/claude/S_barplot_final.png", result$S_plot, width = 6, height = 4, dpi = 120)
  cat("\nSaved final S(z_k) bar plot to S_barplot_final.png\n")
}
