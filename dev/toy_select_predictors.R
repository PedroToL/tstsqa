# ============================================================================
# Toy example for qa_diagnose().
#
# Design: X1 and X3 genuinely predict y; X2 barely predicts y at all but
# strongly predicts z1. This is exactly the pathological configuration
# Section 5.4's S_i(z_k) screen is meant to catch -- X2 should get flagged
# and removed (for a reasonably strict tau), since it enlarges the predictor
# channel for z1 without contributing to first-stage fit for y.
#
# z2 depends only on X1 and X3 (the "good" predictors), so no predictor
# should be flagged with respect to z2 -- a check that the screen doesn't
# over-trigger when there's no real problem for that particular z_k.
# ============================================================================

source("R/qa_fit.R")
source("R/qa_diagnose.R")

set.seed(7)
n <- 3000

X1 <- rnorm(n)
X2 <- rnorm(n)
X3 <- rnorm(n)

# y depends on X1, X3 -- X2's true coefficient is ~0
eps_y <- rnorm(n, sd = 0.8)
y <- exp(1 + 0.6 * X1 + 0.4 * X3 + eps_y)

# z1 leans heavily on X2 (the problem predictor), lightly on X1
eps_z1 <- rnorm(n, sd = 0.5)
z1 <- 0.3 * X1 + 0.9 * X2 + eps_z1

# z2 depends only on the "good" predictors X1, X3 -- no reason to flag anything here
eps_z2 <- rnorm(n, sd = 0.5)
z2 <- 0.5 * X1 + 0.5 * X3 + eps_z2

donor_idx  <- 1:1500
target_idx <- 1501:3000

donor <- data.frame(y = y[donor_idx], X1 = X1[donor_idx], X2 = X2[donor_idx], X3 = X3[donor_idx])
target <- data.frame(X1 = X1[target_idx], X2 = X2[target_idx], X3 = X3[target_idx],
                      z1 = z1[target_idx], z2 = z2[target_idx])

result <- qa_diagnose(
  donor_data = donor, target_data = target,
  y_var = "y", z_vars = c("z1", "z2"), x_vars = c("X1", "X2", "X3"),
  outcome_scale = "log", tau = 10
)

cat("=== Selected predictors ===\n")
print(result$selected_predictors)

cat("\n=== Removal log (should show X2 flagged for z1) ===\n")
if (length(result$removed_predictors) == 0) {
  cat("(nothing removed)\n")
} else {
  for (entry in result$removed_predictors) {
    cat(sprintf("iteration %d: removed '%s' (flagged via z_var = '%s', S = %.2f)\n",
                entry$iteration, entry$predictor, entry$z_var, entry$S))
  }
}

cat("\n=== R^2_y (donor, final predictor set) ===\n")
cat("in-sample:", result$R2_y_donor_insample, " | out-of-sample (5-fold CV):", result$R2_y_donor_oos, "\n")

cat("\n=== R^2_z (target, final predictor set) ===\n")
print(result$R2_z_target)

cat("\n=== rho* (per z_var) ===\n")
print(result$rho_star)
