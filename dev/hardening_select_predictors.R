source("R/qa_fit.R")
source("R/qa_diagnose.R")

set.seed(1)
n <- 400
X1 <- rnorm(n); X2 <- rnorm(n)
y <- exp(1 + 0.5*X1 + rnorm(n))
z1 <- 0.4*X1 + rnorm(n)
donor  <- data.frame(y = y[1:200], X1 = X1[1:200], X2 = X2[1:200])
target <- data.frame(X1 = X1[201:400], X2 = X2[201:400], z1 = z1[201:400])

expect_error <- function(expr, label) {
  tryCatch({ force(expr); cat("FAILED (no error):", label, "\n") },
           error = function(e) cat("OK -", label, "->", conditionMessage(e), "\n"))
}
expect_warning <- function(expr, label) {
  tryCatch({ withCallingHandlers(force(expr),
    warning = function(w) { cat("OK -", label, "->", conditionMessage(w), "\n"); invokeRestart("muffleWarning") }) },
    error = function(e) cat("FAILED (error instead of warning):", label, "->", conditionMessage(e), "\n"))
}

cat("=== Type validation ===\n")
expect_error(qa_diagnose(as.matrix(donor), target, "y", "z1", c("X1","X2")), "donor_data not data.frame")
expect_error(qa_diagnose(donor, as.matrix(target), "y", "z1", c("X1","X2")), "target_data not data.frame")
expect_error(qa_diagnose(donor, target, 5, "z1", c("X1","X2")), "y_var not character")
expect_error(qa_diagnose(donor, target, "y", character(0), c("X1","X2")), "z_vars empty")
expect_error(qa_diagnose(donor, target, "y", "z1", character(0)), "x_vars empty")
expect_error(qa_diagnose(donor, target, "y", c("z1","z1"), c("X1","X2")), "duplicate z_vars")
expect_error(qa_diagnose(donor, target, "y", "z1", c("X1","X1")), "duplicate x_vars")

cat("\n=== Overlap validation ===\n")
expect_error(qa_diagnose(donor, target, "y", "X1", c("X1","X2")), "x_vars/z_vars overlap")
expect_error(qa_diagnose(donor, target, "y", "z1", c("y","X2")), "y_var in x_vars")

cat("\n=== Column existence / numeric ===\n")
expect_error(qa_diagnose(donor, target, "y", "zzz", c("X1","X2")), "z_vars missing from target_data")
expect_error(qa_diagnose(donor, target, "y", "z1", c("X1","XXX")), "x_vars missing from donor_data")
donor_char <- donor; donor_char$y <- as.character(donor_char$y)
expect_error(qa_diagnose(donor_char, target, "y", "z1", c("X1","X2")), "y_var not numeric")

cat("\n=== Missing values ===\n")
donor_na <- donor; donor_na$X1[1] <- NA
expect_error(qa_diagnose(donor_na, target, "y", "z1", c("X1","X2")), "donor_data has NA")
target_na <- target; target_na$z1[1] <- NA
expect_error(qa_diagnose(donor, target_na, "y", "z1", c("X1","X2")), "target_data has NA")

cat("\n=== Log-scale positivity ===\n")
donor_neg <- donor; donor_neg$y[1] <- -1
expect_error(qa_diagnose(donor_neg, target, "y", "z1", c("X1","X2"), outcome_scale = "log"), "y <= 0 under log scale")

cat("\n=== tau / B / k_folds / n_grid ===\n")
expect_error(qa_diagnose(donor, target, "y", "z1", c("X1","X2"), tau = -1), "negative tau")
expect_error(qa_diagnose(donor, target, "y", "z1", c("X1","X2"), B = 1), "B < 2")
expect_error(qa_diagnose(donor, target, "y", "z1", c("X1","X2"), k_folds = 1), "k_folds < 2")
expect_error(qa_diagnose(donor, target, "y", "z1", c("X1","X2"), n_grid = 2), "n_grid < 4")
expect_warning(qa_diagnose(donor, target, "y", "z1", c("X1","X2"), B = 5, verbose = FALSE), "B < 10 warning")

cat("\n=== Sanity: valid call still works ===\n")
res <- qa_diagnose(donor, target, "y", "z1", c("X1","X2"), B = 20, verbose = FALSE)
cat("OK - valid call: selected =", paste(res$selected_predictors, collapse=", "), "\n")
