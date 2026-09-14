# ============================================================================
# qa_diagnose(): Single consolidated function implementing
# Section 5.4 end-to-end:
#
#   VARIABLE SELECTION (Stage 1): for each candidate predictor set, starting
#   from the full model, draw B subsamples of size n* = min(n, subsample_cap)
#   from donor and target (independently; without replacement if n exceeds
#   the cap, with replacement -- classic bootstrap -- otherwise), compute
#   S_i(z_k) (in-sample R^2 within each subsample) on every subsample, and
#   average. If the bootstrap MEAN of any (i, z_k) pair exceeds tau, drop
#   the predictor with the highest mean S and repeat on a fresh set of B
#   subsamples for the reduced predictor set. Stop when no pair's mean S
#   exceeds tau.
#
#   MODEL SELECTION (Stage 2): among the Stage 1 survivors, exhaustively
#   search all non-empty subsets and pick the one maximising mean k-fold
#   cross-validated R^2_{y,d}, bootstrapped over the SAME B subsamples
#   (paired resampling; size n* = min(n_donor, subsample_cap) each), so a
#   95% CI is available for the chosen specification's OOS R^2 alongside
#   the point estimate.
#
#   Finally, fits the selected specification on the FULL data and reports
#   rho* (Eq. 12) via qa_fit().
#
# HARDENED -- current state:
#   - Full input validation: types, required columns, numeric checks, no
#     missing values, y_var/x_vars/z_vars non-overlapping and non-duplicated,
#     outcome positivity under log scale, tau/B/k_folds/n_grid sanity checks.
#   - delta_y_i == 0 (exact) is guarded to avoid NaN/undefined ratios; a
#     tiny-but-nonzero delta_y_i is left as-is (correct, if extreme, paper
#     behavior -- not a bug).
#   - A warning (not an error) fires if Stage 1 survivors exceed 12, since
#     Stage 2's exhaustive 2^m - 1 subset search becomes very expensive.
#   - A warning fires if B < 10, since bootstrap means/CIs get unreliable.
#
# Still-open limitations:
#   - Stage 2 always does EXHAUSTIVE subset search; no non-exhaustive
#     fallback (e.g. stepwise) for large survivor sets, only a warning.
#   - Stage 2's B subsamples are drawn once and reused across all candidate
#     specifications (paired resampling), not literally identical to
#     Stage 1's per-iteration independent subsamples.
#   - z_k ~ X regressions are always plain OLS (unweighted), regardless of
#     outcome_scale, since z is not transformed anywhere in the paper.
#   - Printing is unconditional on verbose = TRUE and not throttled for very
#     large B or very large numbers of Stage 2 candidate specifications.
# ============================================================================

# NOTE: this file depends on qa_fit() (defined in
# qa_fit.R). In the eventual package structure, all R/ files are
# loaded together automatically and no explicit source() is needed. Until
# then, callers of this file must source("R/qa_fit.R") first.

# --- Small internal helper: fit the first stage and compute in-sample R^2 -
.fit_first_stage_r2 <- function(data, y_var, x_vars, outcome_scale) {
  if (outcome_scale == "log") {
    y_model <- log(data[[y_var]])
    fit_data <- cbind(data, y_model)
    f <- stats::reformulate(x_vars, response = "y_model")
    mod <- stats::lm(f, data = fit_data)
    y_hat <- stats::predict(mod, newdata = data)
  } else {
    y_model <- data[[y_var]]
    fit_data <- cbind(data, y_model)
    f <- stats::reformulate(x_vars, response = "y_model")
    mod <- stats::glm(f, data = fit_data, family = stats::gaussian(link = "log"))
    y_hat <- stats::predict(mod, newdata = data, type = "response")
  }
  r2 <- stats::var(y_hat) / stats::var(y_model)
  list(model = mod, y_model = y_model, y_hat = y_hat, r2 = r2)
}

# --- Small internal helper: k-fold out-of-fold predictions -----------------
# Tolerant of a fold's predict() call failing -- most commonly because a
# categorical predictor has a level that, by chance, ended up entirely in
# the held-out fold and never appeared in that fold's training data (R has
# no coefficient for a level it never saw, so predict() errors with "has
# new levels"). Rather than letting one unlucky random split crash the
# whole procedure, the affected fold's predictions are left as NA (with a
# warning), and excluded from the R^2 calculation downstream.
.kfold_oof_predictions <- function(data, fit_fun, predict_fun, k_folds) {
  n <- nrow(data)
  folds <- sample(rep(seq_len(k_folds), length.out = n))
  oof <- rep(NA_real_, n)
  failed_folds <- integer(0)
  for (k in seq_len(k_folds)) {
    train <- data[folds != k, , drop = FALSE]
    test  <- data[folds == k, , drop = FALSE]
    mod_k <- fit_fun(train)
    pred_k <- tryCatch(predict_fun(mod_k, test), error = function(e) NULL)
    if (is.null(pred_k)) {
      failed_folds <- c(failed_folds, k)
    } else {
      oof[folds == k] <- pred_k
    }
  }
  if (length(failed_folds) > 0) {
    warning(sprintf(
      "%d of %d cross-validation fold(s) could not generate predictions ",
      length(failed_folds), k_folds
    ), "(most likely a categorical predictor with a level absent from that ",
    "fold's training data). Affected observations were excluded from this R^2 estimate.")
  }
  oof
}

# --- Small internal helper: OOS R^2 for a given y ~ x_vars specification --
.cv_r2_first_stage <- function(data, y_var, x_vars, outcome_scale, k_folds) {
  fit_fun <- function(train_data) {
    if (outcome_scale == "log") {
      y_model_train <- log(train_data[[y_var]])
      fit_data <- cbind(train_data, y_model_train)
      f <- stats::reformulate(x_vars, response = "y_model_train")
      stats::lm(f, data = fit_data)
    } else {
      y_model_train <- train_data[[y_var]]
      fit_data <- cbind(train_data, y_model_train)
      f <- stats::reformulate(x_vars, response = "y_model_train")
      stats::glm(f, data = fit_data, family = stats::gaussian(link = "log"))
    }
  }
  predict_fun <- function(mod, test_data) {
    if (outcome_scale == "log") {
      stats::predict(mod, newdata = test_data)
    } else {
      stats::predict(mod, newdata = test_data, type = "response")
    }
  }
  y_model <- if (outcome_scale == "log") log(data[[y_var]]) else data[[y_var]]
  oof <- .kfold_oof_predictions(data, fit_fun, predict_fun, k_folds)
  # Restrict both numerator and denominator to the same valid subset, so a
  # partially-failed fold (see .kfold_oof_predictions) doesn't compare
  # predicted variance on a subset against observed variance on everyone.
  valid <- !is.na(oof)
  if (sum(valid) < 2) return(NA_real_)  # not enough valid predictions to estimate R^2 at all
  stats::var(oof[valid]) / stats::var(y_model[valid])
}

# --- Small internal helper: single-pass S_i(z_k) matrix, in-sample --------
.compute_S_matrix <- function(donor_data, target_data, y_var, z_vars, active, outcome_scale) {

  fit_y_full <- .fit_first_stage_r2(donor_data, y_var, active, outcome_scale)
  R2_y_full  <- fit_y_full$r2

  R2_z_full <- sapply(z_vars, function(z_k) {
    f_z <- stats::reformulate(active, response = z_k)
    mod_z <- stats::lm(f_z, data = target_data)
    stats::var(stats::fitted(mod_z)) / stats::var(target_data[[z_k]])
  })
  names(R2_z_full) <- z_vars

  S <- matrix(NA_real_, nrow = length(active), ncol = length(z_vars),
              dimnames = list(active, z_vars))

  if (length(active) > 1) {
    for (i in active) {
      reduced <- setdiff(active, i)

      fit_y_reduced <- .fit_first_stage_r2(donor_data, y_var, reduced, outcome_scale)
      delta_y_i <- R2_y_full - fit_y_reduced$r2

      for (z_k in z_vars) {
        f_z_reduced <- stats::reformulate(reduced, response = z_k)
        mod_z_reduced <- stats::lm(f_z_reduced, data = target_data)
        R2_z_reduced <- stats::var(stats::fitted(mod_z_reduced)) / stats::var(target_data[[z_k]])
        delta_z_ik <- R2_z_full[z_k] - R2_z_reduced

        # Safeguard: delta_y_i == 0 exactly would otherwise give 0/0 = NaN or
        # x/0 = +-Inf. NaN would corrupt max()/which() downstream, and we
        # would rather have a deliberate, well-defined value than a silent
        # propagation of NaN through the whole screening loop. A genuinely
        # tiny-but-nonzero delta_y_i (the "huge S" case discussed and
        # validated earlier) is NOT altered here -- it is correct, if
        # extreme, paper-defined behavior, not a bug.
        if (delta_y_i == 0) {
          S[i, z_k] <- if (delta_z_ik == 0) 0 else sign(delta_z_ik) * Inf
        } else {
          S[i, z_k] <- delta_z_ik / delta_y_i
        }
      }
    }
  }

  list(S = S, R2_y_full = R2_y_full, R2_z_full = R2_z_full)
}

# --- Small internal helper: draw one bootstrap/subsample from a data frame -
# cap = Inf (or any value >= nrow(data)) means every draw uses the FULL
# sample, with replacement (classic bootstrap) -- see subsample_cap in
# qa_diagnose() for why a user might want this.
.draw_subsample <- function(data, cap) {
  n <- nrow(data)
  n_star <- min(n, cap)
  replace <- n <= cap
  idx <- sample(seq_len(n), n_star, replace = replace)
  data[idx, , drop = FALSE]
}

#' Predictor Selection and Regime Diagnostic for TSTS Quantile Adjustment
#'
#' Implements Section 5.4 of Torres-Lopez, "Bias Reduction in TSTSLS
#' Applications: A Quantile Adjustment Method," end to end: bootstrap-based
#' screening of candidate predictors against the external variable(s)
#' \eqn{z_k}, a cross-validated search over the surviving predictors for the
#' best first-stage specification, and the regime diagnostic \eqn{\rho^*}
#' (Eq. 12) evaluated on the selected specification.
#'
#' @details
#' The function proceeds in two stages, followed by a final diagnostic:
#'
#' \strong{Stage 1 (Variable Selection).} Starting from the full candidate
#' set \code{x_vars}, at each iteration \code{B} subsamples of size
#' \eqn{n^* = \min(n, \code{subsample\_cap})} are drawn independently from
#' \code{donor_data} and \code{target_data} (without replacement if the
#' original sample size exceeds \code{subsample_cap}, with replacement --
#' i.e. a classic bootstrap -- otherwise).
#' On every subsample, the leave-one-out ratio
#' \eqn{S_i(z_k) = \Delta R^2_{z_k,t} / \Delta R^2_{y,d}} is computed for
#' every predictor \eqn{i} still in the active set and every external
#' variable \eqn{z_k}, using ordinary in-sample \eqn{R^2 = \mathrm{Var}(\hat
#' y)/\mathrm{Var}(y)} (or the analogous ratio for \eqn{z_k}). If the
#' bootstrap \emph{mean} of any \eqn{(i, z_k)} pair exceeds \code{tau}, the
#' predictor with the largest mean \eqn{S} is dropped and the process
#' repeats on a fresh set of \code{B} subsamples for the reduced predictor
#' set. Screening stops when no pair's mean \eqn{S} exceeds \code{tau} (or
#' when only one predictor remains).
#'
#' \strong{Stage 2 (Model Selection).} Among the Stage 1 survivors, every
#' non-empty subset is treated as a candidate specification (\eqn{2^m - 1}
#' candidates for \eqn{m} survivors). \code{B} subsamples of
#' \code{donor_data} (same \eqn{n^*} and replacement rule as Stage 1) are
#' drawn once and reused across all candidates (paired resampling, so
#' variability reflects the choice of specification rather than independent
#' sampling noise). For each candidate, \code{k_folds}-fold cross-validated
#' \eqn{R^2_{y,d}} is computed on every subsample; the specification with the
#' highest mean cross-validated \eqn{R^2} is selected.
#'
#' \strong{Final diagnostic.} The selected specification is fit once on the
#' full \code{donor_data}/\code{target_data} (not a subsample), and
#' \code{\link{qa_fit}} is called to obtain the realised
#' quantile-gap adjustment \eqn{\tilde\eta_t}, from which \eqn{\rho^*} (the
#' residual correlation at which the adjustment would exactly recover the
#' omitted covariance) is computed for each \code{z_vars} entry.
#'
#' When \code{verbose = TRUE} (the default), progress and results are
#' printed at each step: active predictors and bootstrap progress per
#' iteration, a message whenever a predictor is dropped, the top 5
#' \eqn{S_i(z_k)} pairs (mean and 95\% CI) per iteration, a convergence
#' message, the number of Stage 2 candidates and their evaluation progress,
#' the top 5 specifications by mean cross-validated \eqn{R^2} (with 95\%
#' CI), the selected specification, and the final \eqn{\rho^*} values.
#'
#' @param donor_data Data frame containing the donor sample. Must include
#'   \code{y_var} and all \code{x_vars}, with no missing values in these
#'   columns.
#' @param target_data Data frame containing the target sample. Must include
#'   all \code{x_vars} and \code{z_vars}, with no missing values in these
#'   columns. For any \code{x_vars} column that is a factor or character,
#'   every level present in \code{target_data} must also appear in
#'   \code{donor_data} (checked upfront and rejected with \code{stop()}
#'   otherwise); a level present in \code{donor_data} but absent from
#'   \code{target_data} triggers a \code{warning()} instead.
#' @param y_var Character string naming the outcome column in
#'   \code{donor_data}. Must not also appear in \code{x_vars} or
#'   \code{z_vars}.
#' @param z_vars Character vector naming one or more external-variable
#'   columns in \code{target_data} against which predictors are screened.
#'   Must not overlap with \code{x_vars}.
#' @param x_vars Character vector naming the candidate predictor columns,
#'   present in both \code{donor_data} and \code{target_data}.
#' @param outcome_scale Either \code{"log"} or \code{"level"}; see
#'   \code{\link{qa_fit}} for the distinction. Governs both the
#'   donor-side \eqn{R^2_{y,d}} regressions (Stages 1 and 2) and the final
#'   \code{qa_fit} call.
#' @param tau Numeric threshold for Stage 1 screening: a predictor is
#'   dropped when the bootstrap mean of any \eqn{S_i(z_k)} exceeds
#'   \code{tau}. Default 10, per the paper's suggested default (\code{tau =
#'   5} is more conservative).
#' @param B Number of bootstrap subsamples drawn at each Stage 1 iteration
#'   and reused across all Stage 2 candidates. Default 100. Must be at
#'   least 2; values below 10 trigger a warning, since bootstrap means and
#'   95\% CIs become unreliable with very few replications.
#' @param k_folds Number of folds used for the cross-validated \eqn{R^2} in
#'   Stage 2. Default 5.
#' @param n_grid Passed through to the final \code{\link{qa_fit}}
#'   call; see its documentation.
#' @param subsample_cap Upper limit on the size of each bootstrap/subsample
#'   draw: \eqn{n^* = \min(n, \code{subsample\_cap})}, drawn without
#'   replacement if the original sample exceeds this cap, with replacement
#'   (a classic bootstrap) otherwise. Default 5000, matching the paper's
#'   suggested default. Set to \code{Inf} (or any value at or above your
#'   sample size) to always bootstrap from the FULL sample instead of a
#'   capped subsample. This matters in particular for categorical
#'   predictors with rare levels: capping at a small \eqn{n^*} increases
#'   the chance that a bootstrap draw or cross-validation fold ends up
#'   with zero observations of some level, which causes
#'   \code{predict()} to fail with an "has new levels" error when that
#'   level does appear elsewhere. Raising (or removing) the cap reduces,
#'   but does not entirely eliminate, this risk -- an extremely rare level
#'   can still be excluded from a fold by chance even at the full sample
#'   size.
#' @param verbose Logical, default \code{TRUE}. If \code{TRUE}, prints
#'   progress and results at each step (see Details). Set to \code{FALSE}
#'   for silent operation.
#' @param plotting Logical, default \code{FALSE}. If \code{TRUE}, builds a
#'   \code{ggplot2} horizontal bar chart of the final Stage 1 \eqn{S_i(z_k)}
#'   values (mean and 95\% CI error bars) for the surviving predictors,
#'   with reference lines at \eqn{\tau = 5} and \eqn{\tau = 10}, returned as
#'   \code{S_plot}. Requires the \code{ggplot2} package to be installed;
#'   listed under \code{Suggests} rather than \code{Imports} so it is not a
#'   hard dependency for users who never request a plot.
#'
#' @return A list with components:
#'   \item{selected_predictors}{Character vector of predictors in the final
#'     specification, after both Stage 1 screening and Stage 2 model
#'     selection.}
#'   \item{stage1_survivors}{Character vector of predictors that survived
#'     Stage 1 screening, before the Stage 2 subset search (may differ from
#'     \code{selected_predictors} if a strict subset of the survivors fits
#'     \code{y} better).}
#'   \item{removed_predictors}{A list of removal events from Stage 1, each
#'     with the iteration number, the removed predictor, the \code{z_var}
#'     that triggered removal, and the bootstrap mean/95\% CI of \eqn{S} at
#'     the time of removal.}
#'   \item{spec_search_results}{Data frame of every Stage 2 candidate
#'     specification with its mean cross-validated \eqn{R^2} and 95\% CI.}
#'   \item{R2_y_donor_oos_mean}{Mean cross-validated \eqn{R^2_{y,d}} of the
#'     selected specification.}
#'   \item{R2_y_donor_oos_ci}{Named vector \code{c(lower, upper)}: the 95\%
#'     CI of the selected specification's cross-validated \eqn{R^2_{y,d}}.}
#'   \item{rho_star}{Named vector, one entry per \code{z_vars}, of the
#'     regime diagnostic \eqn{\rho^*} for the selected specification.}
#'   \item{qa_fit}{The full return value of the
#'     \code{\link{qa_fit}} call on the selected specification.}
#'   \item{S_plot}{A \code{ggplot} object (see \code{plotting} above), or
#'     \code{NULL} if \code{plotting = FALSE}.}
#'
#' @examples
#' \dontrun{
#' set.seed(7)
#' n <- 3000
#' X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rnorm(n)
#' y  <- exp(1 + 0.6 * X1 + 0.4 * X3 + rnorm(n, sd = 0.8))
#' z1 <- 0.3 * X1 + 0.9 * X2 + rnorm(n, sd = 0.5)  # X2 barely predicts y
#' z2 <- 0.5 * X1 + 0.5 * X3 + rnorm(n, sd = 0.5)
#' donor  <- data.frame(y = y[1:1500], X1 = X1[1:1500], X2 = X2[1:1500], X3 = X3[1:1500])
#' target <- data.frame(X1 = X1[1501:3000], X2 = X2[1501:3000], X3 = X3[1501:3000],
#'                       z1 = z1[1501:3000], z2 = z2[1501:3000])
#' result <- qa_diagnose(donor, target, y_var = "y",
#'                                     z_vars = c("z1", "z2"),
#'                                     x_vars = c("X1", "X2", "X3"))
#' }
#'
#' @export
qa_diagnose <- function(donor_data, target_data, y_var, z_vars, x_vars,
                                      outcome_scale = c("log", "level"),
                                      tau = 10, B = 100, k_folds = 5, n_grid = 200,
                                      subsample_cap = 5000,
                                      verbose = TRUE, plotting = FALSE) {

  outcome_scale <- match.arg(outcome_scale)

  # --- Input validation: basic types --------------------------------------
  if (!is.data.frame(donor_data))  stop("donor_data must be a data frame.")
  if (!is.data.frame(target_data)) stop("target_data must be a data frame.")
  if (!is.character(y_var) || length(y_var) != 1) {
    stop("y_var must be a single character string naming a column of donor_data.")
  }
  if (!is.character(z_vars) || length(z_vars) < 1) {
    stop("z_vars must be a character vector of at least one column name.")
  }
  if (!is.character(x_vars) || length(x_vars) < 1) {
    stop("x_vars must be a character vector of at least one column name.")
  }
  if (anyDuplicated(z_vars)) stop("z_vars contains duplicate names.")
  if (anyDuplicated(x_vars)) stop("x_vars contains duplicate names.")
  overlap <- intersect(x_vars, z_vars)
  if (length(overlap) > 0) {
    stop(sprintf("x_vars and z_vars must not overlap; found in both: %s",
                 paste(overlap, collapse = ", ")))
  }
  if (y_var %in% x_vars || y_var %in% z_vars) {
    stop("y_var must not also appear in x_vars or z_vars.")
  }

  # --- Input validation: required columns actually present ----------------
  missing_in_donor <- setdiff(c(y_var, x_vars), names(donor_data))
  if (length(missing_in_donor) > 0) {
    stop(sprintf("donor_data is missing column(s): %s", paste(missing_in_donor, collapse = ", ")))
  }
  missing_in_target <- setdiff(c(x_vars, z_vars), names(target_data))
  if (length(missing_in_target) > 0) {
    stop(sprintf("target_data is missing column(s): %s", paste(missing_in_target, collapse = ", ")))
  }
  if (!is.numeric(donor_data[[y_var]])) stop(sprintf("donor_data[['%s']] must be numeric.", y_var))
  # x_vars are only ever used as predictors (right-hand side of a formula),
  # where lm()/glm() handle factors and character columns natively -- no
  # numeric requirement here. z_vars, in contrast, are used as the RESPONSE
  # in z_k ~ x_vars regressions and directly in var()/cov() for the rho*
  # calculation, both of which require a numeric variable.
  for (v in z_vars) {
    if (!is.numeric(target_data[[v]])) stop(sprintf("target_data[['%s']] must be numeric.", v))
  }

  # --- Input validation: missing values -----------------------------------
  donor_cols_needed <- c(y_var, x_vars)
  donor_na <- !stats::complete.cases(donor_data[, donor_cols_needed, drop = FALSE])
  if (any(donor_na)) {
    stop(sprintf("donor_data contains missing values in %d row(s) among columns: %s.",
                 sum(donor_na), paste(donor_cols_needed, collapse = ", ")))
  }
  target_cols_needed <- c(x_vars, z_vars)
  target_na <- !stats::complete.cases(target_data[, target_cols_needed, drop = FALSE])
  if (any(target_na)) {
    stop(sprintf("target_data contains missing values in %d row(s) among columns: %s.",
                 sum(target_na), paste(target_cols_needed, collapse = ", ")))
  }

  # --- Input validation: donor/target categorical level consistency -------
  # Checked once, upfront, on the ORIGINAL x_vars -- a level present in
  # target but absent from donor would otherwise surface as a cryptic
  # "has new levels" error deep inside the bootstrap/CV loops, potentially
  # after substantial computation.
  .check_factor_levels(donor_data, target_data, x_vars)

  # --- Input validation: outcome positivity under log scale ---------------
  if (outcome_scale == "log" && any(donor_data[[y_var]] <= 0)) {
    stop(sprintf("outcome_scale = 'log' requires strictly positive values of '%s' in donor_data, ",
                 y_var), "but non-positive values were found. Use outcome_scale = 'level' instead.")
  }

  # --- Input validation: tau, B, k_folds, n_grid ---------------------------
  if (!is.numeric(tau) || length(tau) != 1 || tau <= 0) {
    stop("tau must be a single positive numeric value.")
  }
  if (!is.numeric(B) || length(B) != 1 || B != round(B) || B < 2) {
    stop("B must be a single integer of at least 2.")
  }
  if (B < 10) {
    warning(sprintf("B = %d is quite small; bootstrap means and 95%% CIs may be unreliable. ",
                     B), "Consider B >= 100 for stable results.")
  }
  if (!is.numeric(k_folds) || length(k_folds) != 1 || k_folds != round(k_folds) || k_folds < 2) {
    stop("k_folds must be a single integer of at least 2.")
  }
  if (!is.numeric(n_grid) || length(n_grid) != 1 || n_grid != round(n_grid) || n_grid < 4) {
    stop("n_grid must be a single integer of at least 4.")
  }
  if (!is.numeric(subsample_cap) || length(subsample_cap) != 1 || subsample_cap < 1) {
    stop("subsample_cap must be a single positive number (use Inf for the full sample).")
  }

  # --- Input validation: plotting dependency ------------------------------
  if (isTRUE(plotting) && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plotting = TRUE requires the 'ggplot2' package. ",
         "Install it with install.packages('ggplot2'), or call with plotting = FALSE.")
  }

  active <- x_vars
  removal_log <- list()
  iteration <- 0

  # ==========================================================================
  # VARIABLE SELECTION (Stage 1)
  # ==========================================================================
  if (verbose) cat("========== Variable Selection ==========\n")

  repeat {
    iteration <- iteration + 1

    # (a) Variables
    if (verbose) {
      cat(sprintf("\n--- Iteration %d ---\n", iteration))
      cat(sprintf("Variables: %s\n", paste(active, collapse = ", ")))
    }

    # (b) Bootstrap progress
    S_boot <- array(NA_real_, dim = c(length(active), length(z_vars), B),
                     dimnames = list(active, z_vars, NULL))
    boot_start_time <- Sys.time()
    for (b in seq_len(B)) {
      if (verbose) {
        elapsed <- as.numeric(difftime(Sys.time(), boot_start_time, units = "secs"))
        eta <- if (b > 1) elapsed / (b - 1) * (B - (b - 1)) else NA
        .print_progress(sprintf("Bootstrap progress: %d/%d (%.0f%%) - ETA: %s",
                                 b, B, 100 * b / B, .format_duration(eta)))
      }
      donor_sub  <- .draw_subsample(donor_data, subsample_cap)
      target_sub <- .draw_subsample(target_data, subsample_cap)
      S_boot[, , b] <- .compute_S_matrix(donor_sub, target_sub, y_var, z_vars,
                                          active, outcome_scale)$S
    }
    if (verbose) cat("\n")

    S_mean     <- apply(S_boot, c(1, 2), mean, na.rm = TRUE)
    S_ci_lower <- apply(S_boot, c(1, 2), stats::quantile, probs = 0.025, na.rm = TRUE)
    S_ci_upper <- apply(S_boot, c(1, 2), stats::quantile, probs = 0.975, na.rm = TRUE)

    # (d) Table of S(z_k) -- top 5
    S_table <- data.frame(
      predictor = rep(active, times = length(z_vars)),
      z_var     = rep(z_vars, each  = length(active)),
      mean_S    = as.vector(S_mean),
      ci_lower  = as.vector(S_ci_lower),
      ci_upper  = as.vector(S_ci_upper)
    )
    S_table <- S_table[order(-S_table$mean_S), ]
    rownames(S_table) <- NULL

    if (verbose) {
      cat("Top", min(5, nrow(S_table)), "S_i(z_k) pairs (mean, 95% CI across", B, "subsamples):\n")
      print(utils::head(S_table, 5), row.names = FALSE)
    }

    # Cannot screen further if only one predictor remains
    if (length(active) == 1) {
      if (verbose) cat("\nOnly one predictor remains; stopping Stage 1.\n")
      break
    }

    # (c) If at least one passes tau
    if (max(S_mean) > tau) {
      max_idx <- which(S_mean == max(S_mean), arr.ind = TRUE)[1, ]
      i_star <- rownames(S_mean)[max_idx["row"]]
      z_star <- colnames(S_mean)[max_idx["col"]]
      max_S  <- S_mean[max_idx["row"], max_idx["col"]]

      if (verbose) {
        cat(sprintf("\nAt least one pair exceeds tau = %.1f: mean S(%s, %s) = %.2f\n",
                     tau, i_star, z_star, max_S))
        cat(sprintf("Removing '%s'.\n", i_star))
      }

      active <- setdiff(active, i_star)
      removal_log[[length(removal_log) + 1]] <- list(
        iteration = iteration, predictor = i_star, z_var = z_star,
        mean_S = max_S, ci_lower = S_ci_lower[max_idx["row"], max_idx["col"]],
        ci_upper = S_ci_upper[max_idx["row"], max_idx["col"]]
      )
    } else {
      # (e) Convergence
      if (verbose) cat(sprintf("\nNo pair exceeds tau = %.1f. Variable selection converged.\n", tau))
      break
    }
  }

  if (verbose) {
    cat(sprintf("\nStage 1 survivors: %s\n", paste(active, collapse = ", ")))
  }

  # --- Optional: bar plot of the FINAL Stage 1 S_i(z_k), mean + 95% CI ----
  # Uses S_mean/S_ci_lower/S_ci_upper exactly as computed on the last loop
  # iteration -- these already correspond to `active` (== stage1 survivors),
  # since the loop breaks without recomputing them further.
  S_plot <- NULL
  if (isTRUE(plotting)) {
    plot_df <- data.frame(
      predictor = rep(active, times = length(z_vars)),
      z_var     = rep(z_vars, each  = length(active)),
      mean_S    = as.vector(S_mean),
      ci_lower  = as.vector(S_ci_lower),
      ci_upper  = as.vector(S_ci_upper)
    )
    plot_df$label <- paste0(plot_df$z_var, "\n", plot_df$predictor)

    # Order bars from highest to lowest mean S (top to bottom after coord_flip)
    plot_df <- plot_df[order(-plot_df$mean_S), ]
    plot_df$label <- factor(plot_df$label, levels = rev(plot_df$label))

    S_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$label, y = .data$mean_S,
                                                     fill = .data$z_var)) +
      ggplot2::geom_col(width = 0.65) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
                              width = 0.2, color = "gray30", linewidth = 0.5) +
      ggplot2::coord_flip() +
      ggplot2::geom_hline(yintercept = c(5, 10), linetype = "dashed", color = "gray50") +
      ggplot2::annotate("text", x = 0.6, y = 5,  label = "tau = 5",  vjust = -0.5,
                        size = 3.3, color = "gray30") +
      ggplot2::annotate("text", x = 0.6, y = 10, label = "tau = 10", vjust = -0.5,
                        size = 3.3, color = "gray30") +
      ggplot2::labs(x = NULL, y = expression(S[i](z[k])), fill = NULL) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.12))) +
      ggplot2::scale_fill_brewer(palette = "Set2") +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     legend.position = "top",
                     plot.background = ggplot2::element_rect(fill = "white", color = NA),
                     panel.background = ggplot2::element_rect(fill = "white", color = NA))
  }

  # ==========================================================================
  # MODEL SELECTION (Stage 2)
  # ==========================================================================
  if (verbose) cat("\n========== Model Selection ==========\n")

  if (length(active) > 12) {
    warning(sprintf(
      "%d predictors survived Stage 1 screening, giving 2^%d - 1 = %d candidate ",
      length(active), length(active), 2^length(active) - 1
    ), "specifications for Stage 2's exhaustive search. This may be very slow. ",
    "Consider a stricter tau in Stage 1, or a non-exhaustive search strategy.")
  }

  candidate_subsets <- unlist(
    lapply(seq_along(active), function(k) utils::combn(active, k, simplify = FALSE)),
    recursive = FALSE
  )
  n_candidates <- length(candidate_subsets)

  # (a) Number of models to test
  if (verbose) {
    cat(sprintf("Testing %d candidate specification(s) via %d-fold CV, bootstrapped over %d subsamples.\n",
                n_candidates, k_folds, B))
  }

  # Draw the B subsamples ONCE and reuse them for every candidate specification
  # (paired resampling): this isolates variability due to the choice of
  # specification, rather than adding extra noise from independent draws
  # per candidate, and is also more efficient (B draws instead of B * n_candidates).
  donor_subs_spec <- lapply(seq_len(B), function(b) .draw_subsample(donor_data, subsample_cap))

  r2_boot <- matrix(NA_real_, nrow = B, ncol = n_candidates)
  total_units <- n_candidates * B
  spec_start_time <- Sys.time()

  for (idx in seq_len(n_candidates)) {
    for (b in seq_len(B)) {
      # (b) Progress
      current_unit <- (idx - 1) * B + b
      if (verbose) {
        elapsed <- as.numeric(difftime(Sys.time(), spec_start_time, units = "secs"))
        eta <- if (current_unit > 1) {
          elapsed / (current_unit - 1) * (total_units - (current_unit - 1))
        } else NA
        .print_progress(sprintf(
          "Model %d/%d - Bootstrap %d/%d (%.0f%% overall) - ETA: %s",
          idx, n_candidates, b, B, 100 * current_unit / total_units, .format_duration(eta)
        ))
      }
      r2_boot[b, idx] <- .cv_r2_first_stage(donor_subs_spec[[b]], y_var,
                                             candidate_subsets[[idx]],
                                             outcome_scale, k_folds)
    }
    gc(verbose = FALSE)  # modest housekeeping between models, not expected to be the main
                         # lever on speed -- the slowdown with larger subsets is primarily
                         # the genuine cost of wider design matrices (see qa_diagnose() Details)
  }
  if (verbose) cat("\n")

  r2_mean     <- colMeans(r2_boot, na.rm = TRUE)
  r2_ci_lower <- apply(r2_boot, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
  r2_ci_upper <- apply(r2_boot, 2, stats::quantile, probs = 0.975, na.rm = TRUE)

  spec_results <- data.frame(
    specification = vapply(candidate_subsets, paste, collapse = " + ", FUN.VALUE = character(1)),
    mean_r2_oos   = r2_mean,
    ci_lower      = r2_ci_lower,
    ci_upper      = r2_ci_upper
  )

  # (c) Table with top 5 OOS R^2
  spec_results_sorted <- spec_results[order(-spec_results$mean_r2_oos), ]
  rownames(spec_results_sorted) <- NULL
  if (verbose) {
    cat("Top", min(5, nrow(spec_results_sorted)), "specifications by mean OOS R^2 (95% CI):\n")
    print(utils::head(spec_results_sorted, 5), row.names = FALSE)
  }

  best_idx <- which.max(spec_results$mean_r2_oos)
  final_predictors <- candidate_subsets[[best_idx]]

  # (d) Selected Model
  if (verbose) {
    cat(sprintf("\nSelected model: %s (mean OOS R^2 = %.4f, 95%% CI [%.4f, %.4f])\n",
                paste(final_predictors, collapse = ", "), spec_results$mean_r2_oos[best_idx],
                spec_results$ci_lower[best_idx], spec_results$ci_upper[best_idx]))
  }

  # ==========================================================================
  # Final fit on the FULL data, using the Stage-2-selected specification
  # ==========================================================================
  fit_y_final <- .fit_first_stage_r2(donor_data, y_var, final_predictors, outcome_scale)
  eps_d <- fit_y_final$y_model - fit_y_final$y_hat

  z_perp_t <- sapply(z_vars, function(z_k) {
    f_z <- stats::reformulate(final_predictors, response = z_k)
    mod_z <- stats::lm(f_z, data = target_data)
    stats::residuals(mod_z)
  })
  colnames(z_perp_t) <- z_vars

  qa_result <- qa_fit(donor_data, target_data, y_var, final_predictors,
                                    outcome_scale = outcome_scale, n_grid = n_grid)
  eta_t <- qa_result$eta_target

  rho_star <- sapply(z_vars, function(z_k) {
    stats::cov(eta_t, target_data[[z_k]]) /
      (stats::sd(eps_d) * stats::sd(z_perp_t[, z_k]))
  })
  names(rho_star) <- z_vars

  # (e) rho*
  if (verbose) {
    cat("\nrho*:\n")
    print(round(rho_star, 4))
  }

  list(
    selected_predictors = final_predictors,   # after BOTH stages
    stage1_survivors    = active,             # survivors of S-screening, before spec search
    removed_predictors  = removal_log,
    spec_search_results = spec_results,
    R2_y_donor_oos_mean = spec_results$mean_r2_oos[best_idx],
    R2_y_donor_oos_ci   = c(lower = spec_results$ci_lower[best_idx],
                            upper = spec_results$ci_upper[best_idx]),
    rho_star            = rho_star,
    qa_fit = qa_result,
    S_plot              = S_plot
  )
}
