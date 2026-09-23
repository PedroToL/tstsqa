# ============================================================================
# qa_diagnose(): Variable selection and regime diagnostic implementing
# Section 5.4's screening procedure, plus the regime diagnostic (Eq. 12):
#
#   VARIABLE SELECTION: for each candidate predictor set, starting from the
#   full model, draw B subsamples of size n* = min(n, subsample_cap) from
#   donor and target (independently; without replacement if n exceeds the
#   cap, with replacement -- classic bootstrap -- otherwise), compute
#   S_i(z_k) (in-sample R^2 within each subsample) on every subsample, and
#   average. If the bootstrap MEAN of any (i, z_k) pair exceeds tau, drop
#   the predictor with the highest mean S and repeat on a fresh set of B
#   subsamples for the reduced predictor set. Stop when no pair's mean S
#   exceeds tau. There is no further subset/specification search: all
#   survivors of this screening are the final predictor set.
#
#   FINAL ESTIMATION: using the fixed set of survivors, draw a fresh B
#   subsamples and, on each, compute in-sample R^2_{y,d} and the regime
#   diagnostic rho* (Eq. 12), giving a bootstrap mean and 95% CI for both
#   rather than a single point estimate. A final qa_fit() call on the FULL
#   data (not a subsample) is also returned, for the actual adjustment
#   applied to the whole target sample.
#
#   Optional donor_weights/target_weights (survey/design weights) are used
#   EVERYWHERE the corresponding sample is involved: the lm()/glm()
#   first-stage fits, R^2_{y,d}, S_i(z_k), the z_k ~ x_vars regressions in
#   the target sample, the quantile step inside every qa_fit() call, and the
#   rho*/theta/channel calculations. Weighted regression is the point: OLS
#   forces Cov(eps, X) = 0 only in the metric the model was fit in, so an
#   unweighted fit paired with weighted moments leaves Cov_weighted(eps, X)
#   nonzero and breaks the lambda = rho*/rho identity. Uniform weights
#   reproduce the unweighted results exactly.
#
# HARDENED -- current state:
#   - Full input validation: types, required columns, numeric checks, no
#     missing values, y_var/x_vars/z_vars non-overlapping and non-duplicated,
#     outcome positivity under log scale, tau/B/n_grid sanity checks.
#   - delta_y_i == 0 (exact) is guarded to avoid NaN/undefined ratios; a
#     tiny-but-nonzero delta_y_i is left as-is (correct, if extreme, paper
#     behavior -- not a bug).
#   - A warning fires if B < 10, since bootstrap means/CIs get unreliable.
#   - Each Final Estimation bootstrap draw's qa_fit() call is tolerant of a
#     donor/target level mismatch specific to that subsample (rare, but
#     possible even when the full donor/target agree overall): that
#     replicate's rho* is skipped (left NA) rather than crashing the run.
#
# Still-open limitations:
#   - z_k ~ X regressions are always plain OLS (unweighted), regardless of
#     outcome_scale, since z is not transformed anywhere in the paper.
#   - Printing is not throttled for very large B, at whatever verbose level
#     is active.
# ============================================================================

# NOTE: this file depends on qa_fit() (defined in
# qa_fit.R). In the eventual package structure, all R/ files are
# loaded together automatically and no explicit source() is needed. Until
# then, callers of this file must source("R/qa_fit.R") first.

# --- Small internal helper: fit the first stage and compute in-sample R^2 -
# The fit IS weighted when weights are supplied, and the R^2 is then the
# weighted one. An earlier version left both unweighted on the grounds that
# the first stage is a prediction device, but that leaves Cov_weighted(eps,
# X) nonzero (OLS only zeroes it in the metric it was fit in), which breaks
# the lambda = rho*/rho identity the diagnostics rely on. Fitting weighted
# keeps the fit, its R^2, and the population moments below all in one
# metric. Uniform weights reproduce the unweighted result exactly.
.fit_first_stage_r2 <- function(data, y_var, x_vars, outcome_scale, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, nrow(data))
  if (outcome_scale == "log") {
    y_model <- log(data[[y_var]])
    fit_data <- cbind(data, y_model)
    f <- stats::reformulate(x_vars, response = "y_model")
    mod <- stats::lm(f, data = fit_data, weights = weights)
    y_hat <- stats::predict(mod, newdata = data)
  } else {
    y_model <- data[[y_var]]
    fit_data <- cbind(data, y_model)
    f <- stats::reformulate(x_vars, response = "y_model")
    mod <- stats::glm(f, data = fit_data, family = stats::gaussian(link = "log"),
                      weights = weights)
    y_hat <- stats::predict(mod, newdata = data, type = "response")
  }
  # Weighted squared correlation between the outcome and its fitted values.
  # With uniform weights this is exactly cor(y_model, y_hat)^2, which in
  # turn equals summary(mod)$r.squared for the OLS branch with an intercept;
  # for the GLM branch, which has no r.squared slot, it is a bounded
  # pseudo-R^2. Using the weighted version keeps R^2 consistent with the
  # weighted fit that produced y_hat.
  r2 <- .wcov(y_model, y_hat, weights)^2 /
    (.wvar(y_model, weights) * .wvar(y_hat, weights))
  list(model = mod, y_model = y_model, y_hat = y_hat, r2 = r2)
}



# --- Small internal helper: single-pass S_i(z_k) matrix, in-sample --------
# donor_weights weights the donor-side first stage, target_weights the
# target-side z_k ~ x_vars regressions, so both sides of S_i(z_k) are
# estimated in the same metric as the population moments elsewhere.
.compute_S_matrix <- function(donor_data, target_data, y_var, z_vars, active, outcome_scale,
                               donor_weights = NULL, target_weights = NULL) {

  if (is.null(target_weights)) target_weights <- rep(1, nrow(target_data))

  fit_y_full <- .fit_first_stage_r2(donor_data, y_var, active, outcome_scale, donor_weights)
  R2_y_full  <- fit_y_full$r2

  R2_z_full <- sapply(z_vars, function(z_k) {
    f_z <- stats::reformulate(active, response = z_k)
    mod_z <- stats::lm(f_z, data = target_data, weights = target_weights)
    .wvar(stats::fitted(mod_z), target_weights) /
      .wvar(target_data[[z_k]], target_weights)
  })
  names(R2_z_full) <- z_vars

  S <- matrix(NA_real_, nrow = length(active), ncol = length(z_vars),
              dimnames = list(active, z_vars))

  if (length(active) > 1) {
    for (i in active) {
      reduced <- setdiff(active, i)

      fit_y_reduced <- .fit_first_stage_r2(donor_data, y_var, reduced, outcome_scale, donor_weights)
      delta_y_i <- R2_y_full - fit_y_reduced$r2

      for (z_k in z_vars) {
        f_z_reduced <- stats::reformulate(reduced, response = z_k)
        mod_z_reduced <- stats::lm(f_z_reduced, data = target_data, weights = target_weights)
        R2_z_reduced <- .wvar(stats::fitted(mod_z_reduced), target_weights) /
          .wvar(target_data[[z_k]], target_weights)
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

# --- Small internal helper: subsample a data frame AND an accompanying
# weights vector using the SAME drawn indices, so weights stay aligned with
# whichever rows were actually sampled. Used only for donor_data when
# donor_weights is supplied; target_data and the no-weights case use the
# plain .draw_subsample() above.
.draw_subsample_with_weights <- function(data, weights, cap) {
  n <- nrow(data)
  n_star <- min(n, cap)
  replace <- n <= cap
  idx <- sample(seq_len(n), n_star, replace = replace)
  list(data = data[idx, , drop = FALSE], weights = weights[idx])
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
#' The function proceeds in two phases:
#'
#' \strong{Variable Selection.} Starting from the full candidate
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
#' when only one predictor remains). There is no further specification
#' search: all survivors form the final predictor set directly.
#'
#' \strong{Final Estimation.} Using the fixed set of survivors, a fresh
#' \code{B} subsamples are drawn (same rule as above) and, on each, both
#' in-sample \eqn{R^2_{y,d}} and the regime diagnostic \eqn{\rho^*} (Eq. 12)
#' are computed -- the latter via a \code{\link{qa_fit}} call on that
#' subsample, from which the realised quantile-gap adjustment
#' \eqn{\tilde\eta_t} is obtained. Both are reported as a bootstrap mean and
#' 95\% CI rather than a single point estimate. Finally, \code{\link{qa_fit}}
#' is called once more on the FULL data (not a subsample), giving the actual
#' adjustment applied to the whole target sample.
#'
#' \code{verbose} controls how much of this is printed as it happens (see
#' the \code{verbose} parameter below); \code{1} shows section headers and
#' bootstrap progress with percentage and ETA, and \code{2} additionally
#' shows active predictors per iteration, a message whenever a predictor is
#' dropped, the top 5 \eqn{S_i(z_k)} pairs (mean and 95\% CI) per iteration,
#' a convergence message, and the final \eqn{R^2_{y,d}} and \eqn{\rho^*}
#' (mean and 95\% CI). Regardless of \code{verbose}, the full results are
#' always available afterward via \code{print()} on the returned object.
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
#'   donor-side \eqn{R^2_{y,d}} regressions and every \code{qa_fit} call.
#' @param tau Numeric threshold for screening: a predictor is
#'   dropped when the bootstrap mean of any \eqn{S_i(z_k)} exceeds
#'   \code{tau}. Default 10, per the paper's suggested default (\code{tau =
#'   5} is more conservative).
#' @param B Number of bootstrap subsamples drawn at each Variable Selection
#'   iteration, and again (fresh) for Final Estimation. Default 100. Must
#'   be at least 2; values below 10 trigger a warning, since bootstrap
#'   means and 95\% CIs become unreliable with very few replications.
#' @param n_grid Passed through to every \code{\link{qa_fit}}
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
#' @param donor_weights Optional numeric vector, length \code{nrow(donor_data)},
#'   of nonnegative survey/design weights for the donor sample. Weighting is
#'   applied neither to the \code{lm()}/\code{glm()} first-stage fits, nor to
#'   \eqn{R^2_{y,d}}, nor to \eqn{S_i(z_k)} -- all three are unweighted
#'   devices. \code{donor_weights} instead enters the quantile estimation
#'   inside every \code{qa_fit} call (including the final one), and the
#'   donor-side standard deviation in \eqn{\rho^*}. Whenever the donor
#'   sample is bootstrapped, the corresponding weights are resampled using
#'   the identical drawn indices, so they stay aligned with whichever rows
#'   were actually sampled. Defaults to \code{NULL} (uniform weights).
#' @param target_weights Optional numeric vector, length
#'   \code{nrow(target_data)}, of nonnegative survey/design weights for the
#'   target sample. Used only in the \eqn{\rho^*} calculation, where both
#'   \eqn{\mathrm{Cov}(\tilde\eta, z_k)} and \eqn{\sigma_{z_\perp}} are
#'   target-sample population moments. Supply this whenever
#'   \code{donor_weights} is supplied and any downstream \eqn{\lambda} is
#'   computed with weights, since \eqn{\lambda = \rho^*/\rho} holds only
#'   when both sides are computed in the same metric. Defaults to
#'   \code{NULL} (uniform weights).
#' @param verbose Integer: \code{0} for silent operation, \code{1} to print
#'   only progress (section headers and bootstrap progress/ETA), or
#'   \code{2} (default) for the full print described in Details (adds
#'   per-iteration variable lists, the top-5 \eqn{S_i(z_k)} table, removal/
#'   convergence messages, and the final \eqn{R^2_{y,d}}/\eqn{\rho^*} table).
#'
#' @return A list (class \code{"qa_diagnose"}) with components:
#'   \item{selected_predictors}{Character vector of predictors that
#'     survived Variable Selection; the final predictor set.}
#'   \item{removed_predictors}{A list of removal events, each with the
#'     iteration number, the removed predictor, the \code{z_var} that
#'     triggered removal, and the bootstrap mean/95\% CI of \eqn{S} at the
#'     time of removal.}
#'   \item{S_table}{Data frame of every (predictor, z_var) pair among the
#'     final survivors, with \code{mean_S}, \code{ci_lower}, and
#'     \code{ci_upper} (from the last Variable Selection iteration), ordered
#'     highest to lowest \code{mean_S}. \code{plot()} shows only the top few
#'     rows; use this directly for the complete table.}
#'   \item{R2_y_donor}{A list with \code{mean}, \code{ci_lower}, and
#'     \code{ci_upper}: the bootstrapped in-sample \eqn{R^2_{y,d}} of the
#'     final predictor set. Always unweighted, regardless of
#'     \code{donor_weights}; see Details.}
#'   \item{rho_star}{A list with \code{mean}, \code{ci_lower}, and
#'     \code{ci_upper}, each a named vector (one entry per \code{z_vars}):
#'     the bootstrapped regime diagnostic \eqn{\rho^*}. Weighted by
#'     \code{donor_weights}/\code{target_weights} whenever either is
#'     supplied.}
#'   \item{cov_yhat_eta}{A list with \code{mean}, \code{ci_lower}, and
#'     \code{ci_upper}: the bootstrapped target-sample covariance between
#'     the first-stage prediction and the realised quantile gap,
#'     \eqn{\mathrm{Cov}(\hat y, \tilde\eta)}, for the final predictor set.
#'     Computed on the model scale (log scale when
#'     \code{outcome_scale = "log"}), which is the scale on which
#'     \eqn{\mathrm{Var}(\hat y + \tilde\eta) = \mathrm{Var}(\hat y) +
#'     \mathrm{Var}(\tilde\eta) + 2\,\mathrm{Cov}(\hat y, \tilde\eta)}
#'     holds, so this is the cross term in the variance the adjustment
#'     restores. Weighted by \code{target_weights} when supplied.}
#'   \item{cov_eta_z}{A list with \code{mean}, \code{ci_lower},
#'     \code{ci_upper} and \code{sign_share}, each a named vector (one entry
#'     per \code{z_vars}): the bootstrapped target-sample covariance between
#'     the realised quantile gap and observed \eqn{z_k},
#'     \eqn{\mathrm{Cov}(\tilde\eta, z_k)}. This is the covariance the
#'     adjustment restores and is \eqn{\rho^*}'s numerator. It requires no
#'     assumption about \eqn{g(\mathbf{X})}: \eqn{\tilde\eta} comes from
#'     the \code{qa_fit} call on the target and \eqn{z_k} is observed
#'     there, so both are moments of the same sample. \code{sign_share} is
#'     the fraction of bootstrap draws agreeing in sign with the mean.
#'     Weighted by \code{target_weights} when supplied.}
#'   \item{linearity_check}{A list of named vectors (one entry per
#'     \code{z_vars}) flagging the two assumptions behind the
#'     predictor-channel factorisation, which fail independently.
#'     \code{gap_gbar_linear_pct} is the percentage gap between
#'     \eqn{\theta_k\,\mathrm{Cov}(\hat y, \tilde\eta)} and
#'     \eqn{\mathrm{Cov}(\tilde\eta, \hat g(\mathbf{X}))}; both use the
#'     same \eqn{\hat g}, so it isolates curvature in \eqn{\bar g}
#'     (condition (iii)) and goes to zero under joint normality of
#'     \eqn{\mathbf{X}}. \code{gap_index_reduction_pct} is the gap between
#'     \eqn{\mathrm{Cov}(\tilde\eta, \hat g(\mathbf{X}))} and
#'     \eqn{\mathrm{Cov}(\tilde\eta, z_k)} with raw \eqn{z_k}, which
#'     measures whether the linear projection stands in for
#'     \eqn{E[z_k \mid \mathbf{X}]} -- typically large for a binary
#'     \eqn{z_k} fitted by OLS. \code{cov_eta_z} is that raw-\eqn{z_k}
#'     covariance, and \code{sign_agrees} records whether the factorised
#'     channel and it share a sign. A large \code{gap_gbar_linear_pct}
#'     invalidates the channel's magnitude but not its sign, since
#'     \eqn{\mathrm{Cov}(\hat y, \tilde\eta) \ge 0}. Neither gap bears on
#'     \eqn{\lambda}, whose denominator \eqn{\mathrm{Cov}(\varepsilon,
#'     z)} is not estimable in a genuine TSTS application.}
#'   \item{theta}{A list with \code{mean}, \code{ci_lower}, \code{ci_upper}
#'     and \code{sign_share}, each a named vector (one entry per
#'     \code{z_vars}): the slope of \eqn{\bar g} on predicted income,
#'     \eqn{\theta_k = \mathrm{Cov}(g(\mathbf{X}), \hat y) /
#'     \mathrm{Var}(\hat y)}, which under the linear-projection conditions
#'     equals \eqn{\gamma'\Sigma\beta / \beta'\Sigma\beta}.
#'     \code{sign_share} is the fraction of bootstrap draws agreeing in sign
#'     with the mean, so values near 1 indicate a stable sign and values near
#'     0.5 indicate the sign is not pinned down.}
#'   \item{predictor_channel}{A list with \code{mean}, \code{ci_lower},
#'     \code{ci_upper} (the factorised channel \eqn{\theta_k\,
#'     \mathrm{Cov}(\hat y, \tilde\eta)}), plus \code{direct_mean},
#'     \code{direct_ci_lower}, \code{direct_ci_upper} (the same quantity
#'     computed directly as \eqn{\mathrm{Cov}(\tilde\eta,
#'     g(\mathbf{X}))}) and \code{sign_agrees}. The two forms coincide only
#'     when \eqn{\bar g} is linear in predicted income, which requires
#'     \eqn{\mathbf{X}} jointly normal or elliptical and generally fails
#'     with factor or bounded predictors; the gap between them measures that
#'     failure, while \code{sign_agrees} records whether the sign rule still
#'     holds. Both are oracle-free, so unlike \eqn{\lambda} they are
#'     available in a genuine TSTS application.}
#'   \item{qa_fit}{The full return value of the
#'     \code{\link{qa_fit}} call on the final predictor set, fit on the
#'     full (non-subsampled) data.}
#'
#' Use \code{plot()} on the returned object for the top \eqn{S_i(z_k)}
#' pairs, and \code{print()} for a formatted summary.
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
                                      tau = 10, B = 100, n_grid = 200,
                                      subsample_cap = 5000, donor_weights = NULL,
                                      target_weights = NULL,
                                      verbose = 2) {

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

  # --- Input validation: tau, B, n_grid ------------------------------------
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
  if (!is.numeric(n_grid) || length(n_grid) != 1 || n_grid != round(n_grid) || n_grid < 4) {
    stop("n_grid must be a single integer of at least 4.")
  }
  if (!is.numeric(subsample_cap) || length(subsample_cap) != 1 || subsample_cap < 1) {
    stop("subsample_cap must be a single positive number (use Inf for the full sample).")
  }

  # --- Input validation: donor_weights (mirrors qa_fit()'s own checks) ----
  if (!is.null(donor_weights)) {
    if (length(donor_weights) != nrow(donor_data)) {
      stop(sprintf(
        "donor_weights has length %d but donor_data has %d rows; they must match.",
        length(donor_weights), nrow(donor_data)
      ))
    }
    if (any(donor_weights < 0)) stop("donor_weights must be nonnegative.")
    if (sum(donor_weights) == 0) stop("donor_weights cannot be all zero.")
  }

  # --- Input validation: target_weights -----------------------------------
  if (!is.null(target_weights)) {
    if (length(target_weights) != nrow(target_data)) {
      stop(sprintf(
        "target_weights has length %d but target_data has %d rows; they must match.",
        length(target_weights), nrow(target_data)
      ))
    }
    if (any(target_weights < 0)) stop("target_weights must be nonnegative.")
    if (sum(target_weights) == 0) stop("target_weights cannot be all zero.")
  }

  # --- Input validation: verbose ------------------------------------------
  if (!is.numeric(verbose) || length(verbose) != 1 || !(verbose %in% c(0, 1, 2))) {
    stop("verbose must be a single value: 0 (silent), 1 (progress only), or 2 (full print).")
  }

  active <- x_vars
  removal_log <- list()
  iteration <- 0

  # ==========================================================================
  # VARIABLE SELECTION (Stage 1)
  # ==========================================================================
  if (verbose >= 1) {
    cat(sprintf("Outcome scale: %s\n", outcome_scale))
    cat("========== Variable Selection ==========\n")
  }

  repeat {
    iteration <- iteration + 1

    # (a) Variables
    if (verbose >= 2) {
      cat(sprintf("\n--- Iteration %d ---\n", iteration))
      cat(sprintf("Variables: %s\n", paste(active, collapse = ", ")))
    }

    # (b) Bootstrap progress
    S_boot <- array(NA_real_, dim = c(length(active), length(z_vars), B),
                     dimnames = list(active, z_vars, NULL))
    boot_start_time <- Sys.time()
    for (b in seq_len(B)) {
      if (verbose >= 1) {
        elapsed <- as.numeric(difftime(Sys.time(), boot_start_time, units = "secs"))
        eta <- if (b > 1) elapsed / (b - 1) * (B - (b - 1)) else NA
        .print_progress(sprintf("Bootstrap progress: %d/%d (%.0f%%) - ETA: %s",
                                 b, B, 100 * b / B, .format_duration(eta)))
      }
      if (is.null(donor_weights)) {
        donor_sub <- .draw_subsample(donor_data, subsample_cap)
        donor_weights_sub <- NULL
      } else {
        drawn <- .draw_subsample_with_weights(donor_data, donor_weights, subsample_cap)
        donor_sub <- drawn$data
        donor_weights_sub <- drawn$weights
      }
      if (is.null(target_weights)) {
        target_sub <- .draw_subsample(target_data, subsample_cap)
        target_weights_sub <- NULL
      } else {
        drawn_t <- .draw_subsample_with_weights(target_data, target_weights, subsample_cap)
        target_sub <- drawn_t$data
        target_weights_sub <- drawn_t$weights
      }
      S_boot[, , b] <- .compute_S_matrix(donor_sub, target_sub, y_var, z_vars,
                                          active, outcome_scale,
                                          donor_weights_sub, target_weights_sub)$S
    }
    if (verbose >= 1) cat("\n")

    S_mean     <- apply(S_boot, c(1, 2), mean, na.rm = TRUE)
    S_ci_lower <- apply(S_boot, c(1, 2), stats::quantile, probs = 0.025, na.rm = TRUE)
    S_ci_upper <- apply(S_boot, c(1, 2), stats::quantile, probs = 0.975, na.rm = TRUE)

    # (d) Full table of S(z_k), all (predictor, z_var) pairs, ordered highest
    # to lowest mean_S. Stored in full on the return value (not just top 5)
    # so callers can inspect every pair, not only what gets printed here.
    S_table <- data.frame(
      predictor = rep(active, times = length(z_vars)),
      z_var     = rep(z_vars, each  = length(active)),
      mean_S    = as.vector(S_mean),
      ci_lower  = as.vector(S_ci_lower),
      ci_upper  = as.vector(S_ci_upper)
    )
    S_table <- S_table[order(-S_table$mean_S), ]
    rownames(S_table) <- NULL

    if (verbose >= 2) {
      cat("Top", min(5, nrow(S_table)), "S_i(z_k) pairs (mean, 95% CI across", B, "subsamples):\n")
      print(utils::head(S_table, 5), row.names = FALSE)
    }

    # Cannot screen further if only one predictor remains
    if (length(active) == 1) {
      if (verbose >= 2) cat("\nOnly one predictor remains; stopping Stage 1.\n")
      break
    }

    # (c) If at least one passes tau
    if (max(S_mean) > tau) {
      max_idx <- which(S_mean == max(S_mean), arr.ind = TRUE)[1, ]
      i_star <- rownames(S_mean)[max_idx["row"]]
      z_star <- colnames(S_mean)[max_idx["col"]]
      max_S  <- S_mean[max_idx["row"], max_idx["col"]]

      if (verbose >= 2) {
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
      if (verbose >= 2) cat(sprintf("\nNo pair exceeds tau = %.1f. Variable selection converged.\n", tau))
      break
    }
  }

  if (verbose >= 1) {
    cat(sprintf("\nStage 1 survivors: %s\n", paste(active, collapse = ", ")))
  }

  # S_table (from the final loop iteration, above) already corresponds to
  # `active` == the final survivors, and is already ordered highest-to-
  # lowest mean_S -- returned as-is (in full, not trimmed to top 5) so
  # callers can inspect every pair. plot() builds the top-5 bar chart from
  # it on demand; see plot.qa_diagnose().

  # ==========================================================================
  # FINAL ESTIMATION: bootstrap R^2_{y,d} and rho* for the final predictors
  # ==========================================================================
  if (verbose >= 1) cat("\n========== Final Estimation ==========\n")

  r2_boot  <- rep(NA_real_, B)
  rho_boot <- matrix(NA_real_, nrow = B, ncol = length(z_vars), dimnames = list(NULL, z_vars))
  cov_yhat_eta_boot <- rep(NA_real_, B)
  theta_boot   <- matrix(NA_real_, nrow = B, ncol = length(z_vars), dimnames = list(NULL, z_vars))
  channel_boot <- matrix(NA_real_, nrow = B, ncol = length(z_vars), dimnames = list(NULL, z_vars))
  channel_direct_boot <- matrix(NA_real_, nrow = B, ncol = length(z_vars), dimnames = list(NULL, z_vars))
  cov_eta_z_boot <- matrix(NA_real_, nrow = B, ncol = length(z_vars), dimnames = list(NULL, z_vars))
  final_start_time <- Sys.time()

  for (b in seq_len(B)) {
    if (verbose >= 1) {
      elapsed <- as.numeric(difftime(Sys.time(), final_start_time, units = "secs"))
      eta <- if (b > 1) elapsed / (b - 1) * (B - (b - 1)) else NA
      .print_progress(sprintf("Bootstrap progress: %d/%d (%.0f%%) - ETA: %s",
                               b, B, 100 * b / B, .format_duration(eta)))
    }

    if (is.null(donor_weights)) {
      donor_sub <- .draw_subsample(donor_data, subsample_cap)
      donor_weights_sub <- NULL
    } else {
      drawn <- .draw_subsample_with_weights(donor_data, donor_weights, subsample_cap)
      donor_sub <- drawn$data
      donor_weights_sub <- drawn$weights
    }
    if (is.null(target_weights)) {
      target_sub <- .draw_subsample(target_data, subsample_cap)
      target_weights_sub <- NULL
    } else {
      drawn_t <- .draw_subsample_with_weights(target_data, target_weights, subsample_cap)
      target_sub <- drawn_t$data
      target_weights_sub <- drawn_t$weights
    }

    fit_sub <- .fit_first_stage_r2(donor_sub, y_var, active, outcome_scale, donor_weights_sub)
    r2_boot[b] <- fit_sub$r2
    eps_d_sub  <- fit_sub$y_model - fit_sub$y_hat

    # z_k ~ x_vars in the target sample, weighted by target_weights so that
    # z_perp is orthogonal to the predictors in the SAME metric the
    # covariances below are computed in. Without this, Cov_weighted(eps, z)
    # and Cov_weighted(eps, z_perp) diverge and rho* stops satisfying
    # lambda = rho*/rho.
    tw_fit <- if (is.null(target_weights_sub)) rep(1, nrow(target_sub)) else target_weights_sub
    z_perp_sub <- sapply(z_vars, function(z_k) {
      f_z <- stats::reformulate(active, response = z_k)
      mod_z <- stats::lm(f_z, data = target_sub, weights = tw_fit)
      stats::residuals(mod_z)
    })
    colnames(z_perp_sub) <- z_vars

    # This subsample's donor/target could, by chance, have a level mismatch
    # even when the FULL donor/target agree overall (checked upfront via
    # .check_factor_levels) -- tolerate that here by skipping this
    # replicate's rho* (left NA) rather than crashing the whole run.
    qa_sub <- tryCatch(
      qa_fit(donor_sub, target_sub, y_var, active, outcome_scale = outcome_scale,
             n_grid = n_grid, donor_weights = donor_weights_sub),
      error = function(e) NULL
    )
    if (is.null(qa_sub)) next
    eta_sub <- qa_sub$eta_target

    dw <- if (is.null(donor_weights_sub))  rep(1, nrow(donor_sub))  else donor_weights_sub
    tw <- if (is.null(target_weights_sub)) rep(1, nrow(target_sub)) else target_weights_sub

    rho_boot[b, ] <- sapply(z_vars, function(z_k) {
      .wcov(eta_sub, target_sub[[z_k]], tw) /
        sqrt(.wvar(eps_d_sub, dw) * .wvar(z_perp_sub[, z_k], tw))
    })

    # Cov(y_hat, eta) in the target sample, on the MODEL scale (log scale
    # when outcome_scale = "log"), which is the scale on which the variance
    # decomposition Var(y_hat + eta) = Var(y_hat) + Var(eta) + 2Cov(y_hat,
    # eta) applies. A value near zero means the adjustment adds variance
    # roughly additively; a large negative value means the correction is
    # partly offsetting the prediction rather than adding to it. Weighted by
    # target_weights, since this is a target-sample population moment.
    cov_yhat_eta_boot[b] <- .wcov(qa_sub$y_hat_target, eta_sub, tw)

    # --- Linear-projection diagnostics (appendix "The Linear Projection
    # Case") ------------------------------------------------------------
    # Under f(X) = beta'X, g(X) = gamma'X and a linear gbar, the slope of
    # gbar on predicted income is
    #     theta_k = gamma' Sigma beta / beta' Sigma beta
    #             = Cov(g(X), y_hat) / Var(y_hat),
    # and the predictor channel factorises as
    #     Cov(eta, g(X)) = theta_k * Cov(y_hat, eta).
    # Both are computed here from quantities the loop already has:
    # g_hat = z_k - z_perp is the fitted part of z_k on the predictors
    # (i.e. gamma'X), and y_hat comes from the same qa_fit() call as eta.
    #
    # This matters because it is ORACLE-FREE. Unlike lambda, whose
    # denominator Cov(epsilon, z) is unobservable in a genuine TSTS
    # application, theta_k and Cov(y_hat, eta) are both estimable from the
    # donor and target samples alone. Since Cov(y_hat, eta) >= 0, the SIGN
    # of the channel is the sign of theta_k, so sign(theta_k) predicts the
    # direction in which the adjustment moves Cov(eta, z) before any
    # validation data is seen. Mis-recovery is exactly the case where that
    # sign opposes the sign of the omitted covariance.
    #
    # Weighted by target_weights throughout: Sigma, Cov(g(X), y_hat) and
    # Var(y_hat) are all target-population moments.
    var_yhat_sub <- .wvar(qa_sub$y_hat_target, tw)
    theta_boot[b, ] <- sapply(z_vars, function(z_k) {
      g_hat <- target_sub[[z_k]] - z_perp_sub[, z_k]
      .wcov(g_hat, qa_sub$y_hat_target, tw) / var_yhat_sub
    })
    channel_boot[b, ] <- theta_boot[b, ] * cov_yhat_eta_boot[b]

    # The channel computed DIRECTLY as Cov(eta, g(X)), without going through
    # theta. The two agree exactly only under condition (iii) of the
    # appendix -- gbar linear in predicted income -- which needs X jointly
    # normal (or elliptical) and generally fails with factor and bounded
    # predictors. Returning both makes the size of that failure measurable:
    # theta captures only the component of gbar that is linear in y_hat, so
    # the gap between channel_direct and predictor_channel is the part of
    # gbar orthogonal to y_hat. Signs can be expected to agree even when
    # magnitudes do not, which is what the sign rule of Equation
    # (eq:app_sign) actually needs.
    channel_direct_boot[b, ] <- sapply(z_vars, function(z_k) {
      g_hat <- target_sub[[z_k]] - z_perp_sub[, z_k]
      .wcov(eta_sub, g_hat, tw)
    })

    # Cov(eta, z) against RAW z, the quantity rho*'s numerator uses. Needed
    # for the second of the two linearity gaps below: comparing it to
    # channel_direct (which uses the fitted g_hat instead of z) isolates
    # whether the linear projection of z on the predictors stands in for
    # E[z | X]. For a binary z_k fitted by OLS it generally does not.
    cov_eta_z_boot[b, ] <- sapply(z_vars, function(z_k) {
      .wcov(eta_sub, target_sub[[z_k]], tw)
    })
  }
  if (verbose >= 1) cat("\n")

  R2_y_donor <- list(
    mean     = mean(r2_boot, na.rm = TRUE),
    ci_lower = unname(stats::quantile(r2_boot, 0.025, na.rm = TRUE)),
    ci_upper = unname(stats::quantile(r2_boot, 0.975, na.rm = TRUE))
  )

  rho_star <- list(
    mean     = colMeans(rho_boot, na.rm = TRUE),
    ci_lower = apply(rho_boot, 2, stats::quantile, probs = 0.025, na.rm = TRUE),
    ci_upper = apply(rho_boot, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
  )

  cov_yhat_eta <- list(
    mean     = mean(cov_yhat_eta_boot, na.rm = TRUE),
    ci_lower = unname(stats::quantile(cov_yhat_eta_boot, 0.025, na.rm = TRUE)),
    ci_upper = unname(stats::quantile(cov_yhat_eta_boot, 0.975, na.rm = TRUE))
  )

  # theta: slope of gbar on predicted income, per z_k. sign_share is the
  # fraction of bootstrap draws whose theta agrees in sign with the mean --
  # a value near 1 means the predicted direction of the predictor channel is
  # stable, near 0.5 means the sign is not pinned down by the data.
  theta <- list(
    mean       = colMeans(theta_boot, na.rm = TRUE),
    ci_lower   = apply(theta_boot, 2, stats::quantile, probs = 0.025, na.rm = TRUE),
    ci_upper   = apply(theta_boot, 2, stats::quantile, probs = 0.975, na.rm = TRUE),
    sign_share = sapply(z_vars, function(z_k) {
      v <- theta_boot[, z_k]
      v <- v[!is.na(v)]
      if (length(v) == 0) return(NA_real_)
      mean(sign(v) == sign(mean(v)))
    })
  )

  # predictor_channel: theta_k * Cov(y_hat, eta), the appendix's prediction
  # for Cov(eta, g(X)). Oracle-free, so its sign is available before any
  # validation data.
  # linearity_check: the two assumptions behind the appendix's factorisation
  # of the predictor channel, each with its own gap, measured separately
  # because they fail independently.
  #
  #   gap_gbar_linear_pct -- condition (iii), gbar linear in predicted
  #     income. Compares theta * Cov(y_hat, eta) against Cov(eta, g_hat).
  #     BOTH sides use the same g_hat, so the only thing that can separate
  #     them is curvature in gbar. Goes to zero when X is jointly normal
  #     (or elliptical) and stays near zero even when z is binary, since
  #     that breaks a different assumption. A small gap licenses reading the
  #     factorised channel as a magnitude; a large one does not, but leaves
  #     the SIGN rule intact, since Cov(y_hat, eta) >= 0 makes
  #     sign(channel) = sign(theta) whatever the shape of gbar.
  #
  #   gap_index_reduction_pct -- whether the linear projection g_hat stands
  #     in for E[z | X]. Compares Cov(eta, g_hat) against Cov(eta, z) with
  #     raw z. The index reduction Cov(eta, z) = Cov(eta, E[z|X]) is exact
  #     when eta is X-measurable, but g_hat is the LINEAR projection, not
  #     the conditional expectation; for a binary z_k fitted by OLS the two
  #     differ. Large values here suggest fitting z_k ~ x_vars by logit
  #     rather than OLS if the channel magnitude is to be reported.
  #
  # Gaps are computed from the bootstrap means rather than averaged over
  # per-draw ratios, which would be unstable whenever a denominator is near
  # zero. NA where the denominator is ~0 (no meaningful relative gap).
  .rel_gap <- function(a, b) {
    ifelse(abs(b) < .Machine$double.eps^0.5, NA_real_, 100 * abs(a - b) / abs(b))
  }
  channel_lin_mean    <- colMeans(channel_boot, na.rm = TRUE)
  channel_direct_mean <- colMeans(channel_direct_boot, na.rm = TRUE)
  cov_eta_z_mean      <- colMeans(cov_eta_z_boot, na.rm = TRUE)

  linearity_check <- list(
    cov_eta_z               = cov_eta_z_mean,
    gap_gbar_linear_pct     = .rel_gap(channel_lin_mean, channel_direct_mean),
    gap_index_reduction_pct = .rel_gap(channel_direct_mean, cov_eta_z_mean),
    sign_agrees             = sign(channel_lin_mean) == sign(cov_eta_z_mean)
  )

  # Cov(eta, z): the covariance the adjustment actually restores, measured
  # directly against observed z_k. This needs no assumption about g(X) at
  # all -- eta comes from the qa_fit() call on the target and z_k is
  # observed in the target, so both sides are sample moments of the same
  # sample. It is also rho*'s numerator. sign_share is the fraction of
  # bootstrap draws agreeing in sign with the mean: 1 means every draw
  # agreed on the direction, 0.5 means the sign is not pinned down.
  cov_eta_z <- list(
    mean       = colMeans(cov_eta_z_boot, na.rm = TRUE),
    ci_lower   = apply(cov_eta_z_boot, 2, stats::quantile, probs = 0.025, na.rm = TRUE),
    ci_upper   = apply(cov_eta_z_boot, 2, stats::quantile, probs = 0.975, na.rm = TRUE),
    sign_share = sapply(z_vars, function(z_k) {
      v <- cov_eta_z_boot[, z_k]
      v <- v[!is.na(v)]
      if (length(v) == 0) return(NA_real_)
      mean(sign(v) == sign(mean(v)))
    })
  )

  predictor_channel <- list(
    mean     = colMeans(channel_boot, na.rm = TRUE),
    ci_lower = apply(channel_boot, 2, stats::quantile, probs = 0.025, na.rm = TRUE),
    ci_upper = apply(channel_boot, 2, stats::quantile, probs = 0.975, na.rm = TRUE),
    direct_mean = colMeans(channel_direct_boot, na.rm = TRUE),
    direct_ci_lower = apply(channel_direct_boot, 2, stats::quantile, probs = 0.025, na.rm = TRUE),
    direct_ci_upper = apply(channel_direct_boot, 2, stats::quantile, probs = 0.975, na.rm = TRUE),
    sign_agrees = sapply(z_vars, function(z_k)
      sign(mean(channel_boot[, z_k], na.rm = TRUE)) ==
        sign(mean(channel_direct_boot[, z_k], na.rm = TRUE)))
  )

  if (verbose >= 2) {
    cat(sprintf("R^2_y,d: mean = %.4f, 95%% CI [%.4f, %.4f]\n",
                R2_y_donor$mean, R2_y_donor$ci_lower, R2_y_donor$ci_upper))
    cat("\nrho*:\n")
    rho_table <- data.frame(
      z_var    = z_vars,
      mean     = round(rho_star$mean, 4),
      ci_lower = round(rho_star$ci_lower, 4),
      ci_upper = round(rho_star$ci_upper, 4)
    )
    print(rho_table, row.names = FALSE)
  }

  # ==========================================================================
  # Final fit on the FULL data, using the selected specification
  # ==========================================================================
  qa_result <- qa_fit(donor_data, target_data, y_var, active,
                       outcome_scale = outcome_scale, n_grid = n_grid,
                       donor_weights = donor_weights)

  structure(
    list(
      selected_predictors = active,
      removed_predictors  = removal_log,
      S_table             = S_table,
      R2_y_donor          = R2_y_donor,
      rho_star            = rho_star,
      cov_yhat_eta        = cov_yhat_eta,
      cov_eta_z           = cov_eta_z,
      theta               = theta,
      predictor_channel   = predictor_channel,
      linearity_check     = linearity_check,
      qa_fit              = qa_result
    ),
    class = "qa_diagnose"
  )
}

#' Plot the Top S(z_k) Pairs From a \code{qa_diagnose} Result
#'
#' A horizontal bar chart of the top \code{n} \eqn{S_i(z_k)} pairs (mean
#' across bootstrap subsamples) among the final survivors, with reference
#' lines at \eqn{\tau = 5} and \eqn{\tau = 10}. Built on demand from the
#' full table already stored in \code{x$S_table}, not computed inside
#' \code{\link{qa_diagnose}} itself -- use \code{x$S_table} directly for
#' every pair, not just the ones shown here.
#'
#' @param x A list returned by \code{\link{qa_diagnose}}.
#' @param n Integer, number of top pairs (by mean \eqn{S}) to show. Default 5.
#' @param ... Ignored; present for S3 method consistency.
#'
#' @return A \code{ggplot} object. Requires the \code{ggplot2} package to be
#'   installed; \code{ggplot2} is listed under \code{Suggests} rather than
#'   \code{Imports} so it is not a hard dependency for users who never plot.
#'
#' @export
plot.qa_diagnose <- function(x, n = 5, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot.qa_diagnose() requires the 'ggplot2' package. ",
         "Install it with install.packages('ggplot2').")
  }

  plot_df <- utils::head(x$S_table, n)
  plot_df$label <- paste0(plot_df$z_var, "\n", plot_df$predictor)
  plot_df$label <- factor(plot_df$label, levels = rev(plot_df$label))

  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$label, y = .data$mean_S,
                                        fill = .data$z_var)) +
    ggplot2::geom_col(width = 0.65) +
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

#' Print a \code{qa_diagnose} Result
#'
#' A formatted summary: the final selected predictors, the bootstrapped
#' \eqn{R^2_{y,d}} and \eqn{\rho^*} with their 95\% CIs.
#'
#' @param x A list returned by \code{\link{qa_diagnose}}.
#' @param ... Ignored; present for S3 method consistency.
#'
#' @export
print.qa_diagnose <- function(x, ...) {
  cat("<qa_diagnose result>\n")
  cat(sprintf("  Selected predictors: %s\n", paste(x$selected_predictors, collapse = ", ")))
  cat(sprintf("  Predictors removed:  %d\n", length(x$removed_predictors)))
  cat(sprintf("  R^2_y,d: mean = %.4f, 95%% CI [%.4f, %.4f]\n",
              x$R2_y_donor$mean, x$R2_y_donor$ci_lower, x$R2_y_donor$ci_upper))
  # Guard for objects created before cov_yhat_eta was added to the return
  # value, so print() still works on a saved result from an older version.
  if (!is.null(x$cov_yhat_eta)) {
    cat(sprintf("  Cov(y_hat, eta): mean = %.4f, 95%% CI [%.4f, %.4f]\n",
                x$cov_yhat_eta$mean, x$cov_yhat_eta$ci_lower, x$cov_yhat_eta$ci_upper))
  }
  cat("\n  rho*:\n")
  rho_table <- data.frame(
    z_var    = names(x$rho_star$mean),
    mean     = round(x$rho_star$mean, 4),
    ci_lower = round(x$rho_star$ci_lower, 4),
    ci_upper = round(x$rho_star$ci_upper, 4)
  )
  print(rho_table, row.names = FALSE)

  # Cov(eta, z): the covariance the adjustment restores, measured directly
  # against observed z_k. No assumption about g(X) enters. theta, the
  # factorised predictor_channel and linearity_check are still computed and
  # available on the object for the linear-case analysis, but are not
  # printed -- Cov(eta, z) is what the adjustment is judged on. Guarded so
  # results saved before this component existed still print.
  if (!is.null(x$cov_eta_z)) {
    cat("\n  Cov(eta, z):\n")
    cov_table <- data.frame(
      z_var      = names(x$cov_eta_z$mean),
      mean       = round(x$cov_eta_z$mean, 5),
      ci_lower   = round(x$cov_eta_z$ci_lower, 5),
      ci_upper   = round(x$cov_eta_z$ci_upper, 5),
      sign_share = round(x$cov_eta_z$sign_share, 2)
    )
    print(cov_table, row.names = FALSE)
  }
  cat("\nUse plot(x) for the top S(z_k) pairs, x$S_table for every pair.\n")
  invisible(x)
}
