# --- Small internal helper: wrap the raw smooth.spline as a plain callable,
# eta_spline(p), rather than requiring predict(eta_spline, p)$y.
#
# A windowed-OLS linear override for the sparse boundary region was tried
# and reverted: it was built to fix R-vs-Python cross-implementation tail
# instability (no longer a live concern -- the Python port is shelved), but
# was found, via Monte Carlo validation, to introduce a systematic
# DOWNWARD bias in recovered variance -- replacing a genuinely curved
# (accelerating) tail shape with a flatter windowed-average linear
# approximation understates how extreme the correction should be exactly
# where it matters most for the distribution's variance. Plain spline
# extrapolation (this version) was empirically closer to unbiased across
# repeated simulation, at the cost of occasional instability in individual
# replications -- a trade-off judged acceptable now that cross-language
# matching is not a goal.
.build_eta_function <- function(eta_spline_raw) {
  function(p) stats::predict(eta_spline_raw, p)$y
}

#' Quantile Adjustment for Two-Sample Two-Stage (TSTS) Imputation
#'
#' Fits a first-stage prediction model in the donor sample, constructs the
#' empirical quantile-gap function \eqn{\hat\eta(p)} between observed and
#' predicted outcomes, and applies the corresponding position-specific
#' correction to the target sample. This implements the quantile adjustment
#' of Section 4.2 (population level) and Section 5.1 (estimation) of
#' Torres-Lopez, "Bias Reduction in TSTSLS Applications: A Quantile
#' Adjustment Method."
#'
#' @details
#' The choice of \code{outcome_scale} determines both the estimation scale
#' and the first-stage model family, not just a cosmetic transform:
#' \itemize{
#'   \item \code{"log"}: the outcome is log-transformed before fitting an
#'     ordinary least squares model, \code{lm(log(y) ~ X)}. The quantile
#'     gap is estimated in log space and the final adjusted values are
#'     exponentiated back to levels.
#'   \item \code{"level"}: the outcome is left untransformed and the first
#'     stage is a Gaussian GLM with a log link,
#'     \code{glm(y ~ X, family = gaussian(link = "log"))}. Predictions are
#'     taken on the response (level) scale via \code{type = "response"}, so
#'     the quantile gap is estimated directly in levels and no
#'     back-transformation is applied. This distinction matters because
#'     the additive variance-restoration property
#'     (\eqn{\mathrm{Var}(\tilde y_q) = \mathrm{Var}(y)}) is derived on
#'     whatever scale the adjustment is applied on; exponentiating a
#'     log-scale additive adjustment does not preserve it in levels.
#' }
#'
#' The estimation grid runs from \eqn{1/(n\_grid+1)} to
#' \eqn{n\_grid/(n\_grid+1)}; a target observation whose predicted rank
#' falls outside this range triggers a common-support warning (Section
#' 5.2) and its adjustment relies on \code{smooth.spline}'s own
#' extrapolation beyond its fitted domain. An earlier version of this
#' function overrode that extrapolation with a windowed linear
#' approximation, built to keep results reproducible against an
#' (since-shelved) Python port; Monte Carlo validation found that override
#' introduced a systematic downward bias in recovered variance, since a
#' linear approximation understates the genuinely accelerating shape
#' \eqn{\hat\eta(p)} takes near the extremes. Plain spline extrapolation
#' (the current behavior) was empirically closer to unbiased across
#' repeated simulation, at the cost of occasional instability in
#' individual replications.
#'
#' Weighting (via \code{donor_weights}) enters only through the empirical
#' quantile estimates used to build \eqn{\hat\eta(p)}. It is the caller's
#' responsibility to have already computed weights appropriate for aligning
#' the donor's covariate distribution with the target's (e.g., via inverse
#' probability weighting or entropy balancing); this function does not
#' estimate weights itself. The spline interpolating the quantile gap is
#' always fit unweighted, since the weighting has already entered through
#' the quantile step.
#'
#' @param donor_data Data frame containing the donor sample. Must include
#'   \code{y_var} and all \code{x_vars}, with no missing values in these
#'   columns.
#' @param target_data Data frame containing the target sample. Must include
#'   all \code{x_vars}, with no missing values in these columns. \code{y_var}
#'   is not required (and is not used) in \code{target_data}, since the
#'   outcome is by construction unobserved in the target sample. For any
#'   \code{x_vars} column that is a factor or character, every level present
#'   in \code{target_data} must also appear in \code{donor_data} (checked
#'   upfront and rejected with \code{stop()} otherwise, since a level never
#'   seen during fitting cannot be predicted); a level present in
#'   \code{donor_data} but absent from \code{target_data} triggers a
#'   \code{warning()} instead, since it does not prevent prediction.
#' @param y_var Character string naming the outcome column in
#'   \code{donor_data}.
#' @param x_vars Character vector naming the predictor columns, present in
#'   both \code{donor_data} and \code{target_data}.
#' @param outcome_scale Either \code{"log"} or \code{"level"}; see Details.
#' @param n_grid Integer, number of interior quantile grid points in
#'   \eqn{(0,1)} used to estimate \eqn{\hat\eta(p)}. Default 200, matching
#'   the paper's applications.
#' @param donor_weights Optional numeric vector, length \code{nrow(donor_data)},
#'   of nonnegative weights used only when computing the weighted empirical
#'   quantiles in the donor sample. Defaults to \code{NULL}, i.e. uniform
#'   weights (unweighted quantiles). Weights need not sum to any particular
#'   value; any positive common scale is equivalent, since quantiles are
#'   scale-invariant to rescaling of the weights.
#'
#' @return A list (class \code{"qa_fit"}) with components:
#'   \item{model}{The fitted first-stage model object (\code{lm} or \code{glm}).}
#'   \item{donor_resid}{Donor in-sample residuals on the model scale.}
#'   \item{eta_spline}{A callable function, \code{eta_spline(p)}, giving the
#'     estimated quantile-gap value at position(s) \code{p}, via the fitted
#'     smoothing spline (using its own extrapolation beyond the estimation
#'     grid; see Details).}
#'   \item{donor_ecdf}{A function giving the weighted empirical CDF of donor
#'     predicted values, used to assign target quantile positions.}
#'   \item{p_grid}{The interior quantile grid in \eqn{(0,1)} used to estimate
#'     \eqn{\hat\eta(p)}.}
#'   \item{Q_y_d}{Weighted empirical quantiles of the observed donor outcome
#'     (model scale) at \code{p_grid}.}
#'   \item{Q_yhat_d}{Weighted empirical quantiles of the donor's predicted
#'     outcome (model scale) at \code{p_grid}.}
#'   \item{p_target}{Quantile position \eqn{p_t} assigned to each target
#'     observation.}
#'   \item{eta_target}{Realised quantile-gap adjustment \eqn{\tilde\eta_t} for
#'     each target observation.}
#'   \item{y_hat_target}{First-stage predicted values in the target sample,
#'     on the model scale (log scale if \code{outcome_scale = "log"}, level
#'     scale if \code{"level"}).}
#'   \item{y_adjusted}{Final quantile-adjusted outcome for the target sample,
#'     always on the level (original) scale.}
#'
#' Use \code{plot()} on the returned object for the \eqn{\hat\eta(p)} curve,
#' and \code{print()} for a formatted one-line summary.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' n <- 5000
#' X <- rnorm(n)
#' eps <- rnorm(n, sd = 1.5)
#' y <- exp(1 + 0.8 * X + eps)
#' donor <- data.frame(y = y[1:2500], X = X[1:2500])
#' target <- data.frame(X = X[2501:5000])
#' result <- qa_fit(donor, target, y_var = "y", x_vars = "X",
#'                                outcome_scale = "log")
#' }
#'
#' @importFrom rlang .data
#' @export
qa_fit <- function(donor_data, target_data, y_var, x_vars,
                                 outcome_scale = c("log", "level"),
                                 n_grid = 200,
                                 donor_weights = NULL) {

  outcome_scale <- match.arg(outcome_scale)

  # --- Input validation: basic types --------------------------------------
  if (!is.data.frame(donor_data)) {
    stop("donor_data must be a data frame.")
  }
  if (!is.data.frame(target_data)) {
    stop("target_data must be a data frame.")
  }
  if (!is.character(y_var) || length(y_var) != 1) {
    stop("y_var must be a single character string naming a column of donor_data.")
  }
  if (!is.character(x_vars) || length(x_vars) < 1) {
    stop("x_vars must be a character vector of at least one column name.")
  }

  # --- Input validation: required columns actually present ----------------
  # Referencing a nonexistent column would otherwise fail deep inside
  # reformulate()/lm()/glm() with a much less informative error.
  missing_in_donor <- setdiff(c(y_var, x_vars), names(donor_data))
  if (length(missing_in_donor) > 0) {
    stop(sprintf("donor_data is missing column(s): %s",
                 paste(missing_in_donor, collapse = ", ")))
  }
  missing_in_target <- setdiff(x_vars, names(target_data))
  if (length(missing_in_target) > 0) {
    stop(sprintf("target_data is missing column(s): %s",
                 paste(missing_in_target, collapse = ", ")))
  }
  if (!is.numeric(donor_data[[y_var]])) {
    stop(sprintf("donor_data[['%s']] must be numeric.", y_var))
  }

  # --- Input validation: n_grid --------------------------------------------
  # smooth.spline() requires at least 4 unique x values to fit a meaningful
  # curve; below that the "smoothing" spline is not doing anything sensible.
  if (!is.numeric(n_grid) || length(n_grid) != 1 || n_grid != round(n_grid) || n_grid < 4) {
    stop("n_grid must be a single integer of at least 4.")
  }

  # --- Input validation: missing values ---------------------------------
  # Silently dropping rows with NA would silently change the sample the
  # quantile gap is estimated on; we stop instead and tell the user where
  # the problem is.
  donor_cols_needed <- c(y_var, x_vars)
  donor_na <- !stats::complete.cases(donor_data[, donor_cols_needed, drop = FALSE])
  if (any(donor_na)) {
    stop(sprintf(
      "donor_data contains missing values in %d row(s) among columns: %s. ",
      sum(donor_na), paste(donor_cols_needed, collapse = ", ")
    ), "Remove or impute these before calling qa_fit().")
  }

  target_na <- !stats::complete.cases(target_data[, x_vars, drop = FALSE])
  if (any(target_na)) {
    stop(sprintf(
      "target_data contains missing values in %d row(s) among columns: %s. ",
      sum(target_na), paste(x_vars, collapse = ", ")
    ), "Remove or impute these before calling qa_fit().")
  }

  # --- Input validation: donor/target categorical level consistency -------
  # A level present in target but absent from donor is guaranteed to crash
  # predict() the moment that row is reached (R has no coefficient for a
  # level it never saw during fitting). Checked once, upfront, rather than
  # letting this surface as a cryptic error mid-computation.
  .check_factor_levels(donor_data, target_data, x_vars)

  # --- Input validation: outcome positivity under log scale --------------
  # log() of a non-positive value produces NaN/-Inf, which would silently
  # corrupt every downstream quantile and spline computation.
  if (outcome_scale == "log" && any(donor_data[[y_var]] <= 0)) {
    stop(sprintf(
      "outcome_scale = 'log' requires strictly positive values of '%s' in donor_data, ",
      y_var
    ), "but non-positive values were found. Use outcome_scale = 'level' instead, ",
    "or transform the outcome before calling this function.")
  }

  # --- Input validation: donor weights ------------------------------------
  if (is.null(donor_weights)) {
    donor_weights <- rep(1, nrow(donor_data))
  } else {
    if (length(donor_weights) != nrow(donor_data)) {
      stop(sprintf(
        "donor_weights has length %d but donor_data has %d rows; they must match.",
        length(donor_weights), nrow(donor_data)
      ))
    }
    if (any(donor_weights < 0)) {
      stop("donor_weights must be nonnegative.")
    }
    if (sum(donor_weights) == 0) {
      stop("donor_weights cannot be all zero.")
    }
  }

  # --- Step 1-2: fit first-stage model on the appropriate scale -----------
  # Log scale: ordinary least squares on log(y). Level scale: Gaussian GLM
  # with a log link, so predictions on the response scale directly target
  # E[y | X] in levels and the additive variance-restoration property holds
  # without a back-transform (see paper's Section 3.1 footnote).
  #
  # donor_weights is intentionally NOT passed to this fit: cross-checking
  # on real survey data found weighting the model fit made no measurable
  # difference to recovery, and the source of a systematic bias traced
  # entirely to the quantile step below. donor_weights is used there only.
  if (outcome_scale == "log") {
    y_model_d <- log(donor_data[[y_var]])
    fit_data <- cbind(donor_data, y_model_d)
    model_formula <- stats::reformulate(x_vars, response = "y_model_d")
    hat_f <- stats::lm(model_formula, data = fit_data)
    y_hat_d <- stats::predict(hat_f, newdata = donor_data)
  } else {
    y_model_d <- donor_data[[y_var]]
    fit_data <- cbind(donor_data, y_model_d)
    model_formula <- stats::reformulate(x_vars, response = "y_model_d")
    hat_f <- stats::glm(model_formula, data = fit_data,
                         family = stats::gaussian(link = "log"))
    y_hat_d <- stats::predict(hat_f, newdata = donor_data, type = "response")
  }

  # --- Step 3: quantile grid on the open interval (0, 1) -------------------
  # Interior-only (never exactly 0 or 1): including the literal sample
  # min/max as spline-fitting knots was tested and found to make things
  # worse (forcing the spline through a single noisy extreme point distorts
  # its shape trying to accommodate that one observation).
  p_grid <- seq(1 / (n_grid + 1), n_grid / (n_grid + 1), length.out = n_grid)

  # --- Step 4: weighted empirical quantile gap at each grid point ----------
  Q_y_d    <- Hmisc::wtd.quantile(y_model_d, weights = donor_weights, probs = p_grid)
  Q_yhat_d <- Hmisc::wtd.quantile(y_hat_d,   weights = donor_weights, probs = p_grid)
  eta_grid <- Q_y_d - Q_yhat_d

  # --- Step 5: smoothing spline over the grid (always unweighted) ---------
  eta_spline_raw <- stats::smooth.spline(p_grid, eta_grid)
  eta_spline <- .build_eta_function(eta_spline_raw)

  # --- Step 6: weighted empirical CDF of donor predicted values ------------
  ord   <- order(y_hat_d)
  cum_w <- cumsum(donor_weights[ord]) / sum(donor_weights)
  F_yhat_d <- stats::approxfun(y_hat_d[ord], cum_w,
                                method = "constant", yleft = 0, yright = 1,
                                ties = "ordered")

  # --- Step 7: apply to target sample --------------------------------------
  if (outcome_scale == "log") {
    y_hat_t <- stats::predict(hat_f, newdata = target_data)
  } else {
    y_hat_t <- stats::predict(hat_f, newdata = target_data, type = "response")
  }
  p_t <- F_yhat_d(y_hat_t)

  # Common-support check (Section 5.2): if target predictions fall outside
  # the donor's predicted range, p_t is pushed to the boundary and eta_hat(p_t)
  # relies on spline extrapolation beyond the estimation grid. We warn and
  # proceed rather than stop, per the agreed design.
  grid_range <- range(p_grid)
  out_of_support <- p_t <= grid_range[1] | p_t >= grid_range[2]
  if (any(out_of_support)) {
    warning(sprintf(
      "%d target observation(s) fall outside the donor's predicted support ",
      sum(out_of_support)
    ), "(common support violation, Section 5.2). The quantile-gap correction ",
    "for these observations relies on extrapolation and may be unreliable.")
  }

  eta_t <- eta_spline(p_t)
  y_tilde_model <- y_hat_t + eta_t

  # --- Step 8: back-transform only when estimation was done in log space --
  y_adjusted <- if (outcome_scale == "log") exp(y_tilde_model) else y_tilde_model

  structure(
    list(
      model        = hat_f,
      donor_resid  = y_model_d - y_hat_d,
      eta_spline   = eta_spline,
      donor_ecdf   = F_yhat_d,
      p_grid       = p_grid,
      Q_y_d        = Q_y_d,
      Q_yhat_d     = Q_yhat_d,
      p_target     = p_t,
      eta_target   = eta_t,
      y_hat_target = y_hat_t,
      y_adjusted   = y_adjusted
    ),
    class = "qa_fit"
  )
}

#' Plot the Estimated Quantile Gap From a \code{qa_fit} Result
#'
#' Shows the raw (noisy) grid-point estimates of
#' \eqn{\hat\eta(p) = Q_y(p) - Q_{\hat y}(p)} as points, against the fitted
#' smoothing spline as a line -- this is the object that actually gets used
#' to adjust the target sample, so plotting it directly (rather than the two
#' underlying quantile-function curves) shows what the correction itself
#' looks like across the distribution. Built on demand from the components
#' already stored in \code{x} (\code{p_grid}, \code{Q_y_d}, \code{Q_yhat_d},
#' \code{eta_spline}), not computed inside \code{\link{qa_fit}} itself.
#'
#' @param x A list returned by \code{\link{qa_fit}}.
#' @param annotate_p Optional single numeric in \eqn{(0,1)}. If supplied,
#'   the plot marks \eqn{\hat\eta(\code{annotate\_p})} (read off the fitted
#'   spline) with a point, a dropline to zero, and a text label. Defaults to
#'   \code{NULL} (no annotation).
#' @param ... Ignored; present for S3 method consistency.
#'
#' @return A \code{ggplot} object. Requires the \code{ggplot2} package to be
#'   installed; \code{ggplot2} is listed under \code{Suggests} rather than
#'   \code{Imports} so it is not a hard dependency for users who never plot.
#'
#' @export
plot.qa_fit <- function(x, annotate_p = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot.qa_fit() requires the 'ggplot2' package. ",
         "Install it with install.packages('ggplot2').")
  }
  if (!is.null(annotate_p)) {
    if (!is.numeric(annotate_p) || length(annotate_p) != 1 ||
        annotate_p <= 0 || annotate_p >= 1) {
      stop("annotate_p must be a single numeric value strictly between 0 and 1.")
    }
  }

  p_grid <- x$p_grid
  eta_grid <- x$Q_y_d - x$Q_yhat_d
  eta_spline <- x$eta_spline

  raw_data <- data.frame(p = p_grid, eta = eta_grid)
  p_fine <- seq(min(p_grid), max(p_grid), length.out = 500)
  smooth_data <- data.frame(p = p_fine, eta = eta_spline(p_fine))

  eta_plot <- ggplot2::ggplot() +
    ggplot2::geom_point(data = raw_data, ggplot2::aes(x = .data$p, y = .data$eta),
                         alpha = 0.35, size = 1, color = "#4C72B0") +
    ggplot2::geom_line(data = smooth_data, ggplot2::aes(x = .data$p, y = .data$eta),
                        color = "#0072B2", linewidth = 1.1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "gray60") +
    ggplot2::labs(x = "Quantile position p", y = expression(hat(eta)(p))) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.background = ggplot2::element_rect(fill = "white", color = NA),
                   panel.background = ggplot2::element_rect(fill = "white", color = NA))

  if (!is.null(annotate_p)) {
    if (annotate_p < min(p_grid) || annotate_p > max(p_grid)) {
      warning(sprintf(
        "annotate_p = %.3f falls outside the estimated grid range [%.3f, %.3f]; ",
        annotate_p, min(p_grid), max(p_grid)
      ), "the annotated value relies on spline extrapolation and may be unreliable.")
    }
    eta_at_p <- eta_spline(annotate_p)
    point_data <- data.frame(p = annotate_p, eta = eta_at_p)

    eta_plot <- eta_plot +
      ggplot2::geom_segment(data = point_data,
                             ggplot2::aes(x = .data$p, xend = .data$p, y = 0, yend = .data$eta),
                             linetype = "dashed", color = "gray40") +
      ggplot2::geom_point(data = point_data, ggplot2::aes(x = .data$p, y = .data$eta),
                           color = "#D55E00", size = 2.8) +
      ggplot2::annotate("text", x = annotate_p, y = eta_at_p,
                         label = sprintf("eta(%.2f) = %.3f", annotate_p, eta_at_p),
                         vjust = -0.8, size = 3.6, color = "gray20")
  }

  eta_plot
}

#' Print a \code{qa_fit} Result
#'
#' A formatted summary: the fitted model's scale and predictors, the size
#' of the donor/target samples, and a brief description of the resulting
#' adjustment's spread on the target sample.
#'
#' @param x A list returned by \code{\link{qa_fit}}.
#' @param ... Ignored; present for S3 method consistency.
#'
#' @export
print.qa_fit <- function(x, ...) {
  outcome_scale <- if (identical(class(x$model)[1], "glm")) "level" else "log"
  cat("<qa_fit result>\n")
  cat(sprintf("  Outcome scale:      %s\n", outcome_scale))
  cat(sprintf("  Predictors:         %s\n", paste(all.vars(stats::formula(x$model))[-1], collapse = ", ")))
  cat(sprintf("  Donor observations: %d\n", length(x$donor_resid)))
  cat(sprintf("  Target observations:%d\n", length(x$y_adjusted)))
  cat(sprintf("  y_adjusted range:   [%.4g, %.4g]\n",
              min(x$y_adjusted), max(x$y_adjusted)))
  cat("\nUse plot(x) for the eta(p) curve; see ?qa_fit for full component list.\n")
  invisible(x)
}
