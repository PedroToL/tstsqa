# tstsqa

Bias reduction for two-sample two-stage (TSTS) survey-to-survey imputation
via quantile adjustment, implementing the method described in Torres-Lopez,
*"Bias Reduction in TSTSLS Applications: A Quantile Adjustment Method."*

## What the method does

Suppose you want to study the relationship between income (`y`) and some
other variable (`z`) — consumption, wealth, parental income — but no single
survey observes both. You have a **donor** sample with `y` and some
predictors `X`, and a **target** sample with `z` and the same `X`, but never
`y` and `z` together. The standard fix is to fit a model for `y` on `X` in
the donor sample and use it to impute `y` in the target sample.

That imputation is biased in two ways:

- **Variance bias**: predicted values are a conditional mean, so they're
  mechanically less spread out than the real `y` — the imputed distribution
  is too compressed.
- **Covariance bias**: whatever part of `y` isn't explained by `X` is lost
  entirely, so the imputed values understate (or distort) their true
  relationship with `z`.

Simply rescaling the variance (stochastic augmentation) fixes the first
problem but not the second. **Quantile adjustment** fixes both at once: it
shifts each imputed value by the gap between the *observed* and *predicted*
quantile functions at that value's rank, `eta(p) = Q_y(p) - Q_yhat(p)`. This
exactly restores the full donor distribution, and — because the correction
tracks each observation's position in the predicted-income ranking — it can
also shift the imputed covariance with `z` back toward the truth.

Whether it helps or hurts depends on how well the chosen predictors align
with both `y` and `z`; `qa_diagnose()` (below) is built to help decide that
before you run the adjustment.

## Installation

```r
# install.packages("devtools")
devtools::install_github("PedroToL/tstsqa")
```

## Toy example

### 1. Fit the adjustment with `qa_fit()`

```r
library(tstsqa)

set.seed(123)
n <- 5000
X <- rnorm(n)
y <- exp(1 + 0.8 * X + rnorm(n, sd = 1.5))   # substantial unexplained variance

donor  <- data.frame(y = y[1:2500], X = X[1:2500])
target <- data.frame(X = X[2501:5000])

result <- qa_fit(donor, target, y_var = "y", x_vars = "X", outcome_scale = "log",
                  plotting = TRUE, annotate_p = 0.9)

result$eta_plot
```

<img src="man/figures/README-eta-plot.png" width="600"/>

`result$y_adjusted` holds the corrected imputed values for the target
sample; `result$eta_plot` shows the estimated quantile-gap function
`eta(p)` that produced them — the raw grid estimates as points, the
smoothed curve used internally, and (since `annotate_p = 0.9`) the exact
adjustment applied at the 90th percentile.

### 2. Choose predictors and check the regime with `qa_diagnose()`

Adding a predictor doesn't always help: one that predicts `z` well but adds
little to predicting `y` can inflate the correction beyond what's actually
justified. `qa_diagnose()` screens candidate predictors against this
failure mode via bootstrap, then searches the survivors for the
best-fitting specification:

```r
set.seed(7)
n2 <- 10000
X1 <- rnorm(n2); X2 <- rnorm(n2); X3 <- rnorm(n2)

y2 <- exp(1 + 0.6 * X1 + 0.4 * X3 + rnorm(n2, sd = 0.8))   # X2 barely predicts y
z1 <- 0.3 * X1 + 0.9 * X2 + rnorm(n2, sd = 0.5)            # ...but strongly predicts z1
z2 <- 0.5 * X1 + 0.5 * X3 + rnorm(n2, sd = 0.5)

donor2  <- data.frame(y = y2[1:5000], X1 = X1[1:5000], X2 = X2[1:5000], X3 = X3[1:5000])
target2 <- data.frame(X1 = X1[5001:10000], X2 = X2[5001:10000], X3 = X3[5001:10000],
                       z1 = z1[5001:10000], z2 = z2[5001:10000])

diagnosis <- qa_diagnose(donor2, target2, y_var = "y", z_vars = c("z1", "z2"),
                          x_vars = c("X1", "X2", "X3"), plotting = TRUE)

diagnosis$S_plot
```

<img src="man/figures/README-s-plot.png" width="600"/>

`X2` gets screened out (it inflates the correction for `z1` far more than
it improves the fit for `y`), leaving `X1` and `X3` — both comfortably
under the screening thresholds shown as dashed lines. `diagnosis$rho_star`
then reports, for each `z`, the residual correlation at which the
adjustment would exactly recover the true covariance — a diagnostic for
whether you're likely under-correcting, over-correcting, or moving in the
wrong direction entirely.

### 3. Refit explicitly on the selected predictors

`qa_diagnose()` tells you *which* predictors to use; running `qa_fit()`
again explicitly on `diagnosis$selected_predictors` gives you both a fresh
`eta_plot` for that final specification and a fitted model. The donor
sample is the only place we can actually check the first stage's
predictions against reality — `y` is never observed in the target sample,
in a real application or in this simulation:

```r
final_fit <- qa_fit(donor2, target2, y_var = "y",
                     x_vars = diagnosis$selected_predictors,
                     outcome_scale = "log", plotting = TRUE, annotate_p = 0.9)

final_fit$eta_plot
```

<img src="man/figures/README-eta-plot-final.png" width="600"/>

```r
observed_log_y  <- log(donor2$y)
predicted_log_y <- observed_log_y - final_fit$donor_resid   # donor_resid = y_model - y_hat_d

density_df <- data.frame(
  log_y = c(observed_log_y, predicted_log_y),
  type  = rep(c("Observed", "Predicted"), each = length(observed_log_y))
)

ggplot2::ggplot(density_df, ggplot2::aes(log_y, color = type, linetype = type)) +
  ggplot2::geom_density(linewidth = 1.1) +
  ggplot2::scale_linetype_manual(values = c(Observed = "solid", Predicted = "dotted"))
```

<img src="man/figures/README-density-plot.png" width="600"/>

The predicted distribution is visibly more peaked and thinner-tailed than
the observed one — the variance bias the whole method is built to correct.
`qa_fit()`'s quantile adjustment (Section 1) is exactly what fixes this:
`result$y_adjusted` restores the full observed shape rather than just
rescaling the variance.

See `?qa_fit` and `?qa_diagnose` for full argument documentation.

## License

MIT © Pedro Torres-Lopez
