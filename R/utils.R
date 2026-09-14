# ============================================================================
# Shared internal helpers used by both qa_fit() and qa_diagnose().
# ============================================================================

# --- Check that donor and target agree on categorical levels for x_vars ---
# For every x_var that is a factor or character column, compares the set of
# levels present in donor_data vs target_data:
#   - a level present in target but NOT donor is a hard problem: a model
#     trained on donor can never have a coefficient for it, so
#     predict(model, newdata = target_data) is guaranteed to fail with a
#     "has new levels" error the moment that row is reached. This stops
#     immediately, before any (possibly expensive) fitting happens, rather
#     than letting the failure surface deep inside a bootstrap loop.
#   - a level present in donor but NOT target is not immediately fatal --
#     target simply never needs a prediction for it -- but is still worth
#     a warning, since it often indicates a coding mismatch worth a second
#     look rather than a genuine feature of the population.
# Purely numeric x_vars are skipped entirely (nothing to compare).
.check_factor_levels <- function(donor_data, target_data, x_vars) {

  extra_in_target <- list()  # x_var -> levels in target but not donor (fatal)
  extra_in_donor  <- list()  # x_var -> levels in donor but not target (warn)

  for (v in x_vars) {
    donor_col  <- donor_data[[v]]
    target_col <- target_data[[v]]

    donor_is_cat  <- is.factor(donor_col)  || is.character(donor_col)
    target_is_cat <- is.factor(target_col) || is.character(target_col)
    if (!donor_is_cat && !target_is_cat) next  # numeric predictor, nothing to check

    donor_levels  <- if (is.factor(donor_col))  levels(donor_col)  else unique(as.character(donor_col))
    target_levels <- if (is.factor(target_col)) levels(target_col) else unique(as.character(target_col))

    missing_from_donor <- setdiff(target_levels, donor_levels)
    missing_from_target <- setdiff(donor_levels, target_levels)

    if (length(missing_from_donor) > 0) extra_in_target[[v]] <- missing_from_donor
    if (length(missing_from_target) > 0) extra_in_donor[[v]] <- missing_from_target
  }

  if (length(extra_in_target) > 0) {
    detail <- paste(sprintf("'%s' (target has: %s)", names(extra_in_target),
                             vapply(extra_in_target, paste, character(1), collapse = ", ")),
                     collapse = "; ")
    stop(sprintf(
      "target_data has categorical level(s) not present in donor_data, so a donor-trained model can never predict them: %s. ",
      detail
    ), "Recode, drop, or combine these levels before calling this function.")
  }

  if (length(extra_in_donor) > 0) {
    detail <- paste(sprintf("'%s' (donor has: %s)", names(extra_in_donor),
                             vapply(extra_in_donor, paste, character(1), collapse = ", ")),
                     collapse = "; ")
    warning(sprintf(
      "donor_data has categorical level(s) not present in target_data: %s. ", detail
    ), "This is not necessarily a problem (target simply never needs these), ",
    "but confirm this reflects a real feature of the population rather than a coding mismatch.")
  }

  invisible(NULL)
}
