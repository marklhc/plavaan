# Variance bounds for random starting values
# Reimplemented locally to avoid depending on unexported lavaan internals.
# Mirrors lavaan's bounded-estimation logic:
# - Lower bound: small fraction of observed variance (prevents near-zero)
# - Upper bound: multiple of observed variance (prevents explosive values)
variance_bounds <- function(obs_var) {
  lower <- max(1e-6, 0.2 * obs_var)
  upper <- 2 * obs_var
  c(lower, upper)
}

#' Generate Random Starting Values for Penalized Estimation
#'
#' Generates random starting vectors for use with [penalized_est()], perturbing
#' lavaan's default starting values by parameter type. The first row is always
#' the unperturbed base vector. Subsequent rows apply type-specific random
#' perturbations mirroring (but not identical to) lavaan's `rstarts` scheme.
#'
#' @param x A fitted lavaan model object.
#' @param n Integer. Number of starting vectors to generate (including the base).
#'   Default is 1, which returns only the base vector.
#'
#' @return A matrix with `n` rows (one per starting vector) and one column per
#'   free parameter. Column order matches `lavaan::coef(x)` /
#'   `lav_export_estimation(x)$starting_values`. The first row is always the
#'   unperturbed base vector.
#'
#' @details
#' Perturbation rules by parameter type (from `parTable`):
#'
#' \describe{
#'   \item{Factor loadings (`op == "=~"`)}: Left at base value. Lavaan itself
#'     does not randomize ordinary loadings in `rstarts`.
#'   \item{Regression coefficients (`op == "~"`, excluding intercepts)}: A random
#'     correlation \(r \sim U(-0.5, 0.5)\) is drawn and converted to the
#'     covariance scale via \(r \times \sqrt{\hat{\sigma}_{lhs}^2 \hat{\sigma}_{rhs}^2}\),
#'     using base residual variance estimates. This deviation from lavaan's
#'     `rstarts` (which fixes regression coefficients at 0) is intentional: in
#'     penalized estimation, zero-starts may collide with the penalty in
#'     unhelpful ways.
#'   \item{Covariances between exogenous variables (`op == "~~"`, `lhs != rhs`)}:
#'     Same treatment as regression coefficients (random correlation scaled by
#'     variances).
#'   \item{(Residual) variances (`op == "~~"`, `lhs == rhs`)}: Drawn uniformly
#'     between bounds computed from the observed variance of the corresponding
#'     variable. Lower bound is ~0.2 * observed variance (min 1e-6), upper bound
#'     is ~2 * observed variance, mirroring lavaan's bounded estimation.
#'   \item{Intercepts/means (`op == "~1"`)}: Left at base value. Lavaan does not
#'     randomize intercepts.
#' }
#'
#' Set a seed with `set.seed()` before calling for reproducibility.
#'
#' @keywords internal
#' @noRd
#' @importFrom stats runif setNames
random_start <- function(x, n = 1) {
  ff <- lavaan::lav_export_estimation(x)
  base <- ff$starting_values
  pt <- lavaan::parTable(x)

  # Filter to free parameters (same order as base vector)
  free_pt <- pt[pt$free > 0, ]
  n_params <- length(base)

  if (nrow(free_pt) != n_params) {
    stop(
      "Number of free parameters (",
      nrow(free_pt),
      ") doesn't match starting values length (",
      n_params,
      ")."
    )
  }

  # Get observed variances for bounds and correlation-to-covariance conversion
  sample_cov <- tryCatch(
    lavaan::lavInspect(x, "cov.ov"),
    error = function(e) lavaan::lavInspect(x, "sample.cov")
  )
  obs_variances <- setNames(diag(sample_cov), rownames(sample_cov))

  # Build a variance lookup by variable name. Include ALL variance entries
  # (both free and fixed) so that latent factor variances (e.g., from std.lv = TRUE)
  # are available for covariance perturbation.
  all_var_mask <- pt$op == "~~" & pt$lhs == pt$rhs
  var_lookup <- setNames(pt$est[all_var_mask], pt$lhs[all_var_mask])

  # For variable names present in both the parameter table and observed covariance,
  # update the variance estimate to the maximum of the two (base estimate and
  # observed variance) for safety against non-positive-definite matrices.
  for (vn in intersect(names(obs_variances), names(var_lookup))) {
    var_lookup[vn] <- max(var_lookup[vn], obs_variances[vn], na.rm = TRUE)
  }

  result <- t(replicate(n, base))

  if (n <= 1) {
    return(result)
  }

  for (i in seq_len(n_params)) {
    op <- free_pt$op[i]
    lhs <- free_pt$lhs[i]
    rhs <- free_pt$rhs[i]

    if (op == "=~") {
      # Factor loadings: leave at base value (no perturbation)
      next
    } else if (op == "~1") {
      # Intercepts/means: leave at base value (no perturbation)
      next
    } else if (op == "~") {
      # Regression coefficients: random correlation scaled by variances
      lhs_sd <- sqrt(abs(var_lookup[lhs]))
      rhs_sd <- sqrt(abs(var_lookup[rhs]))
      # Guard against NA (shouldn't happen, but safety first)
      if (!is.na(lhs_sd) && !is.na(rhs_sd)) {
        result[-1, i] <- base[i] + runif(n - 1, -0.5, 0.5) * lhs_sd * rhs_sd
      } else {
        # Fallback: just jitter by small fraction of base
        result[-1, i] <- base[i] * (1 + runif(n - 1, -0.2, 0.2))
      }
    } else if (op == "~~" & lhs != rhs) {
      # Covariances between exogenous vars: random correlation scaled by variances
      lhs_sd <- sqrt(abs(var_lookup[lhs]))
      rhs_sd <- sqrt(abs(var_lookup[rhs]))
      if (!is.na(lhs_sd) && !is.na(rhs_sd)) {
        result[-1, i] <- base[i] + runif(n - 1, -0.5, 0.5) * lhs_sd * rhs_sd
      } else {
        result[-1, i] <- base[i] * (1 + runif(n - 1, -0.2, 0.2))
      }
    } else if (op == "~~" & lhs == rhs) {
      # Variances: draw uniformly within bounds based on observed variance
      obs_var <- var_lookup[lhs]
      lb <- variance_bounds(obs_var)[1]
      ub <- variance_bounds(obs_var)[2]
      result[-1, i] <- runif(n - 1, lb, ub)
    }
  }

  result
}

#' Multistart Penalized Estimation
#'
#' Repeatedly calls [penalized_est()] with different starting vectors and returns
#' the solution with the lowest penalized objective value. This is recommended
#' when using non-convex penalties (`l0a`, `alf`), which can lead to local optima.
#'
#' @inheritParams penalized_est
#' @param n_starts Integer. Number of random starting vectors to try. The first
#'   start is always lavaan's default (unperturbed), so multistart is never worse
#'   than a single [penalized_est()] call. Default is 10.
#' @param starts Matrix or list of numeric vectors, each with length equal to the
#'   number of free parameters. If supplied, random generation is bypassed entirely
#'   and `n_starts` is ignored. A message is printed noting this.
#' @param keep_all Logical. If `TRUE`, attach all fitted lavaan objects as an
#'   attribute for inspection of alternative local solutions. Default is `FALSE`.
#' @param verbose Logical. If `TRUE`, print progress messages during optimization.
#'   Default is `FALSE`.
#'
#' @return A lavaan model object (same shape as [penalized_est()]'s return), with
#'   additional attributes:
#'
#' \describe{
#'   \item{`multistart`}{A data frame with one row per start, containing columns
#'     `start_id`, `objective` (final penalized objective value), and
#'     `converged` (logical). Rows are sorted by ascending objective.}
#'   \item{`all_fits`}{If `keep_all = TRUE`, a named list of all fitted lavaan
#'     objects, one per starting vector.}
#' }
#'
#' @details
#' Non-convex penalties like `l0a` and `alf` create rugged optimization surfaces
#' where the optimizer can settle in different local solutions depending on
#' starting values. Multistart optimization mitigates this risk by trying several
#' starts and selecting the best.
#'
#' Starting-value generation mirrors lavaan's `rstarts` scheme but with one key
#' deviation: regression coefficients are randomized (not fixed at 0 as in
#' lavaan), because zero-starts may collide with the penalty function in ways
#' that hinder convergence toward the global optimum. Perturbation rules by
#' parameter type follow the same logic as [lavaan::lavOptions]`rstarts`: factor
#' loadings and intercepts stay at base values, variances are drawn within bounds
#' based on observed variance, and regression/covariance parameters receive a
#' random correlation perturbation scaled to the covariance scale.
#'
#' Execution is sequential. For parallel execution, call [penalized_est()] with
#' `start = ...` directly and use `future.apply`, `parallel`, or similar.
#'
#' Users should set their own seed with `set.seed()` for reproducibility. This
#' function does not call `set.seed()` internally.
#'
#' @seealso [penalized_est]
#' @export
penalized_est_multistart <- function(
  x,
  w,
  pen_par_id = NULL,
  pen_diff_id = NULL,
  pen_fn = "l0a",
  pen_gr = NULL,
  se = "none",
  opt_control = list(),
  eps = .01,
  telescoping_control = list(
    eps_1 = 1,
    eps_end = 1e-5,
    eps_steps = 20,
    warm_start = FALSE
  ),
  n_starts = 10,
  starts = NULL,
  keep_all = FALSE,
  verbose = FALSE
) {
  # Generate or validate starting values
  if (is.null(starts)) {
    all_starts <- random_start(x, n = n_starts)
  } else {
    # Convert list to matrix if needed
    if (is.list(starts)) {
      starts <- do.call(rbind, lapply(starts, as.numeric))
    }
    if (!is.matrix(starts)) {
      starts <- matrix(starts, nrow = 1)
    }
    message(
      "Custom starting values provided (",
      nrow(starts),
      " starts). Ignoring n_starts."
    )
    all_starts <- starts
  }

  n_total <- nrow(all_starts)
  results <- vector("list", n_total)
  records <- vector("list", n_total)

  for (i in seq_len(n_total)) {
    if (verbose) {
      message(sprintf("Start %d / %d...", i, n_total))
    }

    fit <- tryCatch(
      penalized_est(
        x = x,
        w = w,
        pen_par_id = pen_par_id,
        pen_diff_id = pen_diff_id,
        pen_fn = pen_fn,
        pen_gr = pen_gr,
        se = se,
        opt_control = opt_control,
        eps = eps,
        telescoping_control = telescoping_control,
        start = all_starts[i, ]
      ),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      records[[i]] <- data.frame(
        start_id = i,
        objective = NA_real_,
        converged = FALSE,
        stringsAsFactors = FALSE
      )
      results[[i]] <- NULL
      if (verbose) {
        message(sprintf("  Start %d: failed", i))
      }
    } else {
      obj_val <- fit@optim$fx
      conv <- fit@optim$converged
      records[[i]] <- data.frame(
        start_id = i,
        objective = obj_val,
        converged = conv,
        stringsAsFactors = FALSE
      )
      results[[i]] <- fit
      if (verbose) {
        status <- if (conv) "converged" else "did not converge"
        message(sprintf("  Start %d: %s, objective = %.6f", i, status, obj_val))
      }
    }
  }

  # Combine records and sort by objective (NA last via default behavior)
  ms_table <- do.call(rbind, records)

  # Find best: among converged runs, pick the lowest objective.
  # Table is sorted by ascending objective afterward.
  converged_idx <- which(ms_table$converged)
  if (length(converged_idx) > 0) {
    best_conv <- converged_idx[which.min(ms_table$objective[converged_idx])]
  } else {
    # No convergence -- warn and pick lowest objective among all runs
    warning(
      "None of the ",
      n_total,
      " optimization runs converged. Returning the run with the lowest objective."
    )
    valid_idx <- which(!is.na(ms_table$objective))
    if (length(valid_idx) == 0) {
      stop("All starts failed. No fit to return.")
    }
    best_conv <- valid_idx[which.min(ms_table$objective[valid_idx])]
  }

  best_start_id <- ms_table$start_id[best_conv]
  ms_table <- ms_table[order(ms_table$objective), ]
  best_fit <- results[[best_start_id]]

  # Attach multistart summary table
  attr(best_fit, "multistart") <- ms_table

  # Optionally attach all fits
  if (keep_all) {
    named_fits <- setNames(results, paste0("start_", seq_len(n_total)))
    attr(best_fit, "all_fits") <- named_fits
  }

  best_fit
}
