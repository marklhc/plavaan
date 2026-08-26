# Penalized objective function
penalized_obj <- function(x, obj_fn, w, pen_fn, pen_par_id, diff_configs, ...) {
  out <- obj_fn(x)

  if (!is.null(pen_par_id)) {
    out <- out + w * sum(pen_fn(x[pen_par_id], ...))
  }

  if (!is.null(diff_configs)) {
    pen_diff <- lapply(diff_configs, function(cfg) {
      x_mat <- matrix(x[cfg$mat], nrow = nrow(cfg$mat), ncol = ncol(cfg$mat))
      # A value-driven NaN here (e.g. log() of a non-positive loading
      # during optimization) is expected and handled via na.rm = TRUE
      # below, so suppress R's low-level "NaNs produced" warning.
      x_trans <- suppressWarnings(as.matrix(cfg$trans(x_mat)))

      # Use pre-computed combn_idx and rescale_val
      diffs <- x_trans[cfg$combn_idx[1, ], , drop = FALSE] -
        x_trans[cfg$combn_idx[2, ], , drop = FALSE]

      sum(pen_fn(diffs, ...), na.rm = TRUE) * cfg$rescale_val
    })
    out <- out + w * sum(unlist(pen_diff))
  }
  out
}

#' Penalized Parameter Estimation for Longitudinal CFA Models
#'
#' Performs penalized estimation on a lavaan model object by optimizing a
#' penalized objective function. The function extracts the objective function
#' from a lavaan model, applies a penalty function to specified parameters
#' or pairwise differences of parameters, and returns an updated model with
#' the optimized parameter estimates.
#'
#' @param x A fitted lavaan model object from which estimation components will
#'   be extracted.
#' @param w Numeric scalar. Penalty weight (multiplier) applied to the penalty
#'   terms.
#' @param pen_par_id Integer vector of parameter IDs to apply the penalty function
#'   directly to, in the same order as returned by `lavaan::coef()` and by
#'   [lavaan::parTable()], with only the free elements.
#' @param pen_diff_id A named list of integer matrices of free-parameter IDs
#'   (same order as `lavaan::coef()` / the `free` column of [lavaan::parTable()]).
#'   Each matrix has one row per group or time point and one column per matched
#'   parameter; the penalty is the sum of pairwise row differences within each
#'   column, rescaled by `(nrow - 1) / ncombn(nrow, 2)`. Structural `NA` entries
#'   mark a parameter absent in that row and are excluded from the differences.
#' @param pen_fn A character string (`"l0a"` or `"alf"`) or a function that computes
#'   the penalty. Default is `"l0a"`.
#' @param pen_gr A function that computes the gradient of the penalty function.
#'   If `pen_fn` is `"l0a"` or `"alf"`, this is automatically set.
#' @param se Character string specifying the type of standard errors to compute.
#'   Options are `"none"` (default; no standard errors) or `"robust.huber.white"`
#'   (robust sandwich estimator using numerical Hessian and first-order information,
#'   which is the same as used in the `"mlr"` estimator).
#' @param test Character string specifying the model test used by the fit
#'   evaluation on the returned object (`fitmeasures()` and the chi-square test
#'   in `summary()`), via an internal "frozen" refit. Fit evaluation for
#'   penalized fits is **experimental**, so it is disabled by default:
#'   `"none"` (default) means no model test is run, `fitmeasures()` is
#'   unavailable, and `summary()` shows no chi-square test. Set to `"Chisq"`
#'   (ML/PML estimators) or `"SatorraBentler"` (WLSMV/MLM/MLR) to enable fit
#'   measures and the chi-square test; an experimental notice is then shown
#'   when `fitmeasures()` or `summary()` is called.
#' @param opt_control A list of control parameters passed to [stats::nlminb()].
#'   Default includes `eval.max = 2e4`, `iter.max = 1e4`, and `abs.tol = 1e-20`.
#' @param start Numeric vector of starting values for the optimizer, or `NULL`
#'   (default) to use lavaan's default starting values. If supplied, its length
#'   must match the number of free parameters in the model.
#' @param eps A positive numeric scalar used by the built-in penalties, or
#'   `"telescoping"` to fit a sequence of decreasing epsilon values. Default is
#'   `.01`. This argument does not alter custom `pen_fn` or `pen_gr` functions.
#' @param telescoping_control A named list controlling telescoping, with
#'   `eps_1` (default `1`), `eps_end` (default `1e-5`), `eps_steps` (default
#'   `20`), and `warm_start` (default `FALSE`). When `warm_start` is `FALSE`,
#'   every epsilon stage uses the original starting values; when `TRUE`, each
#'   stage after the first uses the preceding stage's estimates.
#' @param ... Additional arguments passed to a user-supplied `pen_fn` /
#'   `pen_gr`. Custom penalty functions must accept `...`. Built-in penalties
#'   (`"l0a"`, `"alf"`) ignore it.
#'
#' @section Warning:
#' The returned object is not fitted using standard ML. Standard errors reported
#' by `summary()` or `parameterEstimates()` will be missing unless
#' `se = "robust.huber.white"` was specified. Even then, they are based on an
#' experimental sandwich approximation and should be interpreted with caution.
#'
#' Fit evaluation (`fitmeasures()` and the chi-square test in `summary()`) is
#' also **experimental** and disabled by default (`test = "none"`). Enable it
#' with `test = "Chisq"` (or `"SatorraBentler"`); interpret any resulting fit
#' indices with caution, as they are based on a frozen refit at the penalized
#' estimates with the effective degrees of freedom.
#'
#' @return A lavaan model object updated with the penalized parameter estimates.
#'   The object has S4 class `plavaan` (a subclass of `lavaan`) and a
#'   `penalized` attribute recording the penalty specification, which enables
#'   [effective_df()] and, when `test` is not `"none"`, `fitmeasures()` and
#'   `summary()` with effective degrees of freedom. With `eps = "telescoping"`,
#'   it also includes a `"telescoping"` data frame with per-stage epsilon
#'   values, parameter changes, objective values, and convergence indicators.
#'
#' @details
#' The function uses `nlminb()` to minimize a penalized objective function that
#' combines the standard lavaan objective function with a penalty term. Only the
#' parameter estimates and the log-likelihood should be interpreted. The
#' returned object was not "fitted" (`do.fit = FALSE`) to avoid users
#' interpreting the standard errors, which are generally not valid with
#' penalized estimation. The nominal model degrees of freedom can also be
#' misleading, as the penalized model is often under-identified;
#' [effective_df()] reports the effective number of parameters and the
#' effective model degrees of freedom. When `test` is not `"none"`,
#' `fitmeasures()` / `summary()` on the returned object additionally report
#' fit indices at the effective df (frozen refit at the penalized estimates);
#' this fit evaluation is experimental and disabled by default. If the
#' optimization does not converge (convergence code != 0), a warning is issued.
#'
#' With `eps = "telescoping"`, the model is fit along a log-spaced sequence from
#' `telescoping_control$eps_1` to `telescoping_control$eps_end`. By default,
#' each stage uses the original starting values; set
#' `telescoping_control$warm_start = TRUE` to initialize later stages from the
#' preceding solution. The sequence stops when the largest absolute change
#' between consecutive parameter vectors is at most `5e-4`. The returned object
#' has a `"telescoping"` attribute with stage diagnostics.
#'
#' @seealso \code{\link[lavaan]{lavaan}}, \code{\link[stats]{nlminb}}
#'
#' @importFrom stats nlminb
#' @examples
#' library(lavaan)
#'
#' # Define a longitudinal factor model with PoliticalDemocracy data
#' model <- "
#'   dem60 =~ y1 + y2 + y3 + y4
#'   dem65 =~ y5 + y6 + y7 + y8
#'   dem60 ~~ dem65
#'   dem60 ~~ 1 * dem60
#'   dem65 ~~ NA * dem65
#'   dem60 ~ 0
#'   dem65 ~ NA * 1
#'   y1 ~~ y5
#'   y2 ~~ y6
#'   y3 ~~ y7
#'   y4 ~~ y8
#' "
#'
#' # Fit the model without constraints first to get parameter table
#' fit_un <- cfa(model, data = PoliticalDemocracy, std.lv = TRUE,
#'               meanstructure = TRUE, do.fit = FALSE)
#'
#' # Get parameter IDs
#' pt <- parTable(fit_un)
#' # Loadings
#' load_60 <- pt$free[pt$op == "=~" & pt$lhs == "dem60"]
#' load_65 <- pt$free[pt$op == "=~" & pt$lhs == "dem65"]
#' # Intercepts
#' int_60 <- pt$free[pt$op == "~1" & pt$lhs %in% c("y1", "y2", "y3", "y4")]
#' int_65 <- pt$free[pt$op == "~1" & pt$lhs %in% c("y5", "y6", "y7", "y8")]
#'
#' # Apply penalized estimation to penalize differences in loadings and intercepts
#' pen_fit <- penalized_est(
#'     x = fit_un,
#'     w = 0.03,
#'     pen_diff_id = list(
#'         loadings = rbind(load_60, load_65),
#'         intercepts = rbind(int_60, int_65)
#'     ),
#'     pen_fn = "l0a"
#' )
#'
#' # Compare parameter estimates
#' summary(pen_fit)
#'
#' # Effective number of parameters and degrees of freedom
#' effective_df(pen_fit)
#'
#' @importFrom utils modifyList
#' @export
penalized_est <- function(
  x,
  w,
  pen_par_id = NULL,
  pen_diff_id = NULL,
  pen_fn = "l0a",
  pen_gr = NULL,
  se = "none",
  test = "none",
  opt_control = list(),
  start = NULL,
  eps = .01,
  telescoping_control = list(
    eps_1 = 1,
    eps_end = 1e-5,
    eps_steps = 20,
    warm_start = FALSE
  ),
  ...
) {
  if (is.numeric(eps) && length(eps) == 1 && is.finite(eps) && eps > 0) {
    eps_seq <- eps
  } else if (identical(eps, "telescoping")) {
    telescoping_control <- modifyList(
      list(eps_1 = 1, eps_end = 1e-5, eps_steps = 20, warm_start = FALSE),
      telescoping_control
    )
    eps_1 <- telescoping_control$eps_1
    eps_end <- telescoping_control$eps_end
    eps_steps <- telescoping_control$eps_steps
    warm_start <- telescoping_control$warm_start
    if (
      !is.numeric(eps_1) ||
        length(eps_1) != 1 ||
        !is.finite(eps_1) ||
        eps_1 <= 0 ||
        !is.numeric(eps_end) ||
        length(eps_end) != 1 ||
        !is.finite(eps_end) ||
        eps_end <= 0 ||
        !is.numeric(eps_steps) ||
        length(eps_steps) != 1 ||
        !is.finite(eps_steps) ||
        eps_steps < 1 ||
        eps_steps != as.integer(eps_steps) ||
        eps_1 < eps_end ||
        !is.logical(warm_start) ||
        length(warm_start) != 1 ||
        is.na(warm_start)
    ) {
      stop(
        "telescoping_control must contain positive eps_1 >= eps_end, an ",
        "integer eps_steps >= 1, and a non-missing logical warm_start."
      )
    }
    eps_seq <- exp(seq(log(eps_1), log(eps_end), length.out = eps_steps))
  } else {
    stop("eps must be a positive numeric scalar or 'telescoping'.")
  }

  if (!se %in% c("none", "robust.huber.white")) {
    warning(
      "se must be either 'none' or 'robust.huber.white'. ",
      "Defaulting to 'none'"
    )
    se <- "none"
  }

  if (!is.character(test) || length(test) != 1 || !nzchar(test)) {
    stop(
      "test must be a character string, e.g. 'none', 'Chisq', or ",
      "'SatorraBentler'."
    )
  }

  pen_fn_name <- if (is.character(pen_fn) && length(pen_fn) == 1) {
    pen_fn
  } else {
    NULL
  }
  if (
    !is.function(pen_fn) &&
      !identical(pen_fn_name, "l0a") &&
      !identical(pen_fn_name, "alf")
  ) {
    stop("pen_fn must be 'l0a', 'alf', or a function.")
  }
  if (identical(pen_fn_name, "l0a")) {
    if (!is.null(pen_gr)) {
      message(
        "pen_gr is ignored when pen_fn is 'l0a'; using the built-in gradient function."
      )
    }
    pen_gr <- gr_l0a
  } else if (identical(pen_fn_name, "alf")) {
    if (!is.null(pen_gr)) {
      message(
        "pen_gr is ignored when pen_fn is 'alf'; using the built-in gradient function."
      )
    }
    pen_gr <- gr_alf
  }

  fit_stage <- function(stage_eps, stage_start) {
    pen_fn_stage <- pen_fn
    pen_gr_stage <- pen_gr
    if (identical(pen_fn_name, "l0a")) {
      pen_fn_stage <- function(z, ...) l0a(z, eps = stage_eps)
      pen_gr_stage <- function(z, ...) gr_l0a(z, eps = stage_eps)
    } else if (identical(pen_fn_name, "alf")) {
      pen_fn_stage <- function(z, ...) alf(z, eps = stage_eps)
      pen_gr_stage <- function(z, ...) gr_alf(z, eps = stage_eps)
    }
    # Record the penalty specification so that downstream utilities
    # (effective_df(), fitmeasures(), summary()) can resolve it from the
    # returned object.
    pen_spec <- list(
      w = w,
      pen_par_id = pen_par_id,
      pen_diff_id = pen_diff_id,
      pen_fn = if (is.null(pen_fn_name)) pen_fn else pen_fn_name,
      eps = stage_eps,
      test = test
    )
    penalized_est_stage(
      x = x,
      w = w,
      pen_par_id = pen_par_id,
      pen_diff_id = pen_diff_id,
      pen_fn = pen_fn_stage,
      pen_gr = pen_gr_stage,
      se = se,
      opt_control = opt_control,
      start = stage_start,
      pen_spec = pen_spec,
      ...
    )
  }

  out <- NULL
  stage_count <- 0L
  eps_used <- numeric(length(eps_seq))
  par_changes <- numeric(length(eps_seq))
  objectives <- numeric(length(eps_seq))
  converged <- logical(length(eps_seq))
  original_start <- start
  stage_start <- original_start
  for (stage_eps in eps_seq) {
    stage_count <- stage_count + 1L
    fit <- fit_stage(stage_eps, stage_start)
    if (any(!is.finite(fit@optim$x))) {
      stop(
        "Optimization produced non-finite estimates at eps = ",
        stage_eps,
        "."
      )
    }
    change <- if (is.null(out)) {
      NA_real_
    } else {
      max(abs(fit@optim$x - out@optim$x))
    }
    eps_used[stage_count] <- stage_eps
    par_changes[stage_count] <- change
    objectives[stage_count] <- fit@optim$fx
    converged[stage_count] <- fit@optim$converged
    out <- fit
    if (!is.na(change) && change <= 5e-4) {
      break
    }
    if (identical(eps, "telescoping") && warm_start) {
      stage_start <- fit@optim$x
    } else {
      stage_start <- original_start
    }
  }
  if (identical(eps, "telescoping")) {
    attr(out, "telescoping") <- data.frame(
      eps = eps_used[seq_len(stage_count)],
      max_abs_change = par_changes[seq_len(stage_count)],
      objective = objectives[seq_len(stage_count)],
      converged = converged[seq_len(stage_count)]
    )
  }
  out
}

resolve_penalty_functions <- function(pen_fn, pen_gr = NULL) {
  if (!is.function(pen_fn)) {
    stop("pen_fn must be a function.")
  }
  if (!is.null(pen_gr) && !is.function(pen_gr)) {
    stop("pen_gr must be a function or NULL.")
  }
  list(pen_fn = pen_fn, pen_gr = pen_gr)
}

make_diff_configs <- function(pen_diff_id) {
  if (is.null(pen_diff_id)) {
    return(NULL)
  }
  lapply(pen_diff_id, function(mat) {
    nrow_x <- nrow(mat)
    if (nrow_x < 2) {
      combn_idx <- matrix(integer(), nrow = 2, ncol = 0)
      rescale_val <- 0
    } else {
      combn_idx <- combn(nrow_x, 2)
      rescale_val <- (nrow_x - 1) / ncol(combn_idx)
    }
    grad_idx <- lapply(seq_len(nrow_x), function(i) {
      list(
        idx1 = which(combn_idx[1, ] == i),
        idx2 = which(combn_idx[2, ] == i)
      )
    })
    list(
      mat = mat,
      trans = identity,
      gr_trans = function(x) rep(1, length(x)),
      combn_idx = combn_idx,
      rescale_val = rescale_val,
      grad_idx = grad_idx
    )
  })
}

make_penalized_fit <- function(x, opt, pen_spec = NULL) {
  x_opt <- x@Options
  x_opt$start <- opt$par
  x_opt$do.fit <- FALSE
  x_opt$se <- "none"
  out <- lavaan::lavaan(
    lavaan::parTable(x),
    slotOptions = x_opt,
    slotSampleStats = x@SampleStats,
    slotData = x@Data
  )
  # Note: S4 slot assignment is copy-on-write, so the nlminb information is
  # carried by add_nlminb_info()'s return value (not in-place).
  out <- add_nlminb_info(out, opt)
  if (!is.null(pen_spec)) {
    attr(out, "penalized") <- pen_spec
    # A shared environment: R copies objects on write, so an attribute set
    # inside a later method call would not persist on the caller's object.
    # An environment is a shared reference, so every copy of this fit can
    # read and write the same cache (see plavaan_frozen()).
    attr(out, "plavaan.cache") <- new.env(size = 0)
    class(out) <- "plavaan"
  }
  out
}

penalized_est_stage <- function(
  x,
  w,
  pen_par_id,
  pen_diff_id,
  pen_fn,
  pen_gr,
  se,
  opt_control,
  start,
  pen_spec = NULL,
  ...
) {
  # Define default control parameters
  control_defaults <- list(
    eval.max = 2e4,
    iter.max = 1e4,
    abs.tol = 1e-20
  )

  # Merge with user input
  control <- modifyList(control_defaults, opt_control)

  ff <- lavaan::lav_export_estimation(x)
  if (is.null(start)) {
    start <- ff$starting_values
  } else if (length(start) != length(ff$starting_values)) {
    stop(
      "start must have length ",
      length(ff$starting_values),
      " (number of free parameters), but has length ",
      length(start),
      "."
    )
  }
  penalty <- resolve_penalty_functions(pen_fn, pen_gr)
  pen_fn <- penalty[[1]]
  pen_gr <- penalty[[2]]
  diff_configs <- make_diff_configs(pen_diff_id)
  objective_fn <- function(pars) {
    ff$objective_function(pars, lavaan_model = x)
  }
  gradient_fn <- if (!is.null(pen_gr)) {
    function(pars) ff$gradient_function(pars, lavaan_model = x)
  }
  f1 <- function(v) {
    penalized_obj(
      v,
      obj_fn = objective_fn,
      w = w,
      pen_fn = pen_fn,
      pen_par_id = pen_par_id,
      diff_configs = diff_configs,
      ...
    )
  }
  gr1 <- if (!is.null(pen_gr)) {
    function(v) {
      penalized_gr(
        v,
        gr_fn = gradient_fn,
        w = w,
        pen_gr = pen_gr,
        pen_par_id = pen_par_id,
        diff_configs = diff_configs,
        ...
      )
    }
  } else {
    NULL # Let nlminb compute numerical gradient
  }
  opt <- nlminb(
    start,
    objective = f1,
    gradient = gr1,
    control = control
  )
  if (opt$convergence != 0) {
    warning(
      "Optimization did not converge. Try using better starting values, ",
      "or adjusting optimization control parameters."
    )
  }
  out <- make_penalized_fit(x, opt, pen_spec = pen_spec)
  if (se == "robust.huber.white") {
    hess <- numDeriv::hessian(f1, opt$par)
    attr(out, "hessian") <- hess
    out <- try(add_vcov_pen(out, hess), silent = TRUE)
    if (inherits(out, "try-error")) {
      warning(
        "Computation of robust sandwich estimator standard errors failed ",
        "(likely due to a singular or nearly-singular Hessian). ",
        "Standard errors are not available."
      )
      out <- make_penalized_fit(x, opt, pen_spec = pen_spec)
    }
  }
  out
}

add_nlminb_info <- function(fit, opt) {
  fit@optim$x <- opt$par
  fit@optim$fx <- opt$objective
  fit@optim$iterations <- opt$iterations
  fit@optim$converged <- opt$convergence == 0
  fit@optim$control <- opt$control
  fit@optim$dx <- opt$gradient
  fit@optim$npar <- length(opt$par)
  fit
}

#' @importFrom lavaan lavInspect
add_vcov_pen <- function(fit, hess) {
  meat <- lavInspect(fit, "information.first.order")
  H_inv <- solve(hess)
  vc_out <- H_inv %*% meat %*% H_inv

  fit@vcov$se <- "robust.huber.white"
  fit@vcov$vcov <- vc_out / lavInspect(fit, "nobs")
  fit@vcov$information <- "observed"
  fit@Options$se <- "robust.huber.white"
  fit@Options$information <- rep("observed", 2)
  fit@ParTable$se <- 0 * fit@ParTable$est
  fit@ParTable$se[which(fit@ParTable$free > 0)] <- sqrt(diag(
    fit@vcov$vcov
  ))
  fit
}

penalized_gr <- function(x, gr_fn, w, pen_gr, pen_par_id, diff_configs, ...) {
  out <- gr_fn(x)

  if (!is.null(pen_par_id)) {
    out <- out + w * hot_gr(x, pen_par_id, pen_gr, ...)
  }

  if (!is.null(diff_configs)) {
    pen_diff_gr <- lapply(diff_configs, function(cfg) {
      # As in penalized_obj(): a value-driven NaN from cfg$trans (e.g.
      # log() of a non-positive loading) is expected here and is
      # explicitly handled below, so suppress R's warning.
      x_mat <- suppressWarnings(
        as.matrix(cfg$trans(matrix(x[cfg$mat], nrow = nrow(cfg$mat))))
      )

      diffs <- x_mat[cfg$combn_idx[1, ], , drop = FALSE] -
        x_mat[cfg$combn_idx[2, ], , drop = FALSE]

      grad_contribs <- pen_gr(diffs, ...)
      grad <- matrix(0, nrow = nrow(x_mat), ncol = ncol(x_mat))

      # Loop is now incredibly fast because indices are pre-calculated
      for (i in seq_len(nrow(x_mat))) {
        g1 <- grad_contribs[cfg$grad_idx[[i]]$idx1, , drop = FALSE]
        g2 <- grad_contribs[cfg$grad_idx[[i]]$idx2, , drop = FALSE]

        grad[i, ] <- colSums(g1, na.rm = TRUE) - colSums(g2, na.rm = TRUE)
      }

      # Only mark cells as NA where cfg$mat itself is structurally NA
      # (e.g. an indicator missing from one group/time point). Do NOT
      # use is.na(x_mat) for this: is.na() also matches NaN, which
      # would wrongly overwrite value-driven NaNs (e.g. from log() of a
      # non-positive loading) with plain NA and hide them from the
      # check below.
      grad[is.na(cfg$mat)] <- NA
      grad_vec <- as.vector(grad) * cfg$rescale_val * cfg$gr_trans(x[cfg$mat])

      # Re-map back to the full parameter vector space (hot_gr equivalent inline)
      # `cfg$mat` may contain *structural* NAs (e.g. an indicator missing
      # from one group/time point) -- those positions are dropped here.
      # Separately, `grad_vec` can contain *value-driven* NaNs at
      # positions where `cfg$mat` is valid, e.g. when `cfg$trans` is
      # `log` and the current parameter estimate is negative or zero.
      # These two cases must be tracked separately: silently using
      # `na.omit()` on both and assuming they line up (as before) breaks
      # as soon as a value-driven NaN appears among otherwise valid
      # indices, misaligning the replacement vector.
      idx <- as.numeric(cfg$mat)
      structural_valid <- !is.na(idx)
      value_nan <- is.nan(grad_vec[structural_valid])
      if (any(value_nan)) {
        warning(
          "Gradient of some penalty is undefined ",
          "(NaN) for ",
          sum(value_nan),
          " parameter(s). This could happen when ",
          "for example, log() was applied to a non-positive loading estimate. ",
          "These contributions are set to 0; consider using a different ",
          "transformation function.",
          call. = FALSE
        )
      }
      grad_vec_valid <- grad_vec[structural_valid]
      grad_vec_valid[value_nan] <- 0
      full_grad <- 0 * x
      full_grad[idx[structural_valid]] <- grad_vec_valid
      full_grad
    })
    out <- out + w * Reduce(`+`, pen_diff_gr)
  }
  out
}
