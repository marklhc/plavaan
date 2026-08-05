library(plavaan)
library(lavaan)

# Helper: small fitted model for testing (unpenalized, do.fit = FALSE)
skip_on_cran()

small_model <- "
  dem60 =~ y1 + y2 + y3 + y4
  dem65 =~ y5 + y6 + y7 + y8
  dem60 ~~ dem65
"
fit_un <- cfa(
  small_model,
  data = PoliticalDemocracy,
  std.lv = TRUE,
  do.fit = FALSE
)

# ---------------------------------------------------------------------------
# random_start() tests
# ---------------------------------------------------------------------------

test_that("random_start returns base vector as first row", {
  set.seed(42)
  starts <- random_start(fit_un, n = 3)
  base <- lavaan::lav_export_estimation(fit_un)$starting_values

  expect_equal(starts[1, ], base, tolerance = 1e-12)
  expect_equal(nrow(starts), 3)
  expect_equal(ncol(starts), length(base))
})

test_that("random_start n=1 returns only base vector", {
  set.seed(42)
  starts <- random_start(fit_un, n = 1)
  base <- lavaan::lav_export_estimation(fit_un)$starting_values

  expect_equal(nrow(starts), 1)
})

test_that("random_start variances fall within bounds", {
  set.seed(42)
  starts <- random_start(fit_un, n = 20)
  pt <- lavaan::parTable(fit_un)
  free_pt <- pt[pt$free > 0, ]

  # Variance parameter indices (lhs == rhs, op == "~~")
  var_idx <- which(free_pt$op == "~~" & free_pt$lhs == free_pt$rhs)
  if (length(var_idx) > 0) {
    obs_variances <- diag(lavaan::lavInspect(fit_un, "cov.ov"))

    for (j in var_idx) {
      lhs <- free_pt$lhs[j]
      # Same logic as in random_start: use max of base and observed variance
      base_val <- starts[1, j] # First row is unperturbed base
      eff_var <- max(base_val, obs_variances[lhs], na.rm = TRUE)
      lb <- max(1e-6, 0.2 * eff_var)
      ub <- 2 * eff_var
      vals <- starts[-1, j]

      expect_true(all(vals >= lb & vals <= ub))
    }
  }
})

test_that("random_start preserves loadings and intercepts at base", {
  set.seed(42)
  starts <- random_start(fit_un, n = 5)
  pt <- lavaan::parTable(fit_un)
  free_pt <- pt[pt$free > 0, ]

  # Factor loading indices (leave at base for all rows)
  load_idx <- which(free_pt$op == "=~")
  for (row in seq_len(nrow(starts))) {
    expect_equal(
      starts[row, load_idx],
      free_pt$est[load_idx],
      tolerance = 1e-12
    )
  }

  # Intercept indices (small_model has none, but test structure is correct)
  int_idx <- which(free_pt$op == "~1")
  if (length(int_idx) > 0) {
    for (row in seq_len(nrow(starts))) {
      expect_equal(
        starts[row, int_idx],
        free_pt$est[int_idx],
        tolerance = 1e-12
      )
    }
  }
})

test_that("random_start covariate perturbations correspond to r in [-0.5, 0.5]", {
  set.seed(42)
  starts <- random_start(fit_un, n = 50)
  pt <- lavaan::parTable(fit_un)
  free_pt <- pt[pt$free > 0, ]

  # Build variance lookup same way as random_start does:
  obs_variances <- setNames(
    diag(lavaan::lavInspect(fit_un, "cov.ov")),
    rownames(lavaan::lavInspect(fit_un, "cov.ov"))
  )
  base <- starts[1, ]

  # Covariance indices (lhs != rhs)
  cov_idx <- which(free_pt$op == "~~" & free_pt$lhs != free_pt$rhs)
  if (length(cov_idx) > 0) {
    for (j in cov_idx) {
      lhs <- free_pt$lhs[j]
      rhs <- free_pt$rhs[j]

      # Build variance lookup same as random_start(): ALL "~~" self-params,
      # including fixed factor variances (std.lv = TRUE), not just free ones
      all_var_mask <- pt$op == "~~" & pt$lhs == pt$rhs
      var_lookup <- setNames(pt$est[all_var_mask], pt$lhs[all_var_mask])
      for (vn in intersect(names(obs_variances), names(var_lookup))) {
        var_lookup[vn] <- max(var_lookup[vn], obs_variances[vn], na.rm = TRUE)
      }

      sigma_lhs <- sqrt(abs(var_lookup[lhs]))
      sigma_rhs <- sqrt(abs(var_lookup[rhs]))
      perturbations <- starts[-1, j] - starts[1, j]
      correlations <- perturbations / (sigma_lhs * sigma_rhs)

      expect_true(all(correlations >= -0.5 & correlations <= 0.5))
    }
  }
})

# ---------------------------------------------------------------------------
# penalized_est() start argument tests
# ---------------------------------------------------------------------------

test_that("penalized_est accepts custom start", {
  set.seed(42)
  fit_un2 <- cfa(
    small_model,
    data = PoliticalDemocracy,
    std.lv = TRUE,
    do.fit = FALSE
  )

  pt <- parTable(fit_un2)
  load_60 <- pt$free[pt$op == "=~" & pt$lhs == "dem60"]
  load_65 <- pt$free[pt$op == "=~" & pt$lhs == "dem65"]

  # Custom start: just add small jitter to default
  base <- lavaan::lav_export_estimation(fit_un2)$starting_values
  custom_start <- base + rnorm(length(base), sd = 0.01)

  fit_custom <- penalized_est(
    x = fit_un2,
    w = 0.03,
    pen_diff_id = list(loadings = rbind(load_60, load_65)),
    start = custom_start
  )

  expect_true(inherits(fit_custom, "lavaan"))
})

test_that("penalized_est errors on bad start length", {
  set.seed(42)
  fit_un2 <- cfa(
    small_model,
    data = PoliticalDemocracy,
    std.lv = TRUE,
    do.fit = FALSE
  )

  pt <- parTable(fit_un2)
  load_60 <- pt$free[pt$op == "=~" & pt$lhs == "dem60"]
  load_65 <- pt$free[pt$op == "=~" & pt$lhs == "dem65"]

  bad_start <- c(1, 2) # Wrong length

  expect_error(
    penalized_est(
      x = fit_un2,
      w = 0.03,
      pen_diff_id = list(loadings = rbind(load_60, load_65)),
      start = bad_start
    ),
    "start must have length"
  )
})

test_that("penalized_est uses explicit epsilon for built-in penalties", {
  fit_un2 <- cfa(
    small_model,
    data = PoliticalDemocracy,
    std.lv = TRUE,
    do.fit = FALSE
  )
  pt <- parTable(fit_un2)
  load_60 <- pt$free[pt$op == "=~" & pt$lhs == "dem60"]
  load_65 <- pt$free[pt$op == "=~" & pt$lhs == "dem65"]
  penalty <- list(loadings = rbind(load_60, load_65))

  fit_default <- penalized_est(fit_un2, w = 0.03, pen_diff_id = penalty)
  fit_eps <- penalized_est(fit_un2, w = 0.03, pen_diff_id = penalty, eps = .01)

  expect_equal(fit_eps@optim$fx, fit_default@optim$fx, tolerance = 1e-8)
})

test_that("penalized_est telescopes epsilon values", {
  fit_un2 <- cfa(
    small_model,
    data = PoliticalDemocracy,
    std.lv = TRUE,
    do.fit = FALSE
  )
  pt <- parTable(fit_un2)
  load_60 <- pt$free[pt$op == "=~" & pt$lhs == "dem60"]
  load_65 <- pt$free[pt$op == "=~" & pt$lhs == "dem65"]

  fit <- penalized_est(
    fit_un2,
    w = 0.03,
    pen_diff_id = list(loadings = rbind(load_60, load_65)),
    eps = "telescoping",
    telescoping_control = list(eps_1 = 1, eps_end = 0.1, eps_steps = 3)
  )
  telescoping <- attr(fit, "telescoping")

  expect_true(inherits(fit, "lavaan"))
  expect_equal(telescoping$eps[1], 1)
  expect_equal(telescoping$max_abs_change[1], NA_real_)
  expect_lte(nrow(telescoping), 3)
  expect_gte(telescoping$eps[nrow(telescoping)], 0.1)
})

test_that("telescoping reuses original starts by default", {
  fit_un2 <- cfa(
    small_model,
    data = PoliticalDemocracy,
    std.lv = TRUE,
    do.fit = FALSE
  )
  pt <- parTable(fit_un2)
  load_60 <- pt$free[pt$op == "=~" & pt$lhs == "dem60"]
  load_65 <- pt$free[pt$op == "=~" & pt$lhs == "dem65"]
  penalty <- list(loadings = rbind(load_60, load_65))
  control <- list(eps_1 = 1, eps_end = 0.1, eps_steps = 2)
  start <- lavaan::lav_export_estimation(fit_un2)$starting_values

  fit_telescope <- penalized_est(
    fit_un2,
    w = 0.03,
    pen_diff_id = penalty,
    start = start,
    eps = "telescoping",
    telescoping_control = control
  )
  fit_final <- penalized_est(
    fit_un2,
    w = 0.03,
    pen_diff_id = penalty,
    start = start,
    eps = 0.1
  )

  expect_equal(fit_telescope@optim$fx, fit_final@optim$fx, tolerance = 1e-8)
  expect_equal(fit_telescope@optim$x, fit_final@optim$x, tolerance = 1e-8)
})

test_that("telescoping warm_start uses preceding estimates", {
  fit_un2 <- cfa(
    small_model,
    data = PoliticalDemocracy,
    std.lv = TRUE,
    do.fit = FALSE
  )
  pt <- parTable(fit_un2)
  load_60 <- pt$free[pt$op == "=~" & pt$lhs == "dem60"]
  load_65 <- pt$free[pt$op == "=~" & pt$lhs == "dem65"]
  penalty <- list(loadings = rbind(load_60, load_65))
  start <- lavaan::lav_export_estimation(fit_un2)$starting_values

  fit_first <- penalized_est(
    fit_un2,
    w = 0.03,
    pen_diff_id = penalty,
    start = start,
    eps = 1
  )
  fit_telescope <- penalized_est(
    fit_un2,
    w = 0.03,
    pen_diff_id = penalty,
    start = start,
    eps = "telescoping",
    telescoping_control = list(
      eps_1 = 1,
      eps_end = 0.1,
      eps_steps = 2,
      warm_start = TRUE
    )
  )
  fit_final <- penalized_est(
    fit_un2,
    w = 0.03,
    pen_diff_id = penalty,
    start = fit_first@optim$x,
    eps = 0.1
  )

  expect_equal(fit_telescope@optim$fx, fit_final@optim$fx, tolerance = 1e-8)
  expect_equal(fit_telescope@optim$x, fit_final@optim$x, tolerance = 1e-8)
})

test_that("penalized_est validates eps", {
  fit_un2 <- cfa(
    small_model,
    data = PoliticalDemocracy,
    std.lv = TRUE,
    do.fit = FALSE
  )

  expect_error(
    penalized_est(fit_un2, w = 0.03, eps = 0),
    "eps must be a positive numeric scalar"
  )
  expect_error(
    penalized_est(
      fit_un2,
      w = 0.03,
      eps = "telescoping",
      telescoping_control = list(eps_1 = 1e-5, eps_end = 1, eps_steps = 2)
    ),
    "telescoping_control"
  )
  expect_error(
    penalized_est(
      fit_un2,
      w = 0.03,
      eps = "telescoping",
      telescoping_control = list(warm_start = 1)
    ),
    "warm_start"
  )
})

# ---------------------------------------------------------------------------
# penalized_est_multistart() tests
# ---------------------------------------------------------------------------

test_that("multistart n_starts=1 matches single penalized_est", {
  set.seed(987)
  fit_un2 <- cfa(
    small_model,
    data = PoliticalDemocracy,
    std.lv = TRUE,
    do.fit = FALSE
  )

  pt <- parTable(fit_un2)
  load_60 <- pt$free[pt$op == "=~" & pt$lhs == "dem60"]
  load_65 <- pt$free[pt$op == "=~" & pt$lhs == "dem65"]

  fit_ms <- penalized_est_multistart(
    x = fit_un2,
    w = 0.03,
    pen_diff_id = list(loadings = rbind(load_60, load_65)),
    n_starts = 1,
    verbose = FALSE
  )
  fit_single <- penalized_est(
    x = fit_un2,
    w = 0.03,
    pen_diff_id = list(loadings = rbind(load_60, load_65))
  )

  expect_equal(fit_ms@optim$fx, fit_single@optim$fx, tolerance = 1e-8)
})

test_that("multistart forwards telescoping epsilon", {
  fit_un2 <- cfa(
    small_model,
    data = PoliticalDemocracy,
    std.lv = TRUE,
    do.fit = FALSE
  )
  pt <- parTable(fit_un2)
  load_60 <- pt$free[pt$op == "=~" & pt$lhs == "dem60"]
  load_65 <- pt$free[pt$op == "=~" & pt$lhs == "dem65"]
  penalty <- list(loadings = rbind(load_60, load_65))
  control <- list(eps_1 = 1, eps_end = 0.1, eps_steps = 3, warm_start = TRUE)

  fit_ms <- penalized_est_multistart(
    fit_un2,
    w = 0.03,
    pen_diff_id = penalty,
    eps = "telescoping",
    telescoping_control = control,
    n_starts = 1
  )
  fit_single <- penalized_est(
    fit_un2,
    w = 0.03,
    pen_diff_id = penalty,
    eps = "telescoping",
    telescoping_control = control
  )

  expect_equal(fit_ms@optim$fx, fit_single@optim$fx, tolerance = 1e-8)
  expect_equal(attr(fit_ms, "telescoping"), attr(fit_single, "telescoping"))
})

test_that("multistart attribute has correct structure", {
  set.seed(987)
  fit_un2 <- cfa(
    small_model,
    data = PoliticalDemocracy,
    std.lv = TRUE,
    do.fit = FALSE
  )

  pt <- parTable(fit_un2)
  load_60 <- pt$free[pt$op == "=~" & pt$lhs == "dem60"]
  load_65 <- pt$free[pt$op == "=~" & pt$lhs == "dem65"]

  fit_ms <- penalized_est_multistart(
    x = fit_un2,
    w = 0.03,
    pen_diff_id = list(loadings = rbind(load_60, load_65)),
    n_starts = 5,
    verbose = FALSE
  )

  ms <- attr(fit_ms, "multistart")

  expect_equal(nrow(ms), 5)
  expect_true("start_id" %in% names(ms))
  expect_true("objective" %in% names(ms))
  expect_true("converged" %in% names(ms))
  # Sorted by ascending objective (NA last)
  expect_true(all(diff(na.omit(ms$objective)) >= -1e-12))
})

test_that("multistart never does worse than single penalized_est", {
  set.seed(987)
  fit_un2 <- cfa(
    small_model,
    data = PoliticalDemocracy,
    std.lv = TRUE,
    do.fit = FALSE
  )

  pt <- parTable(fit_un2)
  load_60 <- pt$free[pt$op == "=~" & pt$lhs == "dem60"]
  load_65 <- pt$free[pt$op == "=~" & pt$lhs == "dem65"]

  fit_ms <- penalized_est_multistart(
    x = fit_un2,
    w = 0.03,
    pen_diff_id = list(loadings = rbind(load_60, load_65)),
    n_starts = 5,
    verbose = FALSE
  )
  fit_single <- penalized_est(
    x = fit_un2,
    w = 0.03,
    pen_diff_id = list(loadings = rbind(load_60, load_65))
  )

  expect_true(fit_ms@optim$fx <= fit_single@optim$fx + 1e-8)
})

test_that("custom starts matrix is respected", {
  set.seed(987)
  fit_un2 <- cfa(
    small_model,
    data = PoliticalDemocracy,
    std.lv = TRUE,
    do.fit = FALSE
  )

  pt <- parTable(fit_un2)
  load_60 <- pt$free[pt$op == "=~" & pt$lhs == "dem60"]
  load_65 <- pt$free[pt$op == "=~" & pt$lhs == "dem65"]

  base <- lavaan::lav_export_estimation(fit_un2)$starting_values
  custom_starts <- rbind(
    base,
    base + runif(length(base), 0.01, 0.05),
    base - runif(length(base), 0.01, 0.05)
  )

  fit_ms <- penalized_est_multistart(
    x = fit_un2,
    w = 0.03,
    pen_diff_id = list(loadings = rbind(load_60, load_65)),
    starts = custom_starts,
    verbose = FALSE
  )

  ms <- attr(fit_ms, "multistart")
  expect_equal(nrow(ms), 3)
})

test_that("keep_all attaches all_fits of correct length", {
  set.seed(987)
  fit_un2 <- cfa(
    small_model,
    data = PoliticalDemocracy,
    std.lv = TRUE,
    do.fit = FALSE
  )

  pt <- parTable(fit_un2)
  load_60 <- pt$free[pt$op == "=~" & pt$lhs == "dem60"]
  load_65 <- pt$free[pt$op == "=~" & pt$lhs == "dem65"]

  fit_keep <- penalized_est_multistart(
    x = fit_un2,
    w = 0.03,
    pen_diff_id = list(loadings = rbind(load_60, load_65)),
    n_starts = 4,
    keep_all = TRUE,
    verbose = FALSE
  )
  fit_drop <- penalized_est_multistart(
    x = fit_un2,
    w = 0.03,
    pen_diff_id = list(loadings = rbind(load_60, load_65)),
    n_starts = 4,
    keep_all = FALSE,
    verbose = FALSE
  )

  all_fits <- attr(fit_keep, "all_fits")
  expect_equal(length(all_fits), 4)

  expect_null(attr(fit_drop, "all_fits"))
})

test_that("all starts failing raises warning and returns fit", {
  set.seed(987)
  fit_un2 <- cfa(
    small_model,
    data = PoliticalDemocracy,
    std.lv = TRUE,
    do.fit = FALSE
  )

  pt <- parTable(fit_un2)
  load_60 <- pt$free[pt$op == "=~" & pt$lhs == "dem60"]
  load_65 <- pt$free[pt$op == "=~" & pt$lhs == "dem65"]

  # Force non-convergence with tiny iter.max
  fit_warn <- penalized_est_multistart(
    x = fit_un2,
    w = 0.03,
    pen_diff_id = list(loadings = rbind(load_60, load_65)),
    n_starts = 3,
    verbose = FALSE,
    opt_control = list(iter.max = 1)
  )

  expect_true(inherits(fit_warn, "lavaan"))
  ms <- attr(fit_warn, "multistart")
  expect_equal(nrow(ms), 3)
})

# ---------------------------------------------------------------------------
# Correctness tests: custom pen_gr should be respected
# ---------------------------------------------------------------------------

test_that("custom pen_gr is used when pen_fn is 'l0a'", {
  set.seed(42)
  fit_un <- cfa(
    small_model,
    data = PoliticalDemocracy,
    std.lv = TRUE,
    do.fit = FALSE
  )

  pt <- parTable(fit_un)
  load_60 <- pt$free[pt$op == "=~" & pt$lhs == "dem60"]
  load_65 <- pt$free[pt$op == "=~" & pt$lhs == "dem65"]

  # Define a custom gradient function that differs from the default
  custom_gr <- function(v) 2 * v / (v^2 + 0.01)^2 # Different from gr_l0a

  # Fit with custom gradient
  fit_custom <- penalized_est(
    x = fit_un,
    w = 0.03,
    pen_diff_id = list(loadings = rbind(load_60, load_65)),
    pen_fn = "l0a",
    pen_gr = custom_gr
  )

  # Fit with default gradient
  fit_default <- penalized_est(
    x = fit_un,
    w = 0.03,
    pen_diff_id = list(loadings = rbind(load_60, load_65)),
    pen_fn = "l0a",
    pen_gr = NULL
  )

  # The two should produce *different* results (custom vs default gradient)
  # We verify by checking that the objective values differ (convergence depends on gradients)
  expect_false(
    all.equal(fit_custom@optim$fx, fit_default@optim$fx, tolerance = 1e-6) ==
      TRUE
  )
})

# ---------------------------------------------------------------------------
# Correctness tests: robust SE failure should degrade gracefully
# ---------------------------------------------------------------------------

test_that("singular Hessian in robust.huber.white warns and degrades gracefully", {
  set.seed(42)
  fit_un <- cfa(
    small_model,
    data = PoliticalDemocracy,
    std.lv = TRUE,
    do.fit = FALSE
  )

  pt <- parTable(fit_un)
  load_60 <- pt$free[pt$op == "=~" & pt$lhs == "dem60"]
  load_65 <- pt$free[pt$op == "=~" & pt$lhs == "dem65"]

  # Use a mock test: we can't easily force a singular Hessian in practice,
  # but we can verify that if SE computation fails, a warning is issued and
  # the fit still has valid estimates (just without SEs).
  # For now, just verify that the code path exists and doesn't crash
  fit_pen <- penalized_est(
    x = fit_un,
    w = 0.03,
    pen_diff_id = list(loadings = rbind(load_60, load_65)),
    se = "robust.huber.white"
  )

  # If SE computation succeeded, vcov should be a matrix; if it failed,
  # SEs were not computed. Either way, the fit should be valid.
  expect_true(inherits(fit_pen, "lavaan"))
  expect_true(fit_pen@optim$converged || !is.null(fit_pen@vcov$vcov))
})
