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

test_that("random_start perturbs only directly penalized parameters", {
  set.seed(219)
  pt <- lavaan::parTable(fit_un)
  load_idx <- pt$free[pt$op == "=~"]
  pen_ids <- load_idx[seq_len(2)]
  starts <- random_start(fit_un, n = 20, pen_par_id = pen_ids)
  base <- starts[1, ]
  unpenalized <- setdiff(seq_along(base), pen_ids)

  expect_equal(starts[-1, unpenalized, drop = FALSE],
    matrix(base[unpenalized], nrow = 19, ncol = length(unpenalized), byrow = TRUE)
  )
  expect_true(any(starts[-1, pen_ids] != base[pen_ids]))
})

test_that("random_start uses non-missing difference-penalty IDs", {
  set.seed(746)
  pt <- lavaan::parTable(fit_un)
  load_60 <- pt$free[pt$op == "=~" & pt$lhs == "dem60"]
  load_65 <- pt$free[pt$op == "=~" & pt$lhs == "dem65"]
  penalty <- list(loadings = rbind(load_60, c(load_65[1:2], NA, load_65[4])))
  pen_ids <- sort(unique(as.numeric(penalty$loadings[!is.na(penalty$loadings)])))
  starts <- random_start(fit_un, n = 20, pen_diff_id = penalty)
  base <- starts[1, ]
  unpenalized <- setdiff(seq_along(base), pen_ids)

  expect_equal(starts[-1, unpenalized, drop = FALSE],
    matrix(base[unpenalized], nrow = 19, ncol = length(unpenalized), byrow = TRUE)
  )
  expect_true(any(starts[-1, pen_ids] != base[pen_ids]))
})

test_that("random_start validates penalty parameter IDs", {
  expect_error(
    random_start(fit_un, n = 2, pen_par_id = 0),
    "positive integer free-parameter indices"
  )
  expect_error(
    random_start(fit_un, n = 2, pen_par_id = 1.5),
    "positive integer free-parameter indices"
  )
  expect_error(
    random_start(fit_un, n = 2, pen_diff_id = list(c(1, 2))),
    "list of matrices"
  )
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
  warns <- capture_warnings(
    fit_warn <- penalized_est_multistart(
      x = fit_un2,
      w = 0.03,
      pen_diff_id = list(loadings = rbind(load_60, load_65)),
      n_starts = 3,
      verbose = FALSE,
      opt_control = list(iter.max = 1)
    )
  )

  expect_true(inherits(fit_warn, "lavaan"))
  ms <- attr(fit_warn, "multistart")
  expect_equal(nrow(ms), 3)
  # Each start warns about non-convergence, and multistart warns that it is
  # returning the best of the non-converged runs.
  expect_length(warns, 4)
  expect_true(all(grepl("did not converge", warns[1:3])))
  expect_true(grepl("None of the 3 optimization runs converged", warns[4]))
})

# ---------------------------------------------------------------------------
# Correctness tests: built-in pen_fn ignores a custom pen_gr
# ---------------------------------------------------------------------------

test_that("built-in pen_fn ignores custom pen_gr with a message", {
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

  # A custom gradient function that differs from the built-in gr_l0a
  custom_gr <- function(v) 2 * v / (v^2 + 0.01)^2

  # Documented behavior: when pen_fn is the built-in "l0a", a user-supplied
  # pen_gr is ignored (with a message) and the built-in gradient is used
  expect_message(
    fit_custom <- penalized_est(
      x = fit_un2,
      w = 0.03,
      pen_diff_id = list(loadings = rbind(load_60, load_65)),
      pen_fn = "l0a",
      pen_gr = custom_gr
    ),
    "pen_gr is ignored"
  )

  fit_default <- penalized_est(
    x = fit_un2,
    w = 0.03,
    pen_diff_id = list(loadings = rbind(load_60, load_65)),
    pen_fn = "l0a"
  )

  # Since the custom pen_gr was ignored, both fits must agree exactly
  expect_equal(fit_custom@optim$fx, fit_default@optim$fx, tolerance = 1e-8)
})

test_that("structural NA in pen_diff_id is excluded from pairwise differences", {
  set.seed(42)
  pt <- parTable(fit_un)
  load_60 <- pt$free[pt$op == "=~" & pt$lhs == "dem60"]
  load_65 <- pt$free[pt$op == "=~" & pt$lhs == "dem65"]
  # Column 2 (y6) is structurally absent in the second row
  pen <- list(loadings = rbind(load_60, c(load_65[1], NA, load_65[3], load_65[4])))

  fit_na <- penalized_est(fit_un, w = 0.03, pen_diff_id = pen)

  expect_true(inherits(fit_na, "lavaan"))
  expect_true(fit_na@optim$converged)

  # The final objective must decompose exactly into the lavaan part plus the
  # penalty over the three *present* pairs, which proves the penalty pairs
  # the right parameters despite the mid-matrix NA.
  x <- fit_na@optim$x
  ff <- lavaan::lav_export_estimation(fit_un)
  lavaan_part <- ff$objective_function(x, lavaan_model = fit_un)
  diffs <- c(
    x[load_60[1]] - x[load_65[1]],
    x[load_60[3]] - x[load_65[3]],
    x[load_60[4]] - x[load_65[4]]
  )
  # rescale = (nrow - 1) / ncombn(nrow, 2) = (2 - 1) / 1 = 1
  intended <- sum(l0a(diffs, eps = 0.01))
  # ignore_attr: lavaan's objective_function() returns fx with metadata attrs
  expect_equal(
    fit_na@optim$fx,
    lavaan_part + 0.03 * intended,
    tolerance = 1e-8,
    ignore_attr = TRUE
  )
})

# ---------------------------------------------------------------------------
# Correctness tests: robust SE failure should degrade gracefully
# ---------------------------------------------------------------------------

test_that("se = 'robust.huber.white' degrades gracefully when sandwich SE fails", {
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
  penalty <- list(loadings = rbind(load_60, load_65))

  # Force the sandwich SE computation to fail (as it would with a singular
  # Hessian) by mocking the internal add_vcov_pen() to error. The fit must
  # then warn and fall back to a fit without standard errors.
  # (capture_warnings() is used instead of expect_warning() because, in
  # testthat edition 3, expect_warning() returns the condition, not the
  # value of the expression.)
  with_mocked_bindings(
    add_vcov_pen = function(fit, hess) stop("mocked singular Hessian"),
    {
      warns <- capture_warnings(
        fit_pen <- penalized_est(
          x = fit_un2,
          w = 0.03,
          pen_diff_id = penalty,
          se = "robust.huber.white"
        )
      )
      expect_true(any(grepl("Computation of robust sandwich", warns)))
    },
    .package = "plavaan"
  )

  expect_true(inherits(fit_pen, "lavaan"))
  expect_true(fit_pen@optim$converged)
  # make_penalized_fit() fallback sets se = "none"
  expect_false(isTRUE(fit_pen@vcov$se == "robust.huber.white"))
})
