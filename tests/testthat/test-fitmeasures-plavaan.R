library(plavaan)
library(lavaan)

data(PoliticalDemocracy, package = "lavaan")

# ---------------------------------------------------------------------------
# Fixtures: three penalized fits, fitted once at the top of the file
# (testthat evaluates the top level of each file once, in a fresh
# environment per file). All fits use deterministic starting values
# (no random starts), so no seed is needed.
# ---------------------------------------------------------------------------

# T1: under-identified EFA with a direct penalty (pen_par_id)
mod_efa <- "
  ind60 =~ x1 + x2 + x3 + y1 + y2 + y3 + y4
  dem60 =~ x1 + x2 + x3 + y1 + y2 + y3 + y4
  ind60 ~~ ind60
  x1 ~~ x2 + x3 + y1 + y2 + y3 + y4
  x2 ~~ x3 + y1 + y2 + y3 + y4
  x3 ~~ y1 + y2 + y3 + y4
  y1 ~~ y2 + y3 + y4
  y2 ~~ y3 + y4
  y3 ~~ y4
"
fit_efa <- cfa(mod_efa, data = PoliticalDemocracy, std.lv = TRUE, do.fit = FALSE)
pen_efa <- penalized_est(fit_efa, w = .03, pen_par_id = c(4:10, 15:35), eps = .01, test = "Chisq")

# T2: longitudinal meanstructure model, difference penalties on loadings
# and intercepts
mod_long <- "
  dem60 =~ y1 + y2 + y3 + y4
  dem65 =~ y5 + y6 + y7 + y8
  dem60 ~~ dem65
  dem60 ~~ 1 * dem60
  dem65 ~~ NA * dem65
  dem60 ~ 0
  dem65 ~ NA * 1
  y1 ~~ y5
  y2 ~~ y6
  y3 ~~ y7
  y4 ~~ y8
"
fit_long <- cfa(
  mod_long,
  data = PoliticalDemocracy,
  std.lv = TRUE,
  meanstructure = TRUE,
  do.fit = FALSE
)
pt_long <- parTable(fit_long)
l60 <- pt_long$free[pt_long$op == "=~" & pt_long$lhs == "dem60"]
l65 <- pt_long$free[pt_long$op == "=~" & pt_long$lhs == "dem65"]
i60 <- pt_long$free[pt_long$op == "~1" & pt_long$lhs %in% c("y1", "y2", "y3", "y4")]
i65 <- pt_long$free[pt_long$op == "~1" & pt_long$lhs %in% c("y5", "y6", "y7", "y8")]
pen_long <- penalized_est(
  fit_long,
  w = 0.03,
  pen_diff_id = list(loadings = rbind(l60, l65), intercepts = rbind(i60, i65)),
  eps = .01,
  test = "Chisq"
)

# T3: two-group model; the two groups carry identical data, so the
# optimal loading differences are exactly zero
d2 <- rbind(PoliticalDemocracy, PoliticalDemocracy)
d2$grp <- factor(rep(c("g1", "g2"), each = nrow(d2) / 2))
fit_grp <- cfa(
  "dem =~ y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8",
  data = d2,
  group = "grp",
  std.lv = TRUE,
  do.fit = FALSE
)
pt_grp <- parTable(fit_grp)
gl1 <- pt_grp$free[pt_grp$op == "=~" & pt_grp$group == 1]
gl2 <- pt_grp$free[pt_grp$op == "=~" & pt_grp$group == 2]
pen_grp <- penalized_est(
  fit_grp,
  w = 0.05,
  pen_diff_id = list(loadings = rbind(gl1, gl2)),
  eps = .01,
  test = "Chisq"
)

# Small model for the caching test (kept separate so the cache state of
# the T1-T3 fits is not affected by test ordering)
mod_small <- "
  dem60 =~ y1 + y2 + y3 + y4
  dem65 =~ y5 + y6 + y7 + y8
  dem60 ~~ dem65
"
fit_small <- cfa(mod_small, data = PoliticalDemocracy, std.lv = TRUE, do.fit = FALSE)
pt_small <- parTable(fit_small)
l60s <- pt_small$free[pt_small$op == "=~" & pt_small$lhs == "dem60"]
l65s <- pt_small$free[pt_small$op == "=~" & pt_small$lhs == "dem65"]
pen_small <- penalized_est(
  fit_small,
  w = 0.03,
  pen_diff_id = list(loadings = rbind(l60s, l65s)),
  eps = .01,
  test = "Chisq"
)

# Same small model as pen_small, but fitted with the default test = "none",
# i.e. with the experimental fit evaluation disabled
pen_small_none <- penalized_est(
  fit_small,
  w = 0.03,
  pen_diff_id = list(loadings = rbind(l60s, l65s)),
  eps = .01
)

# ---------------------------------------------------------------------------
# Class and recorded penalty specification
# ---------------------------------------------------------------------------

test_that("penalized_est returns a plavaan object with the recorded spec", {
  expect_true(isS4(pen_efa))
  expect_true(inherits(pen_efa, "plavaan"))
  expect_true(inherits(pen_efa, "lavaan"))

  pen <- attr(pen_efa, "penalized")
  expect_equal(
    names(pen),
    c("w", "pen_par_id", "pen_diff_id", "pen_fn", "eps", "test")
  )
  expect_equal(pen$w, .03)
  expect_equal(pen$pen_par_id, c(4:10, 15:35))
  expect_null(pen$pen_diff_id)
  expect_equal(pen$pen_fn, "l0a")
  # eps records the final stage's epsilon
  expect_equal(pen$eps, .01)
  expect_equal(pen$test, "Chisq")

  # the cache is a (shared) environment, initially empty
  expect_true(is.environment(attr(pen_efa, "plavaan.cache")))
})

test_that("the fitted object keeps the plain lavaan interface", {
  expect_true(is.finite(pen_efa@optim$fx))
  expect_equal(lavInspect(pen_efa, "nobs"), 75)
  expect_equal(length(coef(pen_efa)), 43)
  expect_equal(sum(parTable(pen_efa)$free > 0), 43)
  expect_true(pen_efa@optim$converged)
})

# ---------------------------------------------------------------------------
# fitmeasures() on plavaan objects
# ---------------------------------------------------------------------------

test_that("fitmeasures reports the fit indices at the effective df", {
  eff <- effective_df(pen_efa)
  expect_message(
    fm <- fitmeasures(
      pen_efa,
      c("chisq", "df", "cfi", "rmsea", "aic", "bic", "pvalue", "npar")
    ),
    "experimental"
  )
  # df and npar are the EFFECTIVE values. Assert consistency with
  # effective_df() (platform-independent: both are n_stats - npar_eff
  # computed from the same estimates).
  expect_equal(unname(fm["df"]), attr(eff, "info")$df_model_effective, tolerance = 1e-8)
  expect_equal(unname(fm["npar"]), as.numeric(round(eff["TOTAL", "npar_effective"])))
  # The pvalue is the chi-square tail probability at the reported chisq/df.
  expect_equal(
    unname(fm["pvalue"]),
    pchisq(unname(fm["chisq"]), unname(fm["df"]), lower.tail = FALSE),
    tolerance = 1e-8
  )
  # Reference values. The penalized optimization converges to a slightly
  # different point on different platforms (BLAS/LAPACK), so the fit indices
  # drift at the 2nd-3rd decimal; use a generous relative tolerance. (A real
  # regression such as a frozen refit at the wrong parameters would move the
  # chi-square ~19x, far outside this band.)
  expect_equal(unname(fm["chisq"]), 14.65, tolerance = 0.2)
  expect_equal(unname(fm["rmsea"]), 0.0613, tolerance = 0.25)
  expect_true(unname(fm["cfi"]) > 0.9 && unname(fm["cfi"]) <= 1)
  expect_true(is.finite(unname(fm["aic"])) && is.finite(unname(fm["bic"])))
})

test_that("fitmeasures on a meanstructure fit reports the effective df", {
  eff <- effective_df(pen_long)
  expect_message(fm <- fitmeasures(pen_long, c("chisq", "df")), "experimental")
  expect_equal(unname(fm["df"]), attr(eff, "info")$df_model_effective, tolerance = 1e-8)
  # See the EFA test for why the reference value uses a generous tolerance.
  expect_equal(unname(fm["chisq"]), 27.0, tolerance = 0.2)
})

test_that("fitmeasures on a group fit uses the total sample size", {
  eff <- effective_df(pen_grp)
  expect_message(fm <- fitmeasures(pen_grp, c("chisq", "df", "bic", "npar")), "experimental")
  expect_equal(unname(fm["df"]), attr(eff, "info")$df_model_effective, tolerance = 1e-8)
  expect_equal(unname(fm["npar"]), as.numeric(round(eff["TOTAL", "npar_effective"])))
  # The two groups carry identical data, so the effective df is (very close
  # to) the nominal 48 (n_stats 88 - 40 effective parameters).
  expect_equal(unname(fm["df"]), 48, tolerance = 0.01)
  # See the EFA test for why the reference value uses a generous tolerance.
  expect_equal(unname(fm["chisq"]), 92.07, tolerance = 0.2)
  # regression guard: BIC uses @loglik$ntotal (150, not the per-group nobs)
  expect_true(is.finite(unname(fm["bic"])))
})

# ---------------------------------------------------------------------------
# Robustness of the frozen refit
# ---------------------------------------------------------------------------

test_that("fit evaluation guards against non-positive effective df", {
  # An under-identified model that saves no degrees of freedom has
  # effective df < 0, so the chi-square p-value is undefined.
  pen_np <- penalized_est(fit_efa, w = 0.03, pen_par_id = NULL, eps = .01,
                          test = "Chisq")
  expect_lt(attr(effective_df(pen_np), "info")$df_model_effective, 0)

  ws <- character(0)
  fm <- withCallingHandlers(
    fitmeasures(pen_np, c("chisq", "df", "pvalue")),
    message = function(c) invokeRestart("muffleMessage"),
    warning = function(c) {
      ws <<- c(ws, conditionMessage(c))
      invokeRestart("muffleWarning")
    }
  )
  # pvalue is NA (not NaN) and a warning is issued
  expect_true(is.na(unname(fm["pvalue"])))
  expect_false(is.nan(unname(fm["pvalue"])))
  expect_true(any(grepl("non-positive", ws)))
})

test_that("no non-positive-df warning when fit evaluation is disabled", {
  # The same under-identified fit with the default test = "none": fit
  # evaluation is disabled, so no chi-square p-value is computed and the
  # non-positive-df warning must not fire (e.g. when summary() triggers the
  # frozen refit).
  pen_np_none <- penalized_est(fit_efa, w = 0.03, pen_par_id = NULL, eps = .01)
  expect_lt(attr(effective_df(pen_np_none), "info")$df_model_effective, 0)

  ws <- character(0)
  withCallingHandlers(
    summary(pen_np_none),
    message = function(c) invokeRestart("muffleMessage"),
    warning = function(c) {
      ws <<- c(ws, conditionMessage(c))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(any(grepl("non-positive", ws)))
})

test_that("frozen refit works for fits created from sample statistics", {
  d8 <- PoliticalDemocracy[, paste0("y", 1:8)]
  fit_ss <- lavaan::lavaan(
    "dem60 =~ y1 + y2 + y3 + y4
     dem65 =~ y5 + y6 + y7 + y8
     dem60 ~~ dem65",
    sample_cov = cov(d8), sample_nobs = nrow(d8),
    std.lv = TRUE, do.fit = FALSE
  )
  pt <- parTable(fit_ss)
  l60 <- pt$free[pt$op == "=~" & pt$lhs == "dem60"]
  l65 <- pt$free[pt$op == "=~" & pt$lhs == "dem65"]
  # The 'Sigma.hat is not positive definite' warning during optimization is a
  # benign lavaan artifact of fitting from these sample statistics; it is
  # unrelated to the assertion below, so it is suppressed.
  pen_ss <- suppressWarnings(penalized_est(
    fit_ss, w = 0.03, pen_diff_id = list(loadings = rbind(l60, l65)),
    eps = .01, test = "Chisq"
  ))
  # Regression guard: before slotSampleStats was passed to the refit, this
  # failed with 'subscript out of bounds'.
  frozen <- suppressWarnings(expect_no_error(plavaan_frozen(pen_ss)))
  expect_s4_class(frozen, "lavaan")
})

# ---------------------------------------------------------------------------
# Caching of the frozen refit
# ---------------------------------------------------------------------------

test_that("fitmeasures caches the frozen refit in the shared plavaan.cache", {
  cache <- attr(pen_small, "plavaan.cache")
  expect_true(is.environment(cache))
  expect_null(cache$frozen)

  expect_message(fm1 <- fitmeasures(pen_small, c("chisq", "df")), "experimental")
  frozen <- cache$frozen
  expect_false(is.null(frozen))

  # The second call is a cache hit: identical values, and the very same
  # frozen object (no refit)
  expect_message(fm2 <- fitmeasures(pen_small, c("chisq", "df")), "experimental")
  expect_equal(fm1, fm2)
  expect_identical(frozen, cache$frozen)
})

# ---------------------------------------------------------------------------
# summary() on plavaan objects
# ---------------------------------------------------------------------------

test_that("summary reports the effective npar and df", {
  msgs <- capture_messages(s <- summary(pen_efa))
  expect_true(any(grepl("effective npar", msgs)))
  # the Chisq fit also carries the experimental fit-evaluation notice
  expect_true(any(grepl("experimental", msgs)))
  expect_true(inherits(s, "lavaan.summary"))
})

test_that("summary header reports the penalized run's optimizer iterations", {
  msgs <- capture_messages(s <- summary(pen_efa))
  expect_true(any(grepl("effective npar", msgs)))
  expect_true(any(grepl("experimental", msgs)))
  out <- paste(capture.output(print(s)), collapse = "\n")
  n <- pen_efa@optim$iterations
  unit <- if (n == 1) "iteration" else "iterations"
  expect_match(out, sprintf("ended normally after %d %s", n, unit), fixed = TRUE)
})

test_that("summary reuses the cached frozen refit", {
  cache <- attr(pen_small, "plavaan.cache")
  expect_message(fitmeasures(pen_small, c("chisq", "df")), "experimental")
  frozen <- cache$frozen
  expect_false(is.null(frozen))
  msgs <- capture_messages(summary(pen_small))
  expect_true(any(grepl("effective npar", msgs)))
  expect_true(any(grepl("experimental", msgs)))
  # summary() must not refit or replace the cached object
  expect_identical(frozen, cache$frozen)
})

# ---------------------------------------------------------------------------
# Default (test = "none"): fit evaluation disabled
# ---------------------------------------------------------------------------

test_that("with the default test = 'none', fitmeasures is unavailable", {
  expect_equal(attr(pen_small_none, "penalized")$test, "none")
  expect_message(
    fm <- fitmeasures(pen_small_none, c("chisq", "df")),
    "disabled by default"
  )
  expect_null(fm)
})

test_that("with the default test = 'none', summary skips the chi-square test", {
  msgs <- capture_messages(s <- summary(pen_small_none))
  expect_true(any(grepl("effective npar", msgs)))
  # no experimental notice when fit evaluation is disabled
  expect_false(any(grepl("experimental", msgs)))
  expect_true(inherits(s, "lavaan.summary"))
  # and the printed summary carries no chi-square "Model Test" section
  out <- capture.output(print(s))
  expect_false(any(grepl("Test statistic", out)))
})

# ---------------------------------------------------------------------------
# Validation of the test argument
# ---------------------------------------------------------------------------

test_that("penalized_est validates the test argument", {
  # non-strings, wrong lengths, NAs, and values outside the allow-list
  # ('none', 'Chisq', 'SatorraBentler') are all rejected with the intended
  # message (NA_character_ in particular must hit stop(), not error with an
  # invalid `if` condition)
  for (bad in list(5, "", NA, NA_character_, c("a", "b"), "Foo", "chisq")) {
    expect_error(
      penalized_est(
        fit_small,
        w = 0.03,
        pen_diff_id = list(loadings = rbind(l60s, l65s)),
        eps = .01,
        test = bad
      ),
      "test must be one of"
    )
  }
  # a valid non-default value is accepted and recorded on the object
  pen <- penalized_est(
    fit_small,
    w = 0.03,
    pen_diff_id = list(loadings = rbind(l60s, l65s)),
    eps = .01,
    test = "Chisq"
  )
  expect_equal(attr(pen, "penalized")$test, "Chisq")
})
