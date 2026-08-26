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
pen_efa <- penalized_est(fit_efa, w = .03, pen_par_id = c(4:10, 15:35), eps = .01)

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
  eps = .01
)

# T3: two-group model; the two groups carry identical data, so the
# optimal loading differences are exactly zero
d2 <- rbind(PoliticalDemocracy, PoliticalDemocracy)
d2$grp <- factor(rep(c("g1", "g2"), each = nrow(d2) / 2))
mod_grp <- "dem =~ y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8"
fit_grp <- cfa(mod_grp, data = d2, group = "grp", std.lv = TRUE, do.fit = FALSE)
pt_grp <- parTable(fit_grp)
gl1 <- pt_grp$free[pt_grp$op == "=~" & pt_grp$group == 1]
gl2 <- pt_grp$free[pt_grp$op == "=~" & pt_grp$group == 2]
pen_grp <- penalized_est(
  fit_grp,
  w = 0.05,
  pen_diff_id = list(loadings = rbind(gl1, gl2)),
  eps = .01
)

# Plain (unpenalized) fits for the spec-resolution and validation tests
fit_k4 <- cfa(
  "dem =~ y1 + y2 + y3 + y4",
  data = PoliticalDemocracy,
  std.lv = TRUE,
  do.fit = FALSE
)
fit_k4_ms <- cfa(
  "dem =~ y1 + y2 + y3 + y4",
  data = PoliticalDemocracy,
  std.lv = TRUE,
  meanstructure = TRUE,
  do.fit = FALSE
)

# ---------------------------------------------------------------------------
# effective_df() reference values
# ---------------------------------------------------------------------------

test_that("effective_df reports the reference values for a direct penalty", {
  ed <- effective_df(pen_efa)

  expect_equal(ed["direct penalty", "npar"], 28)
  expect_equal(ed["direct penalty", "npar_effective"], 1.5645, tolerance = 1e-3)

  expect_equal(ed["TOTAL", "npar"], 43)
  expect_equal(ed["TOTAL", "npar_effective"], 16.5645, tolerance = 1e-3)
  expect_equal(ed["TOTAL", "df_saved"], 26.4355, tolerance = 1e-3)

  info <- attr(ed, "info")
  expect_equal(info$n_stats, 28)
  expect_equal(info$df_model, -15) # under-identified nominal model
  expect_equal(info$df_model_effective, 11.4355, tolerance = 1e-3)
  expect_equal(info$w, .03)
  expect_equal(info$eps, .01)
  expect_equal(info$pen_fn, "l0a")
})

test_that("effective_df returns a plavaan_efdf data frame with one row per component", {
  ed <- effective_df(pen_efa)

  expect_true(inherits(ed, "plavaan_efdf"))
  expect_true(inherits(ed, "data.frame"))
  expect_equal(names(ed), c("npar", "npar_effective", "df_saved"))
  expect_equal(rownames(ed), c("direct penalty", "TOTAL"))

  # df_saved is the row-wise complement of npar_effective
  expect_equal(ed$df_saved, ed$npar - ed$npar_effective, tolerance = 1e-12)
  # The TOTAL row counts all free parameters of the fit
  expect_equal(ed["TOTAL", "npar"], sum(parTable(pen_efa)$free > 0))
  expect_equal(ed["TOTAL", "npar"], length(coef(pen_efa)))
})

test_that("effective_df reports one row per difference-penalty block", {
  ed <- effective_df(pen_long)

  expect_equal(rownames(ed), c("loadings", "intercepts", "TOTAL"))
  expect_equal(ed["loadings", "npar"], 8)
  expect_equal(ed["intercepts", "npar"], 8)
  expect_equal(ed["TOTAL", "npar"], 31)
  expect_equal(ed["TOTAL", "npar_effective"], 24.0036, tolerance = 1e-3)

  info <- attr(ed, "info")
  expect_equal(info$n_stats, 44)
  expect_equal(info$df_model, 13)
  expect_equal(info$df_model_effective, 19.9964, tolerance = 1e-3)
})

test_that("block npar_effective is one shared value per column plus the pairwise loss", {
  # Recompute the per-column bookkeeping independently from the estimates
  # (the loadings/intercept IDs are disjoint, so |P| = 16).
  pt <- parTable(pen_long)
  est <- pt$est[pt$free > 0]
  column_loss <- function(ids) {
    sum(vapply(seq_len(ncol(ids)), function(j) {
      composite_pair_loss(est[ids[, j]], fun = l0a, eps = .01)
    }, numeric(1)))
  }
  loadings_loss <- column_loss(rbind(l60, l65))
  intercepts_loss <- column_loss(rbind(i60, i65))

  ed <- effective_df(pen_long)
  # 4 columns -> 4 invariance baseline values, plus the l0a soft count of
  # the pairwise differences in each column.
  expect_equal(ed["loadings", "npar_effective"], 4 + loadings_loss, tolerance = 1e-8)
  expect_equal(ed["intercepts", "npar_effective"], 4 + intercepts_loss, tolerance = 1e-8)
  expect_equal(
    ed["TOTAL", "npar_effective"],
    (length(est) - length(c(l60, l65, i60, i65))) + 4 + loadings_loss + 4 + intercepts_loss,
    tolerance = 1e-8
  )
})

test_that("identical groups save exactly the difference-penalty degrees of freedom", {
  ed <- effective_df(pen_grp)

  expect_equal(ed["loadings", "npar"], 16)
  expect_equal(ed["loadings", "npar_effective"], 8, tolerance = 1e-3)
  expect_equal(ed["loadings", "df_saved"], 8, tolerance = 1e-3)

  expect_equal(ed["TOTAL", "npar"], 48)
  expect_equal(ed["TOTAL", "npar_effective"], 40, tolerance = 1e-3)
  expect_equal(attr(ed, "info")$n_stats, 88)
})

# ---------------------------------------------------------------------------
# plavaan_n_stats()
# ---------------------------------------------------------------------------

test_that("plavaan_n_stats counts the sample moments per group", {
  # k = 4 observed variables: 4 * (4 + 1) / 2 = 10 moments
  expect_equal(plavaan_n_stats(fit_k4), 10)
  # mean structure adds the k means: 10 + 4 = 14
  expect_equal(plavaan_n_stats(fit_k4_ms), 14)
  # two groups with k = 8, mean structure auto-present in group models:
  # 2 * (8 * (8 + 1) / 2 + 8) = 2 * 44 = 88
  expect_equal(plavaan_n_stats(fit_grp), 88)
})

test_that("plavaan_n_stats matches npar + df of a normally fitted group model", {
  fit_grp_fit <- cfa(mod_grp, data = d2, group = "grp", std.lv = TRUE)
  fm <- fitmeasures(fit_grp_fit, c("npar", "df"))
  expect_equal(unname(fm["npar"]), 48)
  expect_equal(unname(fm["df"]), 40)
  expect_equal(unname(fm["npar"]) + unname(fm["df"]), 88)
})

# ---------------------------------------------------------------------------
# Spec resolution
# ---------------------------------------------------------------------------

test_that("args left at their defaults are filled from the recorded spec", {
  ed_def <- effective_df(pen_efa)
  # eps = .01 is indistinguishable from the default, so the recorded value
  # is used and the result is unchanged
  ed_eps <- effective_df(pen_efa, eps = .01)
  expect_equal(ed_eps$npar_effective, ed_def$npar_effective, tolerance = 1e-12)

  # An explicit value at a non-default takes precedence over the recorded
  # one: a larger eps soft-counts fewer parameters
  ed_big <- effective_df(pen_efa, eps = .1)
  expect_lt(
    ed_big["direct penalty", "npar_effective"],
    ed_def["direct penalty", "npar_effective"]
  )
})

test_that("an explicit pen_par_id subset overrides the recorded spec", {
  ed_def <- effective_df(pen_efa)
  ed_sub <- effective_df(pen_efa, pen_par_id = 4:10)

  expect_equal(ed_sub["direct penalty", "npar"], 7)
  expect_lt(
    ed_sub["direct penalty", "npar_effective"],
    ed_def["direct penalty", "npar_effective"]
  )
  expect_equal(ed_sub["TOTAL", "npar"], 43)
})

test_that("a plain lavaan fit without the attribute needs an explicit spec", {
  expect_error(
    effective_df(fit_k4),
    "pen_par_id/pen_diff_id"
  )
  ed <- effective_df(fit_k4, pen_par_id = 1:4)
  expect_equal(ed["direct penalty", "npar"], 4)
  expect_equal(ed["TOTAL", "npar"], sum(parTable(fit_k4)$free > 0))
  # w is not recorded for a non-penalized fit
  expect_true(is.na(attr(ed, "info")$w))
})

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

test_that("an ID in both pen_par_id and pen_diff_id is a hard error", {
  expect_error(
    effective_df(
      fit_k4,
      pen_par_id = 1,
      pen_diff_id = list(b = matrix(c(1, 2), nrow = 2, byrow = TRUE))
    ),
    "appear in both pen_par_id and pen_diff_id"
  )
})

test_that("penalty parameter IDs are validated", {
  # out of range, non-integer, and non-positive
  expect_error(
    effective_df(fit_k4, pen_par_id = 99),
    "positive integer free-parameter indices"
  )
  expect_error(
    effective_df(fit_k4, pen_par_id = 1.5),
    "positive integer free-parameter indices"
  )
  expect_error(
    effective_df(fit_k4, pen_par_id = 0),
    "positive integer free-parameter indices"
  )
  # also through the pen_diff_id matrices
  expect_error(
    effective_df(
      fit_k4,
      pen_diff_id = list(b = matrix(c(1, 99), nrow = 2, byrow = TRUE))
    ),
    "positive integer free-parameter indices"
  )
})

# ---------------------------------------------------------------------------
# Penalty function resolution
# ---------------------------------------------------------------------------

test_that("the alf penalty runs and soft-counts strictly above l0a", {
  ed_l0a <- effective_df(pen_efa)
  ed_alf <- effective_df(pen_efa, pen_fn = "alf")

  expect_true(all(is.finite(ed_alf$npar_effective)))
  expect_true(all(ed_alf$npar_effective > 0))
  # alf(z) = (z^2 + eps)^.25 > 0 at z = 0 while l0a(0) = 0, so every row
  # counts at least as many effective parameters as under l0a
  expect_true(ed_alf["TOTAL", "npar_effective"] > ed_l0a["TOTAL", "npar_effective"])
  expect_equal(attr(ed_alf, "info")$pen_fn, "alf")
})

test_that("a custom penalty function is evaluated on the estimates", {
  # This closure is algebraically identical to l0a(z, eps = .01), so the
  # result must match the default l0a run exactly
  ed_l0a <- effective_df(pen_efa)
  ed_cust <- effective_df(pen_efa, pen_fn = function(z) z^2 / (z^2 + .01))

  expect_equal(ed_cust$npar, ed_l0a$npar)
  expect_equal(ed_cust$npar_effective, ed_l0a$npar_effective, tolerance = 1e-12)
  expect_equal(ed_cust$df_saved, ed_l0a$df_saved, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# print method
# ---------------------------------------------------------------------------

test_that("print.plavaan_efdf prints the nominal and effective model df", {
  expect_output(print(effective_df(pen_efa)), "effective model df")
  expect_output(print(effective_df(pen_efa)), "nominal model df")
  expect_output(
    print(effective_df(pen_long)),
    "n_stats (sample moments)",
    fixed = TRUE
  )
})
