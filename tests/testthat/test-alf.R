test_that("alf is computed correctly", {
    x <- matrix(c(-2, -1, 0.5, 1, 2), ncol = 1)
    result <- alf(x, eps = 1e-6)
    expected <- sqrt(abs(x))
    expect_equal(result, expected, tolerance = 1e-6)
})

test_that("composite_pair_loss computes correct sum", {
    x1 <- rbind(
        c(1, 1.2, 1.2),
        c(1, 0.6, 0.6),
        c(1, 1.2, 0.9)
    )
    result <- composite_pair_loss(x1, fun = alf, eps = 1e-16)
    expected <- sum(
        sqrt(abs(x1[1, ] - x1[2, ])),
        sqrt(abs(x1[1, ] - x1[3, ])),
        sqrt(abs(x1[2, ] - x1[3, ]))
    )
    expect_equal(result, expected * 2 / 3, tolerance = 1e-3)
    res1 <- composite_pair_loss(x1[, 1], fun = alf)
    res2 <- composite_pair_loss(x1[, 2], fun = alf)
    res3 <- composite_pair_loss(x1[, 3], fun = alf)
    expect_true(res3 > res2 & res2 > res1)
})

test_that("composite_pair_loss handles missing data properly", {
    x2 <- rbind(
        c(1, NA, 1.2),
        c(1, 0.6, NA),
        c(1, 1.2, 0.9)
    )
    result <- composite_pair_loss(x2, fun = alf, eps = 1e-16)
    expected <- sum(
        sqrt(abs(x2[1, 1] - x2[2, 1])),
        sqrt(abs(x2[1, 3] - x2[3, 3])),
        sqrt(abs(x2[2, 2] - x2[3, 2]))
    )
    expect_equal(result, expected * 2 / 3, tolerance = 1e-3)
})

# ---------------------------------------------------------------------------
# Production-path gradient tests
#
# The analytic penalized gradient used by penalized_est() (penalized_gr(),
# assembled against penalized_obj() exactly as in penalized_est_stage()) is
# validated against numerical differentiation (numDeriv::grad()).
# ---------------------------------------------------------------------------

small_model <- "
  dem60 =~ y1 + y2 + y3 + y4
  dem65 =~ y5 + y6 + y7 + y8
  dem60 ~~ dem65
"
fit_un <- lavaan::cfa(
  small_model,
  data = lavaan::PoliticalDemocracy,
  std.lv = TRUE,
  do.fit = FALSE
)
pt <- lavaan::parTable(fit_un)
load_60 <- pt$free[pt$op == "=~" & pt$lhs == "dem60"]
load_65 <- pt$free[pt$op == "=~" & pt$lhs == "dem65"]

# Build the penalized objective and its analytic gradient the way
# penalized_est_stage() does, for a given (pen_par_id, pen_diff_id) combo.
make_penalized_pair <- function(pen_par_id, pen_diff_id) {
  ff <- lavaan::lav_export_estimation(fit_un)
  p <- ff$starting_values
  eps <- 0.01
  obj_fn <- function(z) ff$objective_function(z, lavaan_model = fit_un)
  gr_fn <- function(z) ff$gradient_function(z, lavaan_model = fit_un)
  pen_fn <- function(z) l0a(z, eps = eps)
  pen_gr <- function(z) gr_l0a(z, eps = eps)
  diff_configs <- make_diff_configs(pen_diff_id)
  list(
    p = p,
    f1 = function(z) penalized_obj(
      z,
      obj_fn = obj_fn,
      w = 0.03,
      pen_fn = pen_fn,
      pen_par_id = pen_par_id,
      diff_configs = diff_configs
    ),
    g = function(z) penalized_gr(
      z,
      gr_fn = gr_fn,
      w = 0.03,
      pen_gr = pen_gr,
      pen_par_id = pen_par_id,
      diff_configs = diff_configs
    )
  )
}

test_that("penalized gradient matches numDeriv for pen_par_id only", {
  s <- make_penalized_pair(pen_par_id = load_60, pen_diff_id = NULL)
  expect_equal(s$g(s$p), numDeriv::grad(s$f1, s$p), tolerance = 1e-6)
})

test_that("penalized gradient matches numDeriv for pen_diff_id only", {
  s <- make_penalized_pair(
    pen_par_id = NULL,
    pen_diff_id = list(loadings = rbind(load_60, load_65))
  )
  expect_equal(s$g(s$p), numDeriv::grad(s$f1, s$p), tolerance = 1e-6)
})

test_that("penalized gradient matches numDeriv with a structural NA in pen_diff_id", {
  s <- make_penalized_pair(
    pen_par_id = NULL,
    pen_diff_id = list(loadings = rbind(
      load_60,
      c(load_65[1], NA, load_65[3], load_65[4])
    ))
  )
  expect_equal(s$g(s$p), numDeriv::grad(s$f1, s$p), tolerance = 1e-6)
})

test_that("penalized gradient matches numDeriv for combined penalties", {
  s <- make_penalized_pair(
    pen_par_id = load_60,
    pen_diff_id = list(loadings = rbind(load_60, load_65))
  )
  expect_equal(s$g(s$p), numDeriv::grad(s$f1, s$p), tolerance = 1e-6)
})
