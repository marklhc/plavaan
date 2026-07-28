# Multistart Penalized Estimation

Repeatedly calls
[`penalized_est()`](https://marklhc.github.io/plavaan/reference/penalized_est.md)
with different starting vectors and returns the solution with the lowest
penalized objective value. This is recommended when using non-convex
penalties (`l0a`, `alf`), which can lead to local optima.

## Usage

``` r
penalized_est_multistart(
  x,
  w,
  pen_par_id = NULL,
  pen_diff_id = NULL,
  pen_fn = "l0a",
  pen_gr = NULL,
  se = "none",
  opt_control = list(),
  n_starts = 10,
  starts = NULL,
  keep_all = FALSE,
  verbose = FALSE
)
```

## Arguments

- x:

  A fitted lavaan model object from which estimation components will be
  extracted.

- w:

  Numeric scalar. Penalty weight (multiplier) applied to the penalty
  terms.

- pen_par_id:

  Integer vector of parameter IDs to apply the penalty function directly
  to, in the same order as returned by `lavaan::coef()` and by
  [`lavaan::partable()`](https://rdrr.io/pkg/lavaan/man/parTable.html),
  with only the free elements.

- pen_diff_id:

  List of matrices containing parameter IDs. For each matrix, the
  penalty is applied to the pairwise differences of parameters in the
  same column indicated by the IDs.

- pen_fn:

  A character string (`"l0a"` or `"alf"`) or a function that computes
  the penalty. Default is `"l0a"`.

- pen_gr:

  A function that computes the gradient of the penalty function. If
  `pen_fn` is `"l0a"` or `"alf"`, this is automatically set.

- se:

  Character string specifying the type of standard errors to compute.
  Options are `"none"` (default; no standard errors) or
  `"robust.huber.white"` (robust sandwich estimator using numerical
  Hessian and first-order information, which is the same as used in the
  `"mlr"` estimator).

- opt_control:

  A list of control parameters passed to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html). Default
  includes `eval.max = 2e4`, `iter.max = 1e4`, and `abs.tol = 1e-20`.

- n_starts:

  Integer. Number of random starting vectors to try. The first start is
  always lavaan's default (unperturbed), so multistart is never worse
  than a single
  [`penalized_est()`](https://marklhc.github.io/plavaan/reference/penalized_est.md)
  call. Default is 10.

- starts:

  Matrix or list of numeric vectors, each with length equal to the
  number of free parameters. If supplied, random generation is bypassed
  entirely and `n_starts` is ignored. A message is printed noting this.

- keep_all:

  Logical. If `TRUE`, attach all fitted lavaan objects as an attribute
  for inspection of alternative local solutions. Default is `FALSE`.

- verbose:

  Logical. If `TRUE`, print progress messages during optimization.
  Default is `FALSE`.

## Value

A lavaan model object (same shape as \[penalized_est()\]'s return), with
additional attributes:

- \`multistart\`:

  A data frame with one row per start, containing columns \`start_id\`,
  \`objective\` (final penalized objective value), and \`converged\`
  (logical). Rows are sorted by ascending objective.

- \`all_fits\`:

  If \`keep_all = TRUE\`, a named list of all fitted lavaan objects, one
  per starting vector.

## Details

Non-convex penalties like `l0a` and `alf` create rugged optimization
surfaces where the optimizer can settle in different local solutions
depending on starting values. Multistart optimization mitigates this
risk by trying several starts and selecting the best.

Starting-value generation mirrors lavaan's `rstarts` scheme but with one
key deviation: regression coefficients are randomized (not fixed at 0 as
in lavaan), because zero-starts may collide with the penalty function in
ways that hinder convergence toward the global optimum. Perturbation
rules by parameter type follow the same logic as
[lavaan::lavOptions](https://rdrr.io/pkg/lavaan/man/lavOptions.html)`rstarts`:
factor loadings and intercepts stay at base values, variances are drawn
within bounds based on observed variance, and regression/covariance
parameters receive a random correlation perturbation scaled to the
covariance scale.

Execution is sequential. For parallel execution, call
[`penalized_est()`](https://marklhc.github.io/plavaan/reference/penalized_est.md)
with `start = ...` directly and use `future.apply`, `parallel`, or
similar.

Users should set their own seed with
[`set.seed()`](https://rdrr.io/r/base/Random.html) for reproducibility.
This function does not call
[`set.seed()`](https://rdrr.io/r/base/Random.html) internally.

## See also

[penalized_est](https://marklhc.github.io/plavaan/reference/penalized_est.md)
