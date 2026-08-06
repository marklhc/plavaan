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
  eps = 0.01,
  telescoping_control = list(eps_1 = 1, eps_end = 1e-05, eps_steps = 20, warm_start =
    FALSE),
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

- eps:

  A positive numeric scalar used by the built-in penalties, or
  `"telescoping"` to fit a sequence of decreasing epsilon values.
  Default is `.01`. This argument does not alter custom `pen_fn` or
  `pen_gr` functions.

- telescoping_control:

  A named list controlling telescoping, with `eps_1` (default `1`),
  `eps_end` (default `1e-5`), `eps_steps` (default `20`), and
  `warm_start` (default `FALSE`). When `warm_start` is `FALSE`, every
  epsilon stage uses the original starting values; when `TRUE`, each
  stage after the first uses the preceding stage's estimates.

- n_starts:

  Integer. Number of random starting vectors to try. The first start is
  always lavaan's default (unperturbed). Later random starts perturb
  only parameters participating in `pen_par_id` or `pen_diff_id`.
  Default is 10.

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

A lavaan model object (same shape as
[`penalized_est()`](https://marklhc.github.io/plavaan/reference/penalized_est.md)'s
return), with additional attributes:

- `multistart`:

  A data frame with one row per start, containing columns `start_id`,
  `objective` (final penalized objective value), and `converged`
  (logical). Rows are sorted by ascending objective.

- `all_fits`:

  If `keep_all = TRUE`, a named list of all fitted lavaan objects, one
  per starting vector.

## Details

Non-convex penalties like `l0a` and `alf` create rugged optimization
surfaces where the optimizer can settle in different local solutions
depending on starting values. Multistart optimization mitigates this
risk by trying several starts and selecting the best.

Random starts perturb only the free parameters participating in
`pen_par_id` or `pen_diff_id`; all nuisance parameters retain lavaan's
base values. This targets the non-convex penalized surface while
preserving a valid covariance structure from the supplied model start.
Eligible loadings and intercepts are jittered around their base values,
eligible variances are drawn within bounds based on observed variance,
and eligible regression/covariance parameters receive a random
correlation-scale perturbation.

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
