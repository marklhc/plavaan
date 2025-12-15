# Penalized Parameter Estimation for Longitudinal CFA Models

Performs penalized estimation on a lavaan model object by optimizing a
penalized objective function. The function extracts the objective
function from a lavaan model, applies a penalty function to specified
parameters or pairwise differences of parameters, and returns an updated
model with the optimized parameter estimates.

## Usage

``` r
penalized_est(
  x,
  w,
  pen_par_id = NULL,
  pen_diff_id = NULL,
  pen_fn = "l0a",
  pen_gr = NULL,
  se = "none",
  opt_control = list()
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
  same column indicated by the IDs. For matrices with names starting
  with "loading", the log transformation is applied before computing
  differences.

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

## Value

A lavaan model object updated with the penalized parameter estimates.
The returned object includes an attribute `opt_info` containing the
optimization information returned by
[`nlminb()`](https://rdrr.io/r/stats/nlminb.html).

## Details

The function uses [`nlminb()`](https://rdrr.io/r/stats/nlminb.html) to
minimize a penalized objective function that combines the standard
lavaan objective function with a penalty term. Only the parameter
estimates and the log-likelihood should be interpreted. The returned
object was not "fitted" (`do.fit = FALSE`) to avoid users interpreting
the standard errors, which are generally not valid with penalized
estimation. The degrees of freedom may also be inaccurate. If the
optimization does not converge (convergence code != 0), a warning is
issued.

## Warning

The returned object is not fitted using standard ML. Standard errors
reported by [`summary()`](https://rdrr.io/r/base/summary.html) or
[`parameterEstimates()`](https://rdrr.io/pkg/lavaan/man/parameterEstimates.html)
will be missing unless `se = "robust.huber.white"` was specified. Even
then, they are based on an experimental sandwich approximation and
should be interpreted with caution.

## See also

[`lavaan`](https://rdrr.io/pkg/lavaan/man/lavaan.html),
[`nlminb`](https://rdrr.io/r/stats/nlminb.html)

## Examples

``` r
if (FALSE) { # \dontrun{
library(lavaan)

# Fit a longitudinal factor model using PoliticalDemocracy data
ind_mat <- cbind(c("y1", "y2", "y3", "y4"), c("y5", "y6", "y7", "y8"))
fit <- longcfa(ind_mat, lv_names = c("dem60", "dem65"), data = PoliticalDemocracy,
               long_equal = c("loadings", "intercepts"), lag_cov = TRUE)
# Obtain an unidentified model
mod_un <- longcfa_syntax(
    ind_mat, lv_names = c("dem60", "dem65"),
    lag_cov = TRUE,
    free_latvars = TRUE, free_latmeans = TRUE
)
fit_un <- cfa(mod_un, data = PoliticalDemocracy, do.fit = FALSE, std.lv = TRUE,
              start = fit)

# Get parameter IDs for loadings
load_ids <- get_lav_par_id(fit_un, op = "=~", ind_matrix = ind_mat)
int_ids <- get_lav_par_id(fit_un, op = "~1", ind_matrix = ind_mat)

# Apply penalized estimation with alignment loss
pen_fit <- penalized_est(
    x = fit_un,
    w = 0.1,
    pen_diff_id = list(cbind(t(load_ids), t(int_ids))),
    pen_fn = "alf"
)

# Compare parameter estimates
cbind(coef(fit), coef(pen_fit))

# Compare log-likelihoods
c("scalar invariance" = logLik(fit), "penalized" = logLik(pen_fit))
} # }
```
