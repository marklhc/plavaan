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
  opt_control = list(),
  start = NULL,
  eps = 0.01,
  telescoping_control = list(eps_1 = 1, eps_end = 1e-05, eps_steps = 20, warm_start =
    FALSE),
  ...
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
  [`lavaan::parTable()`](https://rdrr.io/pkg/lavaan/man/parTable.html),
  with only the free elements.

- pen_diff_id:

  A named list of integer matrices of free-parameter IDs (same order as
  `lavaan::coef()` / the `free` column of
  [`lavaan::parTable()`](https://rdrr.io/pkg/lavaan/man/parTable.html)).
  Each matrix has one row per group or time point and one column per
  matched parameter; the penalty is the sum of pairwise row differences
  within each column, rescaled by `(nrow - 1) / ncombn(nrow, 2)`.
  Structural `NA` entries mark a parameter absent in that row and are
  excluded from the differences.

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

- start:

  Numeric vector of starting values for the optimizer, or `NULL`
  (default) to use lavaan's default starting values. If supplied, its
  length must match the number of free parameters in the model.

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

- ...:

  Additional arguments passed to a user-supplied `pen_fn` / `pen_gr`.
  Custom penalty functions must accept `...`. Built-in penalties
  (`"l0a"`, `"alf"`) ignore it.

## Value

A lavaan model object updated with the penalized parameter estimates.
With `eps = "telescoping"`, it includes a `"telescoping"` data frame
with per-stage epsilon values, parameter changes, objective values, and
convergence indicators.

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

With `eps = "telescoping"`, the model is fit along a log-spaced sequence
from `telescoping_control$eps_1` to `telescoping_control$eps_end`. By
default, each stage uses the original starting values; set
`telescoping_control$warm_start = TRUE` to initialize later stages from
the preceding solution. The sequence stops when the largest absolute
change between consecutive parameter vectors is at most `5e-4`. The
returned object has a `"telescoping"` attribute with stage diagnostics.

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
library(lavaan)
#> This is lavaan 0.7-2
#> lavaan is FREE software! Please report any bugs.

# Define a longitudinal factor model with PoliticalDemocracy data
model <- "
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

# Fit the model without constraints first to get parameter table
fit_un <- cfa(model, data = PoliticalDemocracy, std.lv = TRUE,
              meanstructure = TRUE, do.fit = FALSE)

# Get parameter IDs
pt <- parTable(fit_un)
# Loadings
load_60 <- pt$free[pt$op == "=~" & pt$lhs == "dem60"]
load_65 <- pt$free[pt$op == "=~" & pt$lhs == "dem65"]
# Intercepts
int_60 <- pt$free[pt$op == "~1" & pt$lhs %in% c("y1", "y2", "y3", "y4")]
int_65 <- pt$free[pt$op == "~1" & pt$lhs %in% c("y5", "y6", "y7", "y8")]

# Apply penalized estimation to penalize differences in loadings and intercepts
pen_fit <- penalized_est(
    x = fit_un,
    w = 0.03,
    pen_diff_id = list(
        loadings = rbind(load_60, load_65),
        intercepts = rbind(int_60, int_65)
    ),
    pen_fn = "l0a"
)

# Compare parameter estimates
summary(pen_fit)
#> lavaan 0.7-2 ended normally after 129 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        31
#> 
#>   Number of observations                            75
#> 
#> 
#> Parameter Estimates:
#> 
#> 
#> Latent Variables:
#>                    Estimate
#>   dem60 =~                 
#>     y1                2.098
#>     y2                2.830
#>     y3                2.561
#>     y4                2.901
#>   dem65 =~                 
#>     y5                2.092
#>     y6                2.825
#>     y7                2.572
#>     y8                2.900
#> 
#> Covariances:
#>                    Estimate
#>   dem60 ~~                 
#>     dem65             0.917
#>  .y1 ~~                    
#>    .y5                0.838
#>  .y2 ~~                    
#>    .y6                1.825
#>  .y3 ~~                    
#>    .y7                1.205
#>  .y4 ~~                    
#>    .y8                0.286
#> 
#> Intercepts:
#>                    Estimate
#>     dem60             0.000
#>     dem65            -0.147
#>    .y1                5.457
#>    .y2                4.252
#>    .y3                6.570
#>    .y4                4.461
#>    .y5                5.456
#>    .y6                3.396
#>    .y7                6.570
#>    .y8                4.462
#> 
#> Variances:
#>                    Estimate
#>     dem60             1.000
#>     dem65             0.950
#>    .y1                2.128
#>    .y2                6.658
#>    .y3                5.385
#>    .y4                2.600
#>    .y5                2.813
#>    .y6                4.000
#>    .y7                3.615
#>    .y8                2.456
#> 
```
