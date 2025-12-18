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
library(lavaan)
#> This is lavaan 0.6-20
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
#> lavaan 0.6-20 ended normally after 103 iterations
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
#>     y1                2.105
#>     y2                2.852
#>     y3                2.531
#>     y4                2.905
#>   dem65 =~                 
#>     y5                2.083
#>     y6                2.819
#>     y7                2.602
#>     y8                2.898
#> 
#> Covariances:
#>                    Estimate
#>   dem60 ~~                 
#>     dem65             0.918
#>  .y1 ~~                    
#>    .y5                0.842
#>  .y2 ~~                    
#>    .y6                1.820
#>  .y3 ~~                    
#>    .y7                1.222
#>  .y4 ~~                    
#>    .y8                0.284
#> 
#> Intercepts:
#>                    Estimate
#>     dem60             0.000
#>     dem65            -0.147
#>    .y1                5.456
#>    .y2                4.251
#>    .y3                6.572
#>    .y4                4.460
#>    .y5                5.454
#>    .y6                3.393
#>    .y7                6.572
#>    .y8                4.460
#> 
#> Variances:
#>                    Estimate
#>     dem60             1.000
#>     dem65             0.951
#>    .y1                2.129
#>    .y2                6.632
#>    .y3                5.388
#>    .y4                2.594
#>    .y5                2.816
#>    .y6                4.003
#>    .y7                3.590
#>    .y8                2.457
#> 
```
