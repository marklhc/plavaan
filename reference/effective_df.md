# Effective Number of Parameters and Degrees of Freedom for Penalized Fits

Computes the effective number of model parameters, and hence the
effective model degrees of freedom, of a penalized lavaan fit using the
soft-count property of the l0a penalty: a directly penalized parameter
with estimate `z` is counted as `l0a(z, eps) = z^2 / (z^2 + eps)`
effective parameters, which is close to 1 for large `|z|` and close to 0
for small `|z|`. For a difference penalty, each column contributes one
shared value (the invariance baseline) plus the l0a soft count of the
pairwise differences in that column.

## Usage

``` r
effective_df(
  x,
  pen_par_id = NULL,
  pen_diff_id = NULL,
  pen_fn = "l0a",
  eps = 0.01,
  ...
)
```

## Arguments

- x:

  A lavaan model object, typically the return value of
  [`penalized_est()`](https://marklhc.github.io/plavaan/reference/penalized_est.md).
  For such objects, the penalty specification is read from the object's
  ‘penalized’ attribute.

- pen_par_id:

  Integer vector of directly penalized free-parameter IDs, in the same
  order as returned by `lavaan::coef()`. For a fit from
  [`penalized_est()`](https://marklhc.github.io/plavaan/reference/penalized_est.md),
  defaults to the value recorded in the object.

- pen_diff_id:

  Named list of integer matrices of free-parameter IDs for difference
  penalties (same structure as in
  [`penalized_est()`](https://marklhc.github.io/plavaan/reference/penalized_est.md)).
  For a fit from
  [`penalized_est()`](https://marklhc.github.io/plavaan/reference/penalized_est.md),
  defaults to the value recorded in the object.

- pen_fn:

  The penalty function: `"l0a"`, `"alf"`, or a function. For a fit from
  [`penalized_est()`](https://marklhc.github.io/plavaan/reference/penalized_est.md),
  defaults to the value recorded in the object.

- eps:

  The epsilon used by the built-in penalties. For a fit from
  [`penalized_est()`](https://marklhc.github.io/plavaan/reference/penalized_est.md),
  defaults to the value recorded in the object.

- ...:

  Additional arguments passed to a user-supplied `pen_fn`.

## Value

A data frame of class `plavaan_efdf` with one row per penalty component
(`"direct penalty"`, and one row per named block of `pen_diff_id`) and a
final `"TOTAL"` row; the row names are the component labels. Columns:

- `npar`:

  Number of (nominal) free parameters covered by the row. For the
  `TOTAL` row, `length(coef(x))`.

- `npar_effective`:

  Effective number of parameters, with penalized parameters soft-counted
  by their penalty values.

- `df_saved`:

  `npar - npar_effective`, the degrees of freedom saved by the penalty
  relative to the nominal model.

The data frame also carries an `info` attribute with `n_stats` (number
of sample statistics), `df_model` (nominal model df, possibly negative),
`df_model_effective`, `w`, `eps`, and `pen_fn`.

## Details

For a fit returned by
[`penalized_est()`](https://marklhc.github.io/plavaan/reference/penalized_est.md),
any argument left at its default (`NULL` for `pen_par_id`/`pen_diff_id`,
`"l0a"` for `pen_fn`, `.01` for `eps`) is filled in from the object's
‘penalized’ attribute; explicitly supplied values take precedence. Note
that a value equal to the default cannot be distinguished from an
omitted one, so passing, e.g., `eps = .01` explicitly will still be
overridden by the recorded value.

An ID that appears in both `pen_par_id` and a `pen_diff_id` matrix is a
hard error, as the parameter would be double-counted.

The effective-df interpretation is calibrated for the l0a penalty, which
takes the value 0 at the origin. With `pen_fn = "alf"`,
`alf(0) = eps^0.25 != 0`, so the soft counts do not reach 0 and the
effective df should be interpreted with care. A custom `pen_fn` is
evaluated as if it had l0a-like 0/1 behavior; the result is approximate.

The nominal model df, `n_stats - length(coef(x))`, may be negative for
penalized (often under-identified) models; the effective model df,
`n_stats - npar_effective`, is the meaningful quantity. The number of
sample statistics `n_stats` is counted per group from the model's
parameter table.

Fit indices based on the effective df (chi-square, CFI, RMSEA, AIC, BIC)
are obtained with
[`fitmeasures()`](https://rdrr.io/pkg/lavaan/man/fitMeasures.html) on a
fit from
[`penalized_est()`](https://marklhc.github.io/plavaan/reference/penalized_est.md),
which freezes the parameters at the penalized estimates and refits. This
fit evaluation is experimental and requires the fit to have been created
with a non-`"none"` `test` (e.g. `test = "Chisq"`); with the default
`test = "none"`,
[`fitmeasures()`](https://rdrr.io/pkg/lavaan/man/fitMeasures.html) is
unavailable.

## Examples

``` r
library(lavaan)
#> This is lavaan 0.7-2
#> lavaan is FREE software! Please report any bugs.

model <- "
  dem60 =~ y1 + y2 + y3 + y4
  dem65 =~ y5 + y6 + y7 + y8
  dem60 ~~ dem65
"
fit <- cfa(model, data = PoliticalDemocracy, std.lv = TRUE, do.fit = FALSE)
pt <- parTable(fit)
load_60 <- pt$free[pt$op == "=~" & pt$lhs == "dem60"]
load_65 <- pt$free[pt$op == "=~" & pt$lhs == "dem65"]

pen <- penalized_est(
  fit, w = 0.03,
  pen_diff_id = list(loadings = rbind(load_60, load_65)),
  eps = .01,
  test = "Chisq"
)
effective_df(pen)
#>          npar npar_effective df_saved
#> loadings    8       4.010374 3.989626
#> TOTAL      17      13.010374 3.989626
#> 
#> n_stats (sample moments):  36
#> nominal model df:  19
#> effective model df:  22.99
#> penalty:  l0a (w = 0.03, eps = 0.01)
fitmeasures(pen, c("chisq", "df", "cfi"))
#> Fit evaluation for penalized fits is experimental; interpret fit indices with caution.
#>                  npar                  fmin                 chisq 
#>                13.000                 0.316                47.429 
#>                    df                pvalue        baseline.chisq 
#>                22.990                 0.002               461.111 
#>           baseline.df       baseline.pvalue                   cfi 
#>                28.000                 0.000                 0.944 
#>                   tli                  nnfi                   rfi 
#>                 0.931                 0.931                 0.875 
#>                   nfi                  pnfi                   ifi 
#>                 0.897                 0.737                 0.944 
#>                   rni                  logl     unrestricted.logl 
#>                 0.944             -1336.287             -1312.572 
#>                   aic                   bic                ntotal 
#>              2698.594              2728.745                75.000 
#>                  bic2                 rmsea        rmsea.ci.lower 
#>              2719.727                 0.119                 0.070 
#>        rmsea.ci.upper        rmsea.ci.level          rmsea.pvalue 
#>                 0.167                 0.900                 0.014 
#>        rmsea.close.h0 rmsea.notclose.pvalue     rmsea.notclose.h0 
#>                 0.050                 0.912                 0.080 
#>                   rmr            rmr_nomean                  srmr 
#>                 0.773                 0.773                 0.071 
#>          srmr_bentler   srmr_bentler_nomean                  crmr 
#>                 0.071                 0.071                 0.060 
#>           crmr_nomean            srmr_mplus     srmr_mplus_nomean 
#>                 0.060                 0.060                 0.060 
#>                   gfi          gfi.ci.lower          gfi.ci.upper 
#>                 0.923                 0.860                 0.971 
#>          gfi.ci.level                 cn_05                 cn_01 
#>                 0.900                56.598                66.821 
#>            gfi_lisrel           agfi_lisrel                  pgfi 
#>                 0.862                 0.784                 0.551 
#>                   mfi                  ecvi 
#>                 0.850                 0.979 
```
