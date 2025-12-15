
<!-- README.md is generated from README.Rmd. Please edit that file -->

# plavaan

<!-- badges: start -->

<!-- badges: end -->

The goal of plavaan is to perform penalized estimation for latent
variable models using lavaan with a differentiable penalty function, as
described in [Robitzsch (2023)](https://www.mdpi.com/1999-4893/16/9/446)
and similar to the approach described in [Asparouhov & Muthén
(2024)](https://www.statmodel.com/download/PML.pdf).

## Installation

You can install the development version of plavaan like so:

``` r
# install.packages("remotes")
remotes::install_github("marklhc/plavaan")
```

## Example

This is a basic example to obtain penalized cross-loadings in a
confirmatory factor analysis model:

``` r
library(lavaan)
#> This is lavaan 0.6-20
#> lavaan is FREE software! Please report any bugs.
library(plavaan)
data(HolzingerSwineford1939)
model <- '
  visual =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
  textual =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
  speed =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
'

# Dry-run to get starting values for the initial un-identified model
fit_dry <- cfa(model, data = HolzingerSwineford1939, std.lv = TRUE, do.fit = FALSE)

# Locate the cross-loadings to be penalized from the parameter table
pt <- parTable(fit_dry)
subset(pt, op == "=~")
#>    id     lhs op rhs user block group free ustart exo label plabel start   est
#> 1   1  visual =~  x1    1     1     1    1     NA   0         .p1. 0.815 0.815
#> 2   2  visual =~  x2    1     1     1    2     NA   0         .p2. 0.515 0.515
#> 3   3  visual =~  x3    1     1     1    3     NA   0         .p3. 0.592 0.592
#> 4   4  visual =~  x4    1     1     1    4     NA   0         .p4. 0.680 0.680
#> 5   5  visual =~  x5    1     1     1    5     NA   0         .p5. 0.769 0.769
#> 6   6  visual =~  x6    1     1     1    6     NA   0         .p6. 0.687 0.687
#> 7   7  visual =~  x7    1     1     1    7     NA   0         .p7. 0.315 0.315
#> 8   8  visual =~  x8    1     1     1    8     NA   0         .p8. 0.335 0.335
#> 9   9  visual =~  x9    1     1     1    9     NA   0         .p9. 0.521 0.521
#> 10 10 textual =~  x1    1     1     1   10     NA   0        .p10. 0.815 0.815
#> 11 11 textual =~  x2    1     1     1   11     NA   0        .p11. 0.515 0.515
#> 12 12 textual =~  x3    1     1     1   12     NA   0        .p12. 0.592 0.592
#> 13 13 textual =~  x4    1     1     1   13     NA   0        .p13. 0.680 0.680
#> 14 14 textual =~  x5    1     1     1   14     NA   0        .p14. 0.769 0.769
#> 15 15 textual =~  x6    1     1     1   15     NA   0        .p15. 0.687 0.687
#> 16 16 textual =~  x7    1     1     1   16     NA   0        .p16. 0.315 0.315
#> 17 17 textual =~  x8    1     1     1   17     NA   0        .p17. 0.335 0.335
#> 18 18 textual =~  x9    1     1     1   18     NA   0        .p18. 0.521 0.521
#> 19 19   speed =~  x1    1     1     1   19     NA   0        .p19. 0.815 0.815
#> 20 20   speed =~  x2    1     1     1   20     NA   0        .p20. 0.515 0.515
#> 21 21   speed =~  x3    1     1     1   21     NA   0        .p21. 0.592 0.592
#> 22 22   speed =~  x4    1     1     1   22     NA   0        .p22. 0.680 0.680
#> 23 23   speed =~  x5    1     1     1   23     NA   0        .p23. 0.769 0.769
#> 24 24   speed =~  x6    1     1     1   24     NA   0        .p24. 0.687 0.687
#> 25 25   speed =~  x7    1     1     1   25     NA   0        .p25. 0.315 0.315
#> 26 26   speed =~  x8    1     1     1   26     NA   0        .p26. 0.335 0.335
#> 27 27   speed =~  x9    1     1     1   27     NA   0        .p27. 0.521 0.521

# Fit the penalized model using plavaan::penalized_est()
fit_pen <- penalized_est(
  fit_dry,
  w = 0.01,  # close to a BIC-like weight for penalty
  pen_par_id = c(4:12, 16:24),  # cross-loadings id in the free column
  pen_fn = "l0a"  # approximate L0 penalty
)
summary(fit_pen)
#> lavaan 0.6-20 ended normally after 113 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        39
#> 
#>   Number of observations                           301
#> 
#> 
#> Parameter Estimates:
#> 
#> 
#> Latent Variables:
#>                    Estimate
#>   visual =~                
#>     x1                0.813
#>     x2                0.545
#>     x3                0.825
#>     x4                0.015
#>     x5               -0.057
#>     x6                0.036
#>     x7               -0.043
#>     x8                0.034
#>     x9                0.382
#>   textual =~               
#>     x1                0.029
#>     x2               -0.012
#>     x3               -0.205
#>     x4                0.975
#>     x5                1.137
#>     x6                0.889
#>     x7                0.016
#>     x8               -0.009
#>     x9               -0.010
#>   speed =~                 
#>     x1                0.011
#>     x2               -0.028
#>     x3                0.007
#>     x4               -0.004
#>     x5                0.003
#>     x6                0.002
#>     x7                0.698
#>     x8                0.773
#>     x9                0.460
#> 
#> Covariances:
#>                    Estimate
#>   visual ~~                
#>     textual           0.481
#>     speed             0.282
#>   textual ~~               
#>     speed             0.203
#> 
#> Variances:
#>                    Estimate
#>    .x1                0.658
#>    .x2                1.091
#>    .x3                0.701
#>    .x4                0.378
#>    .x5                0.412
#>    .x6                0.364
#>    .x7                0.701
#>    .x8                0.408
#>    .x9                0.558
#>     visual            1.000
#>     textual           1.000
#>     speed             1.000
```

Sandwich estimates of standard errors can be obtained using the
`se = "robust.huber.white"` argument in `penalized_est()`, but it is
unclear how valid these are in the presence of penalization.
