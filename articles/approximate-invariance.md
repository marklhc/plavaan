# Approximate Invariance with Penalized Estimation

``` r
library(plavaan)
library(lavaan)
#> This is lavaan 0.6-20
#> lavaan is FREE software! Please report any bugs.
```

``` r
# Load the PoliticalDemocracy data
data("PoliticalDemocracy", package = "lavaan")
# Configural invariance (pretend y7 is not available in dem65),
# and freely estimated latent means and variances
config_mod <- "
  dem60 =~ y1 + y2 + y3 + y4
  dem65 =~ y5 + y6 +      y8
  dem60 ~~ dem65
  dem60 ~~ 1 * dem60
  dem65 ~~ NA * dem65
  dem60 ~ 0
  dem65 ~ NA * 1
  y1 ~~ y5
  y2 ~~ y6
  y4 ~~ y8
"
fit_dry <- cfa(config_mod, data = PoliticalDemocracy, auto.fix.first = FALSE,
               do.fit = FALSE)
```

Penalized estimation can be used to achieve approximate invariance by
penalizing the differences in intercepts and loadings across groups.

``` r
parTable(fit_dry)
#>    id   lhs op   rhs user block group free ustart exo label plabel start   est
#> 1   1 dem60 =~    y1    1     1     1    1     NA   0         .p1. 2.224 2.224
#> 2   2 dem60 =~    y2    1     1     1    2     NA   0         .p2. 2.882 2.882
#> 3   3 dem60 =~    y3    1     1     1    3     NA   0         .p3. 2.347 2.347
#> 4   4 dem60 =~    y4    1     1     1    4     NA   0         .p4. 2.877 2.877
#> 5   5 dem65 =~    y5    1     1     1    5     NA   0         .p5. 1.000 1.000
#> 6   6 dem65 =~    y6    1     1     1    6     NA   0         .p6. 1.545 1.545
#> 7   7 dem65 =~    y8    1     1     1    7     NA   0         .p7. 1.657 1.657
#> 8   8 dem60 ~~ dem65    1     1     1    8     NA   0         .p8. 0.000 0.000
#> 9   9 dem60 ~~ dem60    1     1     1    0      1   0         .p9. 1.000 1.000
#> 10 10 dem65 ~~ dem65    1     1     1    9     NA   0        .p10. 0.050 0.050
#> 11 11 dem60 ~1          1     1     1    0      0   0        .p11. 0.000 0.000
#> 12 12 dem65 ~1          1     1     1   10     NA   0        .p12. 0.000 0.000
#> 13 13    y1 ~~    y5    1     1     1   11     NA   0        .p13. 0.000 0.000
#> 14 14    y2 ~~    y6    1     1     1   12     NA   0        .p14. 0.000 0.000
#> 15 15    y4 ~~    y8    1     1     1   13     NA   0        .p15. 0.000 0.000
#> 16 16    y1 ~~    y1    0     1     1   14     NA   0        .p16. 3.393 3.393
#> 17 17    y2 ~~    y2    0     1     1   15     NA   0        .p17. 7.686 7.686
#> 18 18    y3 ~~    y3    0     1     1   16     NA   0        .p18. 5.310 5.310
#> 19 19    y4 ~~    y4    0     1     1   17     NA   0        .p19. 5.535 5.535
#> 20 20    y5 ~~    y5    0     1     1   18     NA   0        .p20. 3.367 3.367
#> 21 21    y6 ~~    y6    0     1     1   19     NA   0        .p21. 5.612 5.612
#> 22 22    y8 ~~    y8    0     1     1   20     NA   0        .p22. 5.197 5.197
#> 23 23    y1 ~1          0     1     1   21     NA   0        .p23. 5.465 5.465
#> 24 24    y2 ~1          0     1     1   22     NA   0        .p24. 4.256 4.256
#> 25 25    y3 ~1          0     1     1   23     NA   0        .p25. 6.563 6.563
#> 26 26    y4 ~1          0     1     1   24     NA   0        .p26. 4.453 4.453
#> 27 27    y5 ~1          0     1     1   25     NA   0        .p27. 5.136 5.136
#> 28 28    y6 ~1          0     1     1   26     NA   0        .p28. 2.978 2.978
#> 29 29    y8 ~1          0     1     1   27     NA   0        .p29. 4.043 4.043
# Create matrix to indicate the same item loadings across groups/time in columns for
# penalization on the pairwise differences
ld_mat <- rbind(1:4, c(5:6, NA, 7))
int_mat <- rbind(21:24, c(25:26, NA, 27))
fit_pen <- penalized_est(
  fit_dry, w = .03, pen_diff_id = list(loadings = ld_mat, intercepts = int_mat),
  se = "robust.huber.white"
)
summary(fit_pen)
#> lavaan 0.6-20 ended normally after 90 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        27
#> 
#>   Number of observations                            75
#> 
#> 
#> Parameter Estimates:
#> 
#>   Standard errors                             Sandwich
#>   Information bread                           Observed
#>   Observed information based on                Hessian
#> 
#> Latent Variables:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   dem60 =~                                            
#>     y1                2.099    0.213    9.871    0.000
#>     y2                2.971    0.285   10.416    0.000
#>     y3                2.260    0.323    6.991    0.000
#>     y4                2.970    0.219   13.564    0.000
#>   dem65 =~                                            
#>     y5                2.089    0.217    9.617    0.000
#>     y6                2.962    0.292   10.160    0.000
#>     y8                2.993    0.225   13.288    0.000
#> 
#> Covariances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   dem60 ~~                                            
#>     dem65             0.871    0.056   15.616    0.000
#>  .y1 ~~                                               
#>    .y5                0.940    0.454    2.073    0.038
#>  .y2 ~~                                               
#>    .y6                1.748    0.899    1.944    0.052
#>  .y4 ~~                                               
#>    .y8                0.248    0.504    0.491    0.623
#> 
#> Intercepts:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>     dem60             0.000                           
#>     dem65            -0.146    0.073   -2.006    0.045
#>    .y1                5.455    0.290   18.783    0.000
#>    .y2                4.251    0.455    9.351    0.000
#>    .y3                6.563    0.376   17.449    0.000
#>    .y4                4.465    0.376   11.866    0.000
#>    .y5                5.454    0.290   18.803    0.000
#>    .y6                3.411    0.450    7.580    0.000
#>    .y8                4.466    0.376   11.867    0.000
#> 
#> Variances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>     dem60             1.000                           
#>     dem65             0.867    0.094    9.262    0.000
#>    .y1                2.202    0.520    4.236    0.000
#>    .y2                6.470    1.338    4.834    0.000
#>    .y3                5.512    1.087    5.069    0.000
#>    .y4                2.465    0.627    3.930    0.000
#>    .y5                3.007    0.616    4.879    0.000
#>    .y6                3.745    0.801    4.673    0.000
#>    .y8                2.417    0.711    3.401    0.001
```

The penalized estimation finds a solution where the loadings and
intercepts difference are minimized, depending on the penalty weight
`w`.
