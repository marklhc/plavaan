# Penalized Estimation of Cross-Loadings and Unique Covariances

``` r
library(plavaan)
library(lavaan)
#> This is lavaan 0.6-20
#> lavaan is FREE software! Please report any bugs.
data(PoliticalDemocracy)
```

## Penalize cross-loadings

### Two-factor CFA model

``` r
mod0 <- "
  ind60 =~ x1 + x2 + x3
  dem60 =~ y1 + y2 + y3 + y4
  ind60 ~~ dem60
"
fit0 <- cfa(mod0, data = PoliticalDemocracy, std.lv = TRUE)
```

### Two-factor EFA model (unidentified)

``` r
mod <- "
  ind60 =~ x1 + x2 + x3 + y1 + y2 + y3 + y4
  dem60 =~ x1 + x2 + x3 + y1 + y2 + y3 + y4
  ind60 ~~ ind60
"
fit <- cfa(mod, data = PoliticalDemocracy, std.lv = TRUE, do.fit = FALSE)
```

### Two-factor EFA model with penalized cross-loadings

The cross-loadings are the parameters 4 to 10 in the parameter table
(see the `free` column).

``` r
parTable(fit)
#>    id   lhs op   rhs user block group free ustart exo label plabel start   est
#> 1   1 ind60 =~    x1    1     1     1    1     NA   0         .p1. 0.951 0.951
#> 2   2 ind60 =~    x2    1     1     1    2     NA   0         .p2. 2.001 2.001
#> 3   3 ind60 =~    x3    1     1     1    3     NA   0         .p3. 1.687 1.687
#> 4   4 ind60 =~    y1    1     1     1    4     NA   0         .p4. 1.344 1.344
#> 5   5 ind60 =~    y2    1     1     1    5     NA   0         .p5. 1.807 1.807
#> 6   6 ind60 =~    y3    1     1     1    6     NA   0         .p6. 1.778 1.778
#> 7   7 ind60 =~    y4    1     1     1    7     NA   0         .p7. 2.256 2.256
#> 8   8 dem60 =~    x1    1     1     1    8     NA   0         .p8. 0.951 0.951
#> 9   9 dem60 =~    x2    1     1     1    9     NA   0         .p9. 2.001 2.001
#> 10 10 dem60 =~    x3    1     1     1   10     NA   0        .p10. 1.687 1.687
#> 11 11 dem60 =~    y1    1     1     1   11     NA   0        .p11. 1.344 1.344
#> 12 12 dem60 =~    y2    1     1     1   12     NA   0        .p12. 1.807 1.807
#> 13 13 dem60 =~    y3    1     1     1   13     NA   0        .p13. 1.778 1.778
#> 14 14 dem60 =~    y4    1     1     1   14     NA   0        .p14. 2.256 2.256
#> 15 15 ind60 ~~ ind60    1     1     1    0      1   0        .p15. 1.000 1.000
#> 16 16    x1 ~~    x1    0     1     1   15     NA   0        .p16. 0.265 0.265
#> 17 17    x2 ~~    x2    0     1     1   16     NA   0        .p17. 1.126 1.126
#> 18 18    x3 ~~    x3    0     1     1   17     NA   0        .p18. 0.975 0.975
#> 19 19    y1 ~~    y1    0     1     1   18     NA   0        .p19. 3.393 3.393
#> 20 20    y2 ~~    y2    0     1     1   19     NA   0        .p20. 7.686 7.686
#> 21 21    y3 ~~    y3    0     1     1   20     NA   0        .p21. 5.310 5.310
#> 22 22    y4 ~~    y4    0     1     1   21     NA   0        .p22. 5.535 5.535
#> 23 23 dem60 ~~ dem60    0     1     1    0      1   0        .p23. 1.000 1.000
#> 24 24 ind60 ~~ dem60    0     1     1   22     NA   0        .p24. 0.000 0.000
```

``` r
pefa_fit <- penalized_est(
    fit,
    w = .03,
    pen_par_id = 4:10
)
summary(pefa_fit)
#> lavaan 0.6-20 ended normally after 126 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        22
#> 
#>   Number of observations                            75
#> 
#> 
#> Parameter Estimates:
#> 
#> 
#> Latent Variables:
#>                    Estimate
#>   ind60 =~                 
#>     x1                0.658
#>     x2                1.456
#>     x3                1.222
#>     y1               -0.007
#>     y2               -0.608
#>     y3               -0.001
#>     y4                0.006
#>   dem60 =~                 
#>     x1                0.025
#>     x2               -0.002
#>     x3               -0.010
#>     y1                2.071
#>     y2                3.290
#>     y3                2.256
#>     y4                2.999
#> 
#> Covariances:
#>                    Estimate
#>   ind60 ~~                 
#>     dem60             0.481
#> 
#> Variances:
#>                    Estimate
#>     ind60             1.000
#>    .x1                0.079
#>    .x2                0.127
#>    .x3                0.464
#>    .y1                2.493
#>    .y2                6.048
#>    .y3                5.512
#>    .y4                2.017
#>     dem60             1.000
```

## Penalize Cross-loadings and Unique Covariances

### Two-factor EFA model with unique covariances

``` r
mod2 <- "
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
fit2 <- cfa(mod2, data = PoliticalDemocracy, std.lv = TRUE, do.fit = FALSE)
```

### Two-factor EFA model with penalized cross-loadings and unique covariances

The unique covariances are the parameters 15 to 35 in the parameter
table (see the `free` column).

``` r
parTable(fit2)
#>    id   lhs op   rhs user block group free ustart exo label plabel start   est
#> 1   1 ind60 =~    x1    1     1     1    1     NA   0         .p1. 0.951 0.951
#> 2   2 ind60 =~    x2    1     1     1    2     NA   0         .p2. 2.001 2.001
#> 3   3 ind60 =~    x3    1     1     1    3     NA   0         .p3. 1.687 1.687
#> 4   4 ind60 =~    y1    1     1     1    4     NA   0         .p4. 1.344 1.344
#> 5   5 ind60 =~    y2    1     1     1    5     NA   0         .p5. 1.807 1.807
#> 6   6 ind60 =~    y3    1     1     1    6     NA   0         .p6. 1.778 1.778
#> 7   7 ind60 =~    y4    1     1     1    7     NA   0         .p7. 2.256 2.256
#> 8   8 dem60 =~    x1    1     1     1    8     NA   0         .p8. 0.951 0.951
#> 9   9 dem60 =~    x2    1     1     1    9     NA   0         .p9. 2.001 2.001
#> 10 10 dem60 =~    x3    1     1     1   10     NA   0        .p10. 1.687 1.687
#> 11 11 dem60 =~    y1    1     1     1   11     NA   0        .p11. 1.344 1.344
#> 12 12 dem60 =~    y2    1     1     1   12     NA   0        .p12. 1.807 1.807
#> 13 13 dem60 =~    y3    1     1     1   13     NA   0        .p13. 1.778 1.778
#> 14 14 dem60 =~    y4    1     1     1   14     NA   0        .p14. 2.256 2.256
#> 15 15 ind60 ~~ ind60    1     1     1    0      1   0        .p15. 1.000 1.000
#> 16 16    x1 ~~    x2    1     1     1   15     NA   0        .p16. 0.000 0.000
#> 17 17    x1 ~~    x3    1     1     1   16     NA   0        .p17. 0.000 0.000
#> 18 18    x1 ~~    y1    1     1     1   17     NA   0        .p18. 0.000 0.000
#> 19 19    x1 ~~    y2    1     1     1   18     NA   0        .p19. 0.000 0.000
#> 20 20    x1 ~~    y3    1     1     1   19     NA   0        .p20. 0.000 0.000
#> 21 21    x1 ~~    y4    1     1     1   20     NA   0        .p21. 0.000 0.000
#> 22 22    x2 ~~    x3    1     1     1   21     NA   0        .p22. 0.000 0.000
#> 23 23    x2 ~~    y1    1     1     1   22     NA   0        .p23. 0.000 0.000
#> 24 24    x2 ~~    y2    1     1     1   23     NA   0        .p24. 0.000 0.000
#> 25 25    x2 ~~    y3    1     1     1   24     NA   0        .p25. 0.000 0.000
#> 26 26    x2 ~~    y4    1     1     1   25     NA   0        .p26. 0.000 0.000
#> 27 27    x3 ~~    y1    1     1     1   26     NA   0        .p27. 0.000 0.000
#> 28 28    x3 ~~    y2    1     1     1   27     NA   0        .p28. 0.000 0.000
#> 29 29    x3 ~~    y3    1     1     1   28     NA   0        .p29. 0.000 0.000
#> 30 30    x3 ~~    y4    1     1     1   29     NA   0        .p30. 0.000 0.000
#> 31 31    y1 ~~    y2    1     1     1   30     NA   0        .p31. 0.000 0.000
#> 32 32    y1 ~~    y3    1     1     1   31     NA   0        .p32. 0.000 0.000
#> 33 33    y1 ~~    y4    1     1     1   32     NA   0        .p33. 0.000 0.000
#> 34 34    y2 ~~    y3    1     1     1   33     NA   0        .p34. 0.000 0.000
#> 35 35    y2 ~~    y4    1     1     1   34     NA   0        .p35. 0.000 0.000
#> 36 36    y3 ~~    y4    1     1     1   35     NA   0        .p36. 0.000 0.000
#> 37 37    x1 ~~    x1    0     1     1   36     NA   0        .p37. 0.265 0.265
#> 38 38    x2 ~~    x2    0     1     1   37     NA   0        .p38. 1.126 1.126
#> 39 39    x3 ~~    x3    0     1     1   38     NA   0        .p39. 0.975 0.975
#> 40 40    y1 ~~    y1    0     1     1   39     NA   0        .p40. 3.393 3.393
#> 41 41    y2 ~~    y2    0     1     1   40     NA   0        .p41. 7.686 7.686
#> 42 42    y3 ~~    y3    0     1     1   41     NA   0        .p42. 5.310 5.310
#> 43 43    y4 ~~    y4    0     1     1   42     NA   0        .p43. 5.535 5.535
#> 44 44 dem60 ~~ dem60    0     1     1    0      1   0        .p44. 1.000 1.000
#> 45 45 ind60 ~~ dem60    0     1     1   43     NA   0        .p45. 0.000 0.000
```

``` r
pefa_fit2 <- penalized_est(
    fit2,
    w = .03,
    pen_par_id = c(4:10, 15:35)
)
summary(pefa_fit2)
#> lavaan 0.6-20 ended normally after 195 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        43
#> 
#>   Number of observations                            75
#> 
#> 
#> Parameter Estimates:
#> 
#> 
#> Latent Variables:
#>                    Estimate
#>   ind60 =~                 
#>     x1                0.665
#>     x2                1.449
#>     x3                1.225
#>     y1                0.003
#>     y2               -0.005
#>     y3                0.003
#>     y4                0.455
#>   dem60 =~                 
#>     x1                0.020
#>     x2                0.001
#>     x3               -0.012
#>     y1                2.119
#>     y2                3.018
#>     y3                2.307
#>     y4                2.743
#> 
#> Covariances:
#>                    Estimate
#>  .x1 ~~                    
#>    .x2               -0.004
#>    .x3               -0.009
#>    .y1                0.054
#>    .y2               -0.050
#>    .y3                0.001
#>    .y4                0.018
#>  .x2 ~~                    
#>    .x3                0.006
#>    .y1               -0.002
#>    .y2                0.005
#>    .y3                0.010
#>    .y4               -0.012
#>  .x3 ~~                    
#>    .y1               -0.010
#>    .y2                0.003
#>    .y3               -0.008
#>    .y4                0.006
#>  .y1 ~~                    
#>    .y2               -0.002
#>    .y3                0.012
#>    .y4               -0.009
#>  .y2 ~~                    
#>    .y3               -0.006
#>    .y4                0.008
#>  .y3 ~~                    
#>    .y4               -0.003
#>   ind60 ~~                 
#>     dem60             0.391
#> 
#> Variances:
#>                    Estimate
#>     ind60             1.000
#>    .x1                0.069
#>    .x2                0.143
#>    .x3                0.454
#>    .y1                2.193
#>    .y2                6.148
#>    .y3                5.246
#>    .y4                2.295
#>     dem60             1.000
```

The unique covariances were all estimated close to zero. One can
approximate the “effective” number of cross-loadings and unique
covariances by:

``` r
pen_ests <- as.numeric(coef(pefa_fit2)[c(4:10, 15:35)])
sum(l0a(pen_ests))
#> [1] 1.564463
```

So out of 28 parameters penalized, only about 1.6 (or close to 2) are
effectively non-zero.

## Penalize Cross-Loadings, Unique Covariances, and Difference in Loadings and Intercepts Across Time

First, the model without cross-loadings and concurrent unique
covariances

``` r
mod3 <- "
    ind60 =~ NA * x1 + x2 + x3
    dem60 =~ NA * l1 * y1 + l2 * y2 + l3 * y3 + l4 * y4
    dem65 =~ NA * l1 * y5 + l2 * y6 + l3 * y7 + l4 * y8
    dem60 ~ ind60
    dem65 ~ ind60 + dem60
    ind60 ~~ 1 * ind60
    dem60 ~~ 1 * dem60
    dem65 ~~ NA * dem65
    ind60 ~ 0 * 1
    dem60 ~ 0 * 1
    dem65 ~ NA * 1
    x1 + x2 + x3 ~ NA * 1
    y1 ~ i1 * 1
    y2 ~ i2 * 1
    y3 ~ i3 * 1
    y4 ~ i4 * 1
    y5 ~ i1 * 1
    y6 ~ i2 * 1
    y7 ~ i3 * 1
    y8 ~ i4 * 1
    y1 ~~ y5
    y2 ~~ y6
    y3 ~~ y7
    y4 ~~ y8
"
fit3_base <- cfa(mod3, data = PoliticalDemocracy)
```

``` r
# Lavaan example of Political Democracy
mod3_un <- "
    ind60 =~ NA * x1 + x2 + x3 + y1 + y2 + y3 + y4
    dem60 =~ NA * x1 + x2 + x3 + y1 + y2 + y3 + y4
    dem65 =~ NA * y5 + y6 + y7 + y8
    dem60 ~ ind60
    dem65 ~ ind60 + dem60
    ind60 ~~ 1 * ind60
    dem60 ~~ 1 * dem60
    dem65 ~~ NA * dem65
    ind60 ~ 0 * 1
    dem60 ~ 0 * 1
    dem65 ~ NA * 1
    x1 + x2 + x3 + y1 + y2 + y3 + y4 ~ NA * 1
    y5 + y6 + y7 + y8 ~ NA * 1
    x1 ~~ x2 + x3 + y1 + y2 + y3 + y4
    x2 ~~ x3 + y1 + y2 + y3 + y4
    x3 ~~ y1 + y2 + y3 + y4
    y1 ~~ y2 + y3 + y4
    y2 ~~ y3 + y4
    y3 ~~ y4
    y1 ~~ y5
    y2 ~~ y6
    y3 ~~ y7
    y4 ~~ y8
"
fit3 <- cfa(
    mod3_un,
    data = PoliticalDemocracy,
    do.fit = FALSE,
    start = fit3_base
)
```

``` r
pt3 <- parTable(fit3)
# Provide better starting values
pt3$start[c(4:10, 35:55)] <- 0
fit3_2 <- lavaan::cfa(
    pt3,
    data = PoliticalDemocracy,
    do.fit = FALSE
)
```

Parameter IDs:

- (Concurrent) Cross-loadings: 4 to 10
- (Concurrent) Unique covariances: 35 to 55
- Loadings across time: 11 to 18
- Intercepts across time: 27 to 34

``` r
pefa_fit3 <- penalized_est(
    fit3_2,
    w = .03,
    pen_par_id = c(4:10, 35:55),
    pen_diff_id = list(
        loadings = rbind(11:14, 15:18),
        intercepts = rbind(27:30, 31:34)
    )
)
#> Warning in trans(x): NaNs produced
summary(pefa_fit3, standardized = TRUE)
#> lavaan 0.6-20 ended normally after 219 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        70
#> 
#>   Number of observations                            75
#> 
#> 
#> Parameter Estimates:
#> 
#> 
#> Latent Variables:
#>                    Estimate   Std.lv  Std.all
#>   ind60 =~                                   
#>     x1                0.753    0.753    1.038
#>     x2                1.660    1.660    1.108
#>     x3                1.421    1.421    1.020
#>     y1               -0.001   -0.001   -0.000
#>     y2               -0.001   -0.001   -0.000
#>     y3                0.479    0.479    0.141
#>     y4                0.623    0.623    0.191
#>   dem60 =~                                   
#>     x1                0.574    0.702    0.967
#>     x2                1.205    1.473    0.983
#>     x3                0.983    1.202    0.863
#>     y1                1.772    2.167    0.836
#>     y2                2.446    2.992    0.763
#>     y3                2.205    2.697    0.796
#>     y4                2.572    3.145    0.965
#>   dem65 =~                                   
#>     y5                1.768    1.992    0.772
#>     y6                2.439    2.747    0.801
#>     y7                2.250    2.535    0.798
#>     y8                2.535    2.856    0.875
#> 
#> Regressions:
#>                    Estimate   Std.lv  Std.all
#>   dem60 ~                                    
#>     ind60            -0.704   -0.576   -0.576
#>   dem65 ~                                    
#>     ind60             0.316    0.281    0.281
#>     dem60             1.005    1.091    1.091
#> 
#> Covariances:
#>                    Estimate   Std.lv  Std.all
#>  .x1 ~~                                      
#>    .x2                0.001    0.001    0.008
#>    .x3               -0.007   -0.007   -0.037
#>    .y1                0.024    0.024    0.062
#>    .y2               -0.040   -0.040   -0.057
#>    .y3                0.006    0.006    0.010
#>    .y4                0.015    0.015    0.034
#>  .x2 ~~                                      
#>    .x3                0.001    0.001    0.004
#>    .y1               -0.009   -0.009   -0.016
#>    .y2                0.009    0.009    0.009
#>    .y3                0.006    0.006    0.007
#>    .y4               -0.010   -0.010   -0.017
#>  .x3 ~~                                      
#>    .y1               -0.006   -0.006   -0.007
#>    .y2                0.002    0.002    0.001
#>    .y3               -0.010   -0.010   -0.007
#>    .y4                0.009    0.009    0.008
#>  .y1 ~~                                      
#>    .y2               -0.008   -0.008   -0.002
#>    .y3                0.011    0.011    0.003
#>    .y4               -0.009   -0.009   -0.004
#>  .y2 ~~                                      
#>    .y3               -0.004   -0.004   -0.001
#>    .y4                0.010    0.010    0.002
#>  .y3 ~~                                      
#>    .y4                0.000    0.000    0.000
#>  .y1 ~~                                      
#>    .y5                0.902    0.902    0.387
#>  .y2 ~~                                      
#>    .y6                1.594    1.594    0.307
#>  .y3 ~~                                      
#>    .y7                1.257    1.257    0.281
#>  .y4 ~~                                      
#>    .y8                0.175    0.175    0.069
#> 
#> Intercepts:
#>                    Estimate   Std.lv  Std.all
#>     ind60             0.000    0.000    0.000
#>    .dem60             0.000    0.000    0.000
#>    .dem65            -0.235   -0.208   -0.208
#>    .x1                5.072    5.072    6.990
#>    .x2                4.820    4.820    3.216
#>    .x3                3.584    3.584    2.572
#>    .y1                5.458    5.458    2.106
#>    .y2                3.764    3.764    0.960
#>    .y3                6.628    6.628    1.956
#>    .y4                4.503    4.503    1.381
#>    .y5                5.464    5.464    2.118
#>    .y6                3.748    3.748    1.093
#>    .y7                6.633    6.633    2.089
#>    .y8                4.510    4.510    1.381
#> 
#> Variances:
#>                    Estimate   Std.lv  Std.all
#>     ind60             1.000    1.000    1.000
#>    .dem60             1.000    0.668    0.668
#>    .dem65             0.106    0.084    0.084
#>    .x1                0.076    0.076    0.144
#>    .x2                0.137    0.137    0.061
#>    .x3                0.445    0.445    0.229
#>    .y1                2.017    2.017    0.300
#>    .y2                6.413    6.413    0.417
#>    .y3                5.473    5.473    0.476
#>    .y4                2.601    2.601    0.245
#>    .y5                2.690    2.690    0.404
#>    .y6                4.211    4.211    0.358
#>    .y7                3.653    3.653    0.362
#>    .y8                2.503    2.503    0.235
```

We can again compute the “effective” number of cross-loadings and unique
covariances that are non-zero:

``` r
pen_ests2 <- as.numeric(coef(pefa_fit3)[c(4:10, 35:55)])
sum(l0a(pen_ests2))
#> [1] 5.1989
```

And the “effective” number of loadings and intercepts that differ across
time:

``` r
ld_ests <- as.numeric(coef(pefa_fit3)[11:18])
int_ests <- as.numeric(coef(pefa_fit3)[27:34])
ld_mat <- matrix(ld_ests, nrow = 2, byrow = TRUE)
int_mat <- matrix(int_ests, nrow = 2, byrow = TRUE)
composite_pair_loss(ld_mat, fun = l0a) +
    composite_pair_loss(int_mat, fun = l0a)
#> [1] 0.3263673
```
