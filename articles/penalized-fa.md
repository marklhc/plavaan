# Penalized Estimation of Cross-Loadings and Unique Covariances

``` r

library(plavaan)
library(lavaan)
#> This is lavaan 0.7-2
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
    pen_par_id = 4:10,
    test = "Chisq"
)
summary(pefa_fit)
#> Penalized fit (w = 0.03, eps = 0.01, penalty = l0a): effective npar = 16.03, effective df = 11.97 (nominal df = 6).
#> Fit evaluation for penalized fits is experimental; interpret the chi-square test and fit indices with caution.
#> lavaan 0.7-2 ended normally after 128 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        16
#> 
#>   Number of observations                            75
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                                19.803
#>   Degrees of freedom                            11.966
#>   P-value (Chi-square)                           0.070
#> 
#> Parameter Estimates:
#> 
#> 
#> Latent Variables:
#>                    Estimate
#>   ind60 =~                 
#>     x1                0.657
#>     x2                1.458
#>     x3                1.221
#>     y1                0.002
#>     y2               -0.005
#>     y3                0.003
#>     y4                0.423
#>   dem60 =~                 
#>     x1                0.027
#>     x2               -0.001
#>     x3               -0.011
#>     y1                2.129
#>     y2                2.980
#>     y3                2.314
#>     y4                2.752
#> 
#> Covariances:
#>                    Estimate
#>   ind60 ~~                 
#>     dem60             0.394
#> 
#> Variances:
#>                    Estimate
#>     ind60             1.000
#>    .x1                0.081
#>    .x2                0.120
#>    .x3                0.463
#>    .y1                2.228
#>    .y2                6.458
#>    .y3                5.229
#>    .y4                2.341
#>     dem60             1.000
```

Fit indices can be obtained directly from the penalized fit. Fit
evaluation is experimental and disabled by default, so the fit above was
created with `test = "Chisq"`; with the default `test = "none"`,
[`fitmeasures()`](https://rdrr.io/pkg/lavaan/man/fitMeasures.html) is
unavailable and [`summary()`](https://rdrr.io/r/base/summary.html) shows
no chi-square test. When enabled,
[`fitmeasures()`](https://rdrr.io/pkg/lavaan/man/fitMeasures.html)
refits the model with all parameters frozen at the penalized estimates
and reports the indices at the effective degrees of freedom (11.97 here,
versus the nominal 6), and an experimental notice is shown. For the
under-identified models later in this vignette (negative nominal df),
the effective df is the meaningful value.

``` r

fitmeasures(pefa_fit, c("chisq", "df", "cfi", "rmsea"))
#> Fit evaluation for penalized fits is experimental; interpret fit indices with caution.
#>                  npar                  fmin                 chisq 
#>                16.000                 0.132                19.803 
#>                    df                pvalue        baseline.chisq 
#>                11.966                 0.070               406.880 
#>           baseline.df       baseline.pvalue                   cfi 
#>                21.000                 0.000                 0.980 
#>                   tli                  nnfi                   rfi 
#>                 0.964                 0.964                 0.915 
#>                   nfi                  pnfi                   ifi 
#>                 0.951                 0.542                 0.980 
#>                   rni                  logl     unrestricted.logl 
#>                 0.980              -936.115              -926.214 
#>                   aic                   bic                ntotal 
#>              1904.298              1941.458                75.000 
#>                  bic2                 rmsea        rmsea.ci.lower 
#>              1930.344                 0.093                 0.000 
#>        rmsea.ci.upper        rmsea.ci.level          rmsea.pvalue 
#>                 0.164                 0.900                 0.159 
#>        rmsea.close.h0 rmsea.notclose.pvalue     rmsea.notclose.h0 
#>                 0.050                 0.661                 0.080 
#>                   rmr            rmr_nomean                  srmr 
#>                 0.317                 0.317                 0.038 
#>          srmr_bentler   srmr_bentler_nomean                  crmr 
#>                 0.038                 0.038                 0.044 
#>           crmr_nomean            srmr_mplus     srmr_mplus_nomean 
#>                 0.044                 0.038                 0.038 
#>                   gfi          gfi.ci.lower          gfi.ci.upper 
#>                 0.972                 0.917                 1.000 
#>          gfi.ci.level                 cn_05                 cn_01 
#>                 0.900                80.460               100.102 
#>            gfi_lisrel           agfi_lisrel                  pgfi 
#>                 0.931                 0.838                 0.398 
#>                   mfi                  ecvi 
#>                 0.949                 0.691
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
#> Penalized fit (w = 0.03, eps = 0.01, penalty = l0a): effective npar = 16.56, effective df = 11.44 (nominal df = -15).
#> lavaan 0.7-2 ended normally after 171 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        17
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

The unique covariances were all estimated close to zero. The effective
number of cross-loadings and unique covariances can be reported with
[`effective_df()`](https://marklhc.github.io/plavaan/reference/effective_df.md):

``` r

effective_df(pefa_fit2)
#>                npar npar_effective df_saved
#> direct penalty   28       1.564467 26.43553
#> TOTAL            43      16.564467 26.43553
#> 
#> n_stats (sample moments):  28
#> nominal model df:  -15 (negative: the nominal model is under-identified; the effective df is the meaningful quantity)
#> effective model df:  11.44
#> penalty:  l0a (w = 0.03, eps = 0.01)
```

So out of 28 penalized parameters, only about 1.6 (or close to 2) are
effectively non-zero. The table also shows the effective model degrees
of freedom (11.44), where the nominal df (-15) is negative because the
model is under-identified.

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
summary(pefa_fit3, standardized = TRUE)
#> Penalized fit (w = 0.03, eps = 0.01, penalty = l0a): effective npar = 35.4, effective df = 41.6 (nominal df = 7).
#> lavaan 0.7-2 ended normally after 189 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        35
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
#>     x1                0.664    0.664    0.915
#>     x2                1.451    1.451    0.968
#>     x3                1.226    1.226    0.879
#>     y1               -0.002   -0.002   -0.001
#>     y2               -0.000   -0.000   -0.000
#>     y3                0.001    0.001    0.000
#>     y4                0.448    0.448    0.134
#>   dem60 =~                                   
#>     x1                0.022    0.024    0.033
#>     x2                0.002    0.002    0.001
#>     x3               -0.012   -0.014   -0.010
#>     y1                1.944    2.111    0.829
#>     y2                2.620    2.846    0.743
#>     y3                2.306    2.505    0.734
#>     y4                2.527    2.744    0.819
#>   dem65 =~                                   
#>     y5                1.936    2.079    0.789
#>     y6                2.615    2.809    0.805
#>     y7                2.316    2.488    0.792
#>     y8                2.530    2.717    0.858
#> 
#> Regressions:
#>                    Estimate   Std.lv  Std.all
#>   dem60 ~                                    
#>     ind60             0.424    0.390    0.390
#>   dem65 ~                                    
#>     ind60             0.232    0.216    0.216
#>     dem60             0.844    0.854    0.854
#> 
#> Covariances:
#>                    Estimate   Std.lv  Std.all
#>  .x1 ~~                                      
#>    .x2               -0.002   -0.002   -0.025
#>    .x3               -0.008   -0.008   -0.046
#>    .y1                0.027    0.027    0.069
#>    .y2               -0.039   -0.039   -0.056
#>    .y3                0.009    0.009    0.014
#>    .y4                0.016    0.016    0.037
#>  .x2 ~~                                      
#>    .x3                0.004    0.004    0.016
#>    .y1               -0.010   -0.010   -0.018
#>    .y2                0.008    0.008    0.008
#>    .y3                0.011    0.011    0.012
#>    .y4               -0.010   -0.010   -0.017
#>  .x3 ~~                                      
#>    .y1               -0.009   -0.009   -0.010
#>    .y2                0.000    0.000    0.000
#>    .y3               -0.009   -0.009   -0.005
#>    .y4                0.008    0.008    0.007
#>  .y1 ~~                                      
#>    .y2               -0.006   -0.006   -0.002
#>    .y3                0.010    0.010    0.003
#>    .y4               -0.010   -0.010   -0.004
#>  .y2 ~~                                      
#>    .y3               -0.005   -0.005   -0.001
#>    .y4                0.010    0.010    0.002
#>  .y3 ~~                                      
#>    .y4               -0.000   -0.000   -0.000
#>  .y1 ~~                                      
#>    .y5                0.829    0.829    0.360
#>  .y2 ~~                                      
#>    .y6                1.683    1.683    0.317
#>  .y3 ~~                                      
#>    .y7                1.131    1.131    0.255
#>  .y4 ~~                                      
#>    .y8                0.218    0.218    0.084
#> 
#> Intercepts:
#>                    Estimate   Std.lv  Std.all
#>     ind60             0.000    0.000    0.000
#>    .dem60             0.000    0.000    0.000
#>    .dem65            -0.232   -0.216   -0.216
#>    .x1                5.070    5.070    6.985
#>    .x2                4.814    4.814    3.213
#>    .x3                3.577    3.577    2.566
#>    .y1                5.525    5.525    2.169
#>    .y2                3.849    3.849    1.005
#>    .y3                6.685    6.685    1.960
#>    .y4                4.558    4.558    1.361
#>    .y5                5.532    5.532    2.100
#>    .y6                3.834    3.834    1.099
#>    .y7                6.690    6.690    2.131
#>    .y8                4.564    4.564    1.440
#> 
#> Variances:
#>                    Estimate   Std.lv  Std.all
#>     ind60             1.000    1.000    1.000
#>    .dem60             1.000    0.848    0.848
#>    .dem65             0.093    0.081    0.081
#>    .x1                0.072    0.072    0.137
#>    .x2                0.138    0.138    0.061
#>    .x3                0.453    0.453    0.233
#>    .y1                2.031    2.031    0.313
#>    .y2                6.570    6.570    0.448
#>    .y3                5.358    5.358    0.461
#>    .y4                2.524    2.524    0.225
#>    .y5                2.616    2.616    0.377
#>    .y6                4.289    4.289    0.352
#>    .y7                3.666    3.666    0.372
#>    .y8                2.655    2.655    0.264
```

[`effective_df()`](https://marklhc.github.io/plavaan/reference/effective_df.md)
reports the effective number of cross-loadings and unique covariances
that are non-zero, and the effective number of loadings and intercepts
that differ across time:

``` r

effective_df(pefa_fit3)
#>                npar npar_effective  df_saved
#> direct penalty   28       1.343754 26.656246
#> loadings          8       4.020045  3.979955
#> intercepts        8       4.032270  3.967730
#> TOTAL            70      35.396069 34.603931
#> 
#> n_stats (sample moments):  77
#> nominal model df:  7
#> effective model df:  41.6
#> penalty:  l0a (w = 0.03, eps = 0.01)
```

For a difference-penalty block, `npar_effective` is the number of
columns (one shared invariance baseline per parameter) plus the
effective number of non-invariant values. Here the `loadings` and
`intercepts` rows (4.02 and 4.03, versus 4 baseline values each)
indicate that the loadings and intercepts are effectively invariant
across time, while only about 1.3 of the 28 cross-loadings and unique
covariances are effectively non-zero.
