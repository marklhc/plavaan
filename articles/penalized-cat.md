# Penalized Estimation with Ordinal Data and Multiple Groups

``` r

library(plavaan)
library(lavaan)
#> This is lavaan 0.7-2
#> lavaan is FREE software! Please report any bugs.
data(HolzingerSwineford1939)
```

## Prepare Ordinal Data

First, we convert the continuous variables to ordinal with 3 categories:

``` r

# Select the 9 cognitive test variables
hs_data <- HolzingerSwineford1939[, c("school", "x1", "x2", "x3", "x4", 
                                       "x5", "x6", "x7", "x8", "x9")]

# Convert to ordinal with 3 points
for (i in 2:10) {
    hs_data[[i]] <- cut(
        hs_data[[i]], 
        breaks = 3, 
        labels = FALSE,
        include.lowest = TRUE
    )
    hs_data[[i]] <- ordered(hs_data[[i]])
}

head(hs_data)
#>    school x1 x2 x3 x4 x5 x6 x7 x8 x9
#> 1 Pasteur  2  3  1  2  3  1  2  2  2
#> 2 Pasteur  2  2  2  1  1  1  2  2  3
#> 3 Pasteur  2  2  2  1  1  1  1  1  1
#> 4 Pasteur  2  3  2  2  2  2  1  1  1
#> 5 Pasteur  2  2  1  2  2  2  2  2  2
#> 6 Pasteur  2  2  2  1  1  1  2  2  3
```

## Multiple-Group CFA Model

Now fit the model across the two schools:

``` r

mod_base <- "
  visual =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed =~ x7 + x8 + x9
  x1 ~~ 1 * x1
  x2 ~~ 1 * x2
  x3 ~~ 1 * x3
  x4 ~~ 1 * x4
  x5 ~~ 1 * x5
  x6 ~~ 1 * x6
  x7 ~~ 1 * x7
  x8 ~~ 1 * x8
  x9 ~~ 1 * x9
"
fit_mg <- cfa(mod_base, data = hs_data, ordered = TRUE, std.lv = TRUE,
              parameterization = "theta", group = "school")
summary(fit_mg, fit.measures = TRUE)
#> lavaan 0.7-2 ended normally after 179 iterations
#> 
#>   Estimator                                       DWLS
#>   Optimization method                           NLMINB
#>   Number of model parameters                        60
#> 
#>   Number of observations per group:                   
#>     Pasteur                                        156
#>     Grant-White                                    145
#> 
#> Model Test User Model:
#>                                               Standard      Scaled
#>   Test Statistic                                71.775      99.182
#>   Degrees of freedom                                48          48
#>   P-value (Unknown)                                 NA       0.000
#>   Scaling correction factor                                  0.798
#>   Shift parameter                                            9.213
#>     simple second-order correction                                
#>   Test statistic for each group:
#>     Pasteur                                     56.128      56.128
#>     Grant-White                                 43.055      43.055
#> 
#> Model Test Baseline Model:
#> 
#>   Test statistic                              1706.441    1175.245
#>   Degrees of freedom                                72          72
#>   P-value                                           NA       0.000
#>   Scaling correction factor                                  1.481
#> 
#> User Model versus Baseline Model:
#> 
#>   Comparative Fit Index (CFI)                    0.985       0.954
#>   Tucker-Lewis Index (TLI)                       0.978       0.930
#>                                                                   
#>   Robust Comparative Fit Index (CFI)                         0.863
#>   Robust Tucker-Lewis Index (TLI)                            0.794
#> 
#> Root Mean Square Error of Approximation:
#> 
#>   RMSEA                                          0.058       0.084
#>   90 Percent confidence interval - lower         0.026       0.061
#>   90 Percent confidence interval - upper         0.084       0.108
#>   P-value H_0: RMSEA <= 0.050                    0.308       0.011
#>   P-value H_0: RMSEA >= 0.080                    0.084       0.641
#>                                                                   
#>   Robust RMSEA                                               0.143
#>   90 Percent confidence interval - lower                     0.098
#>   90 Percent confidence interval - upper                     0.186
#>   P-value H_0: Robust RMSEA <= 0.050                         0.001
#>   P-value H_0: Robust RMSEA >= 0.080                         0.987
#> 
#> Standardized Root Mean Square Residual:
#> 
#>   SRMR                                           0.087       0.087
#> 
#> Goodness of Fit Index:
#> 
#>   Goodness of Fit Index (GFI)                    0.983            
#>   90 Percent confidence interval - lower         0.964            
#>   90 Percent confidence interval - upper         0.996            
#>                                                                   
#>   Robust GFI                                                 0.902
#>   90 Percent confidence interval - lower                     0.844
#>   90 Percent confidence interval - upper                     0.951
#> 
#> Parameter Estimates:
#> 
#>   Parameterization                               Theta
#>   Standard errors                           Robust.sem
#>   Information                                 Expected
#>   Information saturated (h1) model        Unstructured
#> 
#> 
#> Group 1 [Pasteur]:
#> 
#> Latent Variables:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   visual =~                                           
#>     x1                2.581    2.738    0.943    0.346
#>     x2                0.536    0.162    3.315    0.001
#>     x3                0.683    0.177    3.866    0.000
#>   textual =~                                          
#>     x4                1.504    0.340    4.424    0.000
#>     x5                2.586    1.010    2.560    0.010
#>     x6                1.766    0.443    3.988    0.000
#>   speed =~                                            
#>     x7                0.727    0.243    2.992    0.003
#>     x8                0.938    0.353    2.654    0.008
#>     x9                0.792    0.280    2.830    0.005
#> 
#> Covariances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   visual ~~                                           
#>     textual           0.426    0.099    4.292    0.000
#>     speed             0.187    0.113    1.654    0.098
#>   textual ~~                                          
#>     speed             0.293    0.107    2.738    0.006
#> 
#> Thresholds:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>     x1|t1            -4.074    3.773   -1.080    0.280
#>     x1|t2             2.407    2.214    1.087    0.277
#>     x2|t1            -1.523    0.181   -8.415    0.000
#>     x2|t2             0.909    0.138    6.586    0.000
#>     x3|t1            -0.699    0.139   -5.011    0.000
#>     x3|t2             0.522    0.130    3.998    0.000
#>     x4|t1            -0.940    0.231   -4.077    0.000
#>     x4|t2             2.106    0.376    5.594    0.000
#>     x5|t1            -1.444    0.555   -2.601    0.009
#>     x5|t2             2.042    0.759    2.690    0.007
#>     x6|t1             0.874    0.287    3.045    0.002
#>     x6|t2             3.956    0.801    4.941    0.000
#>     x7|t1            -1.330    0.211   -6.309    0.000
#>     x7|t2             1.018    0.179    5.683    0.000
#>     x8|t1            -0.110    0.139   -0.793    0.428
#>     x8|t2             2.538    0.502    5.054    0.000
#>     x9|t1            -0.505    0.146   -3.465    0.001
#>     x9|t2             2.164    0.352    6.155    0.000
#> 
#> Variances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>    .x1                1.000                           
#>    .x2                1.000                           
#>    .x3                1.000                           
#>    .x4                1.000                           
#>    .x5                1.000                           
#>    .x6                1.000                           
#>    .x7                1.000                           
#>    .x8                1.000                           
#>    .x9                1.000                           
#>     visual            1.000                           
#>     textual           1.000                           
#>     speed             1.000                           
#> 
#> Scales y*:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>     x1                0.361                           
#>     x2                0.881                           
#>     x3                0.826                           
#>     x4                0.554                           
#>     x5                0.361                           
#>     x6                0.493                           
#>     x7                0.809                           
#>     x8                0.730                           
#>     x9                0.784                           
#> 
#> 
#> Group 2 [Grant-White]:
#> 
#> Latent Variables:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   visual =~                                           
#>     x1                1.149    0.383    2.997    0.003
#>     x2                0.575    0.151    3.812    0.000
#>     x3                0.802    0.184    4.370    0.000
#>   textual =~                                          
#>     x4                1.629    0.358    4.551    0.000
#>     x5                2.303    0.698    3.299    0.001
#>     x6                1.166    0.207    5.632    0.000
#>   speed =~                                            
#>     x7                0.977    0.218    4.492    0.000
#>     x8                0.995    0.231    4.314    0.000
#>     x9                3.035    2.806    1.082    0.279
#> 
#> Covariances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   visual ~~                                           
#>     textual           0.596    0.086    6.964    0.000
#>     speed             0.553    0.108    5.100    0.000
#>   textual ~~                                          
#>     speed             0.423    0.091    4.661    0.000
#> 
#> Thresholds:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>     x1|t1            -1.923    0.406   -4.736    0.000
#>     x1|t2             1.245    0.293    4.242    0.000
#>     x2|t1            -2.212    0.258   -8.564    0.000
#>     x2|t2             0.784    0.136    5.769    0.000
#>     x3|t1            -0.167    0.135   -1.239    0.215
#>     x3|t2             1.111    0.170    6.522    0.000
#>     x4|t1            -2.272    0.413   -5.502    0.000
#>     x4|t2             1.471    0.312    4.716    0.000
#>     x5|t1            -3.168    0.850   -3.726    0.000
#>     x5|t2             0.817    0.331    2.465    0.014
#>     x6|t1            -0.093    0.161   -0.579    0.563
#>     x6|t2             1.939    0.274    7.063    0.000
#>     x7|t1            -0.832    0.177   -4.711    0.000
#>     x7|t2             1.820    0.263    6.908    0.000
#>     x8|t1            -0.257    0.151   -1.707    0.088
#>     x8|t2             2.705    0.421    6.430    0.000
#>     x9|t1            -1.458    1.266   -1.151    0.250
#>     x9|t2             6.519    5.519    1.181    0.238
#> 
#> Variances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>    .x1                1.000                           
#>    .x2                1.000                           
#>    .x3                1.000                           
#>    .x4                1.000                           
#>    .x5                1.000                           
#>    .x6                1.000                           
#>    .x7                1.000                           
#>    .x8                1.000                           
#>    .x9                1.000                           
#>     visual            1.000                           
#>     textual           1.000                           
#>     speed             1.000                           
#> 
#> Scales y*:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>     x1                0.656                           
#>     x2                0.867                           
#>     x3                0.780                           
#>     x4                0.523                           
#>     x5                0.398                           
#>     x6                0.651                           
#>     x7                0.715                           
#>     x8                0.709                           
#>     x9                0.313
```

One should also consider the magnitude of the objective function with
the DWLS estimator, which is scaled differently than ML-based functions
and is generally smaller based on experience.

``` r

fit_mg@optim$fx
#> [1] 0.1192268
```

### Strict and partial invariance models

``` r

# Strict invariance: constrain loadings, thresholds, and residual variances
fit_strict <- cfa(mod_base, data = hs_data, ordered = TRUE, std.lv = TRUE,
                  parameterization = "theta", group = "school",
                  group.equal = c("loadings", "thresholds", "residuals"))

# Score test
lavTestScore(fit_strict)
#> $test
#> 
#> total score test:
#> 
#>             test     X2     df p.value
#> 1          score 18.231 27.000   0.896
#> 2   score.scaled 17.168 27.000   0.927
#> 3 score.adjusted 12.843 20.199   0.891
#> 4   score.robust 23.327 27.000   0.667
#> 
#> $uni
#> 
#> univariate score tests:
#> 
#>      lhs op   rhs    X2 df p.value X2.scaled p.value.scaled
#> 1   .p1. == .p64. 0.345  1   0.557     0.305          0.580
#> 2   .p2. == .p65. 0.730  1   0.393     0.443          0.506
#> 3   .p3. == .p66. 0.027  1   0.869     0.024          0.877
#> 4   .p4. == .p67. 0.110  1   0.740     0.177          0.674
#> 5   .p5. == .p68. 0.310  1   0.578     0.599          0.439
#> 6   .p6. == .p69. 0.049  1   0.824     0.059          0.808
#> 7   .p7. == .p70. 0.995  1   0.319     1.150          0.284
#> 8   .p8. == .p71. 0.003  1   0.960     0.003          0.960
#> 9   .p9. == .p72. 0.769  1   0.380     0.694          0.405
#> 10 .p19. == .p82. 1.308  1   0.253     1.491          0.222
#> 11 .p20. == .p83. 1.308  1   0.253     1.491          0.222
#> 12 .p21. == .p84. 2.439  1   0.118     2.954          0.086
#> 13 .p22. == .p85. 2.439  1   0.118     2.954          0.086
#> 14 .p23. == .p86. 0.000  1   0.986     0.001          0.981
#> 15 .p24. == .p87. 0.000  1   0.986     0.001          0.981
#> 16 .p25. == .p88. 0.082  1   0.774     0.149          0.700
#> 17 .p26. == .p89. 0.082  1   0.774     0.149          0.700
#> 18 .p27. == .p90. 0.415  1   0.519     0.750          0.386
#> 19 .p28. == .p91. 0.415  1   0.519     0.750          0.386
#> 20 .p29. == .p92. 1.428  1   0.232     2.014          0.156
#> 21 .p30. == .p93. 1.428  1   0.232     2.014          0.156
#> 22 .p31. == .p94. 0.158  1   0.691     0.219          0.640
#> 23 .p32. == .p95. 0.158  1   0.691     0.219          0.640
#> 24 .p33. == .p96. 0.802  1   0.371     0.965          0.326
#> 25 .p34. == .p97. 0.802  1   0.371     0.965          0.326
#> 26 .p35. == .p98. 4.157  1   0.041     4.138          0.042
#> 27 .p36. == .p99. 4.157  1   0.041     4.138          0.042
```

Only item 9 showed non-invariant thresholds, based on the score test.

## Penalized Multiple-Group Model

We’ll penalize differences in loadings and thresholds across groups.
First, set up an over-specified (unidentified) model with the latent
mean and variance only identified in the first group, without fitting
(`do.fit = FALSE`):

``` r

mod_un <- "
  visual =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed =~ x7 + x8 + x9
  visual ~~ c(1, NA) * visual
  textual ~~ c(1, NA) * textual
  speed ~~ c(1, NA) * speed
  visual ~ c(0, NA) * 1
  textual ~ c(0, NA) * 1
  speed ~ c(0, NA) * 1
  x1 ~~ 1 * x1
  x2 ~~ 1 * x2
  x3 ~~ 1 * x3
  x4 ~~ 1 * x4
  x5 ~~ 1 * x5
  x6 ~~ 1 * x6
  x7 ~~ 1 * x7
  x8 ~~ 1 * x8
  x9 ~~ 1 * x9
"
fit_mg_nofit <- cfa(mod_un, data = hs_data, ordered = TRUE, std.lv = TRUE,
                    auto.fix.first = FALSE,
                    parameterization = "theta", group = "school", do.fit = FALSE)
```

Examine the parameter table to identify loadings and thresholds:

``` r

pt <- parTable(fit_mg_nofit)
# Show loadings
pt[pt$op == "=~", c("lhs", "op", "rhs", "group", "free")]
#>        lhs op rhs group free
#> 1   visual =~  x1     1    1
#> 2   visual =~  x2     1    2
#> 3   visual =~  x3     1    3
#> 4  textual =~  x4     1    4
#> 5  textual =~  x5     1    5
#> 6  textual =~  x6     1    6
#> 7    speed =~  x7     1    7
#> 8    speed =~  x8     1    8
#> 9    speed =~  x9     1    9
#> 64  visual =~  x1     2   31
#> 65  visual =~  x2     2   32
#> 66  visual =~  x3     2   33
#> 67 textual =~  x4     2   34
#> 68 textual =~  x5     2   35
#> 69 textual =~  x6     2   36
#> 70   speed =~  x7     2   37
#> 71   speed =~  x8     2   38
#> 72   speed =~  x9     2   39
```

``` r

# Show thresholds
pt[pt$op == "|", c("lhs", "op", "rhs", "group", "free")]
#>     lhs op rhs group free
#> 25   x1  |  t1     1   10
#> 26   x1  |  t2     1   11
#> 27   x2  |  t1     1   12
#> 28   x2  |  t2     1   13
#> 29   x3  |  t1     1   14
#> 30   x3  |  t2     1   15
#> 31   x4  |  t1     1   16
#> 32   x4  |  t2     1   17
#> 33   x5  |  t1     1   18
#> 34   x5  |  t2     1   19
#> 35   x6  |  t1     1   20
#> 36   x6  |  t2     1   21
#> 37   x7  |  t1     1   22
#> 38   x7  |  t2     1   23
#> 39   x8  |  t1     1   24
#> 40   x8  |  t2     1   25
#> 41   x9  |  t1     1   26
#> 42   x9  |  t2     1   27
#> 88   x1  |  t1     2   46
#> 89   x1  |  t2     2   47
#> 90   x2  |  t1     2   48
#> 91   x2  |  t2     2   49
#> 92   x3  |  t1     2   50
#> 93   x3  |  t2     2   51
#> 94   x4  |  t1     2   52
#> 95   x4  |  t2     2   53
#> 96   x5  |  t1     2   54
#> 97   x5  |  t2     2   55
#> 98   x6  |  t1     2   56
#> 99   x6  |  t2     2   57
#> 100  x7  |  t1     2   58
#> 101  x7  |  t2     2   59
#> 102  x8  |  t1     2   60
#> 103  x8  |  t2     2   61
#> 104  x9  |  t1     2   62
#> 105  x9  |  t2     2   63
```

Identify parameter IDs for loadings and thresholds in each group:

``` r

# Loadings: group 1 (Pasteur) and group 2 (Grant-White)
load_g1 <- pt$free[pt$op == "=~" & pt$group == 1 & pt$free > 0]
load_g2 <- pt$free[pt$op == "=~" & pt$group == 2 & pt$free > 0]

# Thresholds: group 1 and group 2
thresh_g1 <- pt$free[pt$op == "|" & pt$group == 1 & pt$free > 0]
thresh_g2 <- pt$free[pt$op == "|" & pt$group == 2 & pt$free > 0]

print(list(
    loadings_g1 = load_g1,
    loadings_g2 = load_g2,
    thresholds_g1 = thresh_g1,
    thresholds_g2 = thresh_g2
))
#> $loadings_g1
#> [1] 1 2 3 4 5 6 7 8 9
#> 
#> $loadings_g2
#> [1] 31 32 33 34 35 36 37 38 39
#> 
#> $thresholds_g1
#>  [1] 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
#> 
#> $thresholds_g2
#>  [1] 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63
```

Fit the penalized model with penalties on differences in loadings and
thresholds:

``` r

fit_pen_mg <- penalized_est(
    fit_mg_nofit,
    w = 0.03,
    pen_diff_id = list(
        loadings = rbind(load_g1, load_g2),
        thresholds = rbind(thresh_g1, thresh_g2)
    )
)
summary(fit_pen_mg)
#> lavaan 0.7-2 ended normally after 207 iterations
#> 
#>   Estimator                                       DWLS
#>   Optimization method                           NLMINB
#>   Number of model parameters                        66
#> 
#>   Number of observations per group:                   
#>     Pasteur                                        156
#>     Grant-White                                    145
#> 
#> 
#> Parameter Estimates:
#> 
#>   Parameterization                               Theta
#> 
#> 
#> Group 1 [Pasteur]:
#> 
#> Latent Variables:
#>                    Estimate
#>   visual =~                
#>     x1                1.332
#>     x2                0.540
#>     x3                0.785
#>   textual =~               
#>     x4                1.716
#>     x5                2.640
#>     x6                1.502
#>   speed =~                 
#>     x7                0.787
#>     x8                0.809
#>     x9                1.189
#> 
#> Covariances:
#>                    Estimate
#>   visual ~~                
#>     textual           0.445
#>     speed             0.205
#>   textual ~~               
#>     speed             0.275
#> 
#> Intercepts:
#>                    Estimate
#>     visual            0.000
#>     textual           0.000
#>     speed             0.000
#> 
#> Thresholds:
#>                    Estimate
#>     x1|t1            -2.520
#>     x1|t2             1.248
#>     x2|t1            -1.534
#>     x2|t2             0.771
#>     x3|t1            -0.571
#>     x3|t2             0.691
#>     x4|t1            -1.063
#>     x4|t2             2.366
#>     x5|t1            -1.511
#>     x5|t2             2.256
#>     x6|t1             0.788
#>     x6|t2             3.112
#>     x7|t1            -1.158
#>     x7|t2             1.063
#>     x8|t1            -0.217
#>     x8|t2             2.475
#>     x9|t1            -0.763
#>     x9|t2             2.883
#> 
#> Variances:
#>                    Estimate
#>     visual            1.000
#>     textual           1.000
#>     speed             1.000
#>    .x1                1.000
#>    .x2                1.000
#>    .x3                1.000
#>    .x4                1.000
#>    .x5                1.000
#>    .x6                1.000
#>    .x7                1.000
#>    .x8                1.000
#>    .x9                1.000
#> 
#> 
#> Group 2 [Grant-White]:
#> 
#> Latent Variables:
#>                    Estimate
#>   visual =~                
#>     x1                1.331
#>     x2                0.540
#>     x3                0.785
#>   textual =~               
#>     x4                1.717
#>     x5                2.639
#>     x6                1.503
#>   speed =~                 
#>     x7                0.792
#>     x8                0.807
#>     x9                1.187
#> 
#> Covariances:
#>                    Estimate
#>   visual ~~                
#>     textual           0.519
#>     speed             0.728
#>   textual ~~               
#>     speed             0.477
#> 
#> Intercepts:
#>                    Estimate
#>     visual           -0.298
#>     textual           0.596
#>     speed            -0.139
#> 
#> Thresholds:
#>                    Estimate
#>     x1|t1            -2.520
#>     x1|t2             1.245
#>     x2|t1            -2.352
#>     x2|t2             0.767
#>     x3|t1            -0.566
#>     x3|t2             0.695
#>     x4|t1            -1.064
#>     x4|t2             2.367
#>     x5|t1            -1.511
#>     x5|t2             2.257
#>     x6|t1             0.788
#>     x6|t2             3.110
#>     x7|t1            -1.153
#>     x7|t2             1.693
#>     x8|t1            -0.221
#>     x8|t2             2.475
#>     x9|t1            -0.766
#>     x9|t2             2.884
#> 
#> Variances:
#>                    Estimate
#>     visual            1.130
#>     textual           0.746
#>     speed             1.563
#>    .x1                1.000
#>    .x2                1.000
#>    .x3                1.000
#>    .x4                1.000
#>    .x5                1.000
#>    .x6                1.000
#>    .x7                1.000
#>    .x8                1.000
#>    .x9                1.000
```

## Evaluate Invariance

Here are the estimated loadings and thresholds, and we can calculate the
effective number of parameters that differ across groups:

``` r

# Loadings
load_ests_g1 <- as.numeric(coef(fit_pen_mg)[load_g1])
load_ests_g2 <- as.numeric(coef(fit_pen_mg)[load_g2])
load_mat <- rbind(load_ests_g1, load_ests_g2)
colnames(load_mat) <- names(coef(fit_pen_mg))[load_g1]
eff_load_diff <- composite_pair_loss(load_mat, fun = l0a)

# Thresholds
thresh_ests_g1 <- as.numeric(coef(fit_pen_mg)[thresh_g1])
thresh_ests_g2 <- as.numeric(coef(fit_pen_mg)[thresh_g2])
thresh_mat <- rbind(thresh_ests_g1, thresh_ests_g2)
colnames(thresh_mat) <- names(coef(fit_pen_mg))[thresh_g1]
eff_thresh_diff <- composite_pair_loss(thresh_mat, fun = l0a)

cat("Penalized Loading Estimates:\n")
#> Penalized Loading Estimates:
print(load_mat, digits = 3)
#>              visual=~x1 visual=~x2 visual=~x3 textual=~x4 textual=~x5
#> load_ests_g1       1.33       0.54      0.785        1.72        2.64
#> load_ests_g2       1.33       0.54      0.785        1.72        2.64
#>              textual=~x6 speed=~x7 speed=~x8 speed=~x9
#> load_ests_g1         1.5     0.787     0.809      1.19
#> load_ests_g2         1.5     0.792     0.807      1.19

cat("\nPenalized Threshold Estimates:\n")
#> 
#> Penalized Threshold Estimates:
print(thresh_mat, digits = 3)
#>                x1|t1 x1|t2 x2|t1 x2|t2  x3|t1 x3|t2 x4|t1 x4|t2 x5|t1 x5|t2
#> thresh_ests_g1 -2.52  1.25 -1.53 0.771 -0.571 0.691 -1.06  2.37 -1.51  2.26
#> thresh_ests_g2 -2.52  1.24 -2.35 0.767 -0.566 0.695 -1.06  2.37 -1.51  2.26
#>                x6|t1 x6|t2 x7|t1 x7|t2  x8|t1 x8|t2  x9|t1 x9|t2
#> thresh_ests_g1 0.788  3.11 -1.16  1.06 -0.217  2.47 -0.763  2.88
#> thresh_ests_g2 0.788  3.11 -1.15  1.69 -0.221  2.48 -0.766  2.88

cat("Effective number of non-invariant loadings:", eff_load_diff, "\n")
#> Effective number of non-invariant loadings: 0.002959255
cat("Effective number of non-invariant thresholds:", eff_thresh_diff, "\n")
#> Effective number of non-invariant thresholds: 1.973607
```

The penalized estimation approach identifies which loadings and
thresholds substantively differ across groups, providing an efficient,
data-driven assessment of measurement invariance for ordinal data.
