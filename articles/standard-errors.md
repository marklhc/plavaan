# Standard Errors

This vignette demonstrates the calculation of standard errors for
penalized estimates, using the sandwich estimator approach as in the
“MLR” estimator of `lavaan`. The standard errors are obtained as the
square roots of the diagonal elements of the sandwich variance
estimator:

``` math
\mathrm{B}_\text{bread}^{-1} \; \mathrm{B}_\text{meat} \; \mathrm{B}_\text{bread}^{-1},
```

where $`\mathrm{B}_\text{bread}`$ is the observed information matrix
(i.e., the Hessian of the penalized objective function) and
$`\mathrm{B}_\text{meat}`$ is the first-order information matrix
(obtained using `lavInspect(..., information.first.order)`).

``` r

library(lavaan)
#> This is lavaan 0.7-2
#> lavaan is FREE software! Please report any bugs.
library(plavaan)
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
fit0 <- cfa(mod0, data = PoliticalDemocracy, std.lv = TRUE, estimator = "MLR")
```

Penalized

``` r

mod <- "
  ind60 =~ x1 + x2 + x3 + y1 + y2 + y3 + y4
  dem60 =~ x1 + x2 + x3 + y1 + y2 + y3 + y4
  ind60 ~~ ind60
"
fit <- cfa(mod, data = PoliticalDemocracy, std.lv = TRUE, do.fit = FALSE)
```

``` r

pefa_fit <- penalized_est(
    fit,
    w = .03,
    pen_par_id = 4:10,
    se = "robust.huber.white"
)
summary(pefa_fit)
#> Penalized fit (w = 0.03, eps = 0.01, penalty = l0a): effective npar = 16.05, effective df = 11.95 (nominal df = 6).
#> lavaan 0.7-2 ended normally after 126 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        16
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

``` r

# Quick simulation to check SEs
set.seed(1234)
R <- 250
est_res <- matrix(NA, nrow = R, ncol = length(coef(pefa_fit)))
se_res <- matrix(NA, nrow = R, ncol = length(coef(pefa_fit)))

# Use the simple structure model as population
pop_model <- parTable(pefa_fit)

for (i in 1:R) {
    # Simulate data
    dat_sim <- simulateData(pop_model, sample.nobs = 200)

    # Fit penalized model
    fit_sim <- cfa(mod, data = dat_sim, std.lv = TRUE, do.fit = FALSE)
    pefa_sim <- try(
        penalized_est(
            fit_sim,
            w = .03,
            pen_par_id = 4:10,
            se = "robust.huber.white"
        ),
        silent = TRUE
    )

    if (!inherits(pefa_sim, "try-error")) {
        est_res[i, ] <- coef(pefa_sim)
        se_res[i, ] <- sqrt(diag(vcov(pefa_sim)))
    }
}

# Compare empirical SD vs mean SE for a few parameters
# (e.g., first few loadings)
res_summary <- data.frame(
    param = names(coef(pefa_fit)),
    emp_sd = apply(est_res, 2, sd, na.rm = TRUE),
    mean_se = apply(se_res, 2, mean, na.rm = TRUE)
)
```

``` r

print(res_summary, digits = 2)
#>           param emp_sd mean_se
#> 1     ind60=~x1 0.0380  0.0393
#> 2     ind60=~x2 0.0750  0.0772
#> 3     ind60=~x3 0.0787  0.0777
#> 4     ind60=~y1 0.0291  0.0060
#> 5     ind60=~y2 0.0704  0.0087
#> 6     ind60=~y3 0.0532  0.0069
#> 7     ind60=~y4 0.2631  0.1266
#> 8     dem60=~x1 0.0154  0.0152
#> 9     dem60=~x2 0.0093  0.0091
#> 10    dem60=~x3 0.0109  0.0106
#> 11    dem60=~y1 0.1586  0.1575
#> 12    dem60=~y2 0.2616  0.2445
#> 13    dem60=~y3 0.2160  0.2088
#> 14    dem60=~y4 0.2340  0.2057
#> 15       x1~~x1 0.0106  0.0116
#> 16       x2~~x2 0.0388  0.0436
#> 17       x3~~x3 0.0540  0.0536
#> 18       y1~~y1 0.3085  0.3133
#> 19       y2~~y2 0.8018  0.7797
#> 20       y3~~y3 0.5965  0.6039
#> 21       y4~~y4 0.4302  0.4333
#> 22 ind60~~dem60 0.0763  0.0686
```

``` r

# meat <- lavInspect(pefa_fit, "information.first.order")
# bread <- attr(pefa_fit, "hessian")
# vc_pefa <- solve(bread) %*% meat %*% solve(bread) / 75
# pefa_fit@vcov$vcov <- vc_pefa
# se_pefa <- sqrt(diag(vc_pefa))
# pefa_fit@ParTable$se <- 0 * pefa_fit@ParTable$est
# pefa_fit@ParTable$se[which(pefa_fit@ParTable$free > 0)] <- se_pefa
# cbind(coef(pefa_fit), sqrt(diag(vc_pefa)))
```

## Penalize non-invariance

``` r

lconfig_mod_un <- "
    # Time 1
    dem60 =~ .l1_1 * y1 + .l2_1 * y2 + .l3_1 * y3 + .l4_1 * y4
    y1 ~ .i1_1 * 1
    y2 ~ .i2_1 * 1
    y3 ~ .i3_1 * 1
    y4 ~ .i4_1 * 1
    y1 ~~ .u1_1 * y1
    y2 ~~ .u2_1 * y2
    y3 ~~ .u3_1 * y3
    y4 ~~ .u4_1 * y4
    
    # Time 2
    dem65 =~ .l1_2 * y5 + .l2_2 * y6 + .l3_2 * y7 + .l4_2 * y8
    y5 ~ .i1_2 * 1
    y6 ~ .i2_2 * 1
    y7 ~ .i3_2 * 1
    y8 ~ .i4_2 * 1
    y5 ~~ .u1_2 * y5
    y6 ~~ .u2_2 * y6
    y7 ~~ .u3_2 * y7
    y8 ~~ .u4_2 * y8
    
    # Latent variances
    dem60 ~~ 1 * dem60
    dem65 ~~ NA * dem65
    
    # Latent means
    dem60 ~ 0 * 1
    dem65 ~ NA * 1
    
    # Lag Covariances
    y1 ~~ y5
    y2 ~~ y6
    y3 ~~ y7
    y4 ~~ y8
"
# Specify the under-identified model
lconfig_fit_un <- cfa(
    lconfig_mod_un,
    data = PoliticalDemocracy,
    do.fit = FALSE,
    std.lv = TRUE,
    missing = "fiml",
    estimator = "mlr"
)
ld_id <- rbind(1:4, 13:16)
int_id <- rbind(5:8, 17:20)
pen_fit <- penalized_est(
    lconfig_fit_un,
    w = 0.03,
    pen_fn = "l0a",
    pen_diff_id = list(loadings = ld_id, intercepts = int_id),
    se = "robust.huber.white"
)
parameterEstimates(pen_fit)
#>      lhs op   rhs label    est
#> 1  dem60 =~    y1 .l1_1  2.098
#> 2  dem60 =~    y2 .l2_1  2.830
#> 3  dem60 =~    y3 .l3_1  2.561
#> 4  dem60 =~    y4 .l4_1  2.901
#> 5     y1 ~1       .i1_1  5.457
#> 6     y2 ~1       .i2_1  4.252
#> 7     y3 ~1       .i3_1  6.570
#> 8     y4 ~1       .i4_1  4.461
#> 9     y1 ~~    y1 .u1_1  2.128
#> 10    y2 ~~    y2 .u2_1  6.658
#> 11    y3 ~~    y3 .u3_1  5.385
#> 12    y4 ~~    y4 .u4_1  2.600
#> 13 dem65 =~    y5 .l1_2  2.092
#> 14 dem65 =~    y6 .l2_2  2.825
#> 15 dem65 =~    y7 .l3_2  2.572
#> 16 dem65 =~    y8 .l4_2  2.900
#> 17    y5 ~1       .i1_2  5.456
#> 18    y6 ~1       .i2_2  3.396
#> 19    y7 ~1       .i3_2  6.570
#> 20    y8 ~1       .i4_2  4.462
#> 21    y5 ~~    y5 .u1_2  2.813
#> 22    y6 ~~    y6 .u2_2  4.000
#> 23    y7 ~~    y7 .u3_2  3.615
#> 24    y8 ~~    y8 .u4_2  2.456
#> 25 dem60 ~~ dem60        1.000
#> 26 dem65 ~~ dem65        0.950
#> 27 dem60 ~1              0.000
#> 28 dem65 ~1             -0.147
#> 29    y1 ~~    y5        0.838
#> 30    y2 ~~    y6        1.825
#> 31    y3 ~~    y7        1.205
#> 32    y4 ~~    y8        0.286
#> 33 dem60 ~~ dem65        0.917
```

Compared to scalar invariance model

``` r

lscalar_mod <- "
    # Time 1
    dem60 =~ .l1 * y1 + .l2 * y2 + .l3 * y3 + .l4 * y4
    y1 ~ .i1 * 1
    y2 ~ .i2 * 1
    y3 ~ .i3 * 1
    y4 ~ .i4 * 1
    y1 ~~ .u1_1 * y1
    y2 ~~ .u2_1 * y2
    y3 ~~ .u3_1 * y3
    y4 ~~ .u4_1 * y4
    
    # Time 2
    dem65 =~ .l1 * y5 + .l2 * y6 + .l3 * y7 + .l4 * y8
    y5 ~ .i1 * 1
    y6 ~ .i2 * 1
    y7 ~ .i3 * 1
    y8 ~ .i4 * 1
    y5 ~~ .u1_2 * y5
    y6 ~~ .u2_2 * y6
    y7 ~~ .u3_2 * y7
    y8 ~~ .u4_2 * y8
    
    # Latent variances
    dem60 ~~ 1 * dem60
    dem65 ~~ NA * dem65
    
    # Latent means
    dem60 ~ 0 * 1
    dem65 ~ NA * 1
    
    # Lag Covariances
    y1 ~~ y5
    y2 ~~ y6
    y3 ~~ y7
    y4 ~~ y8
"
lscalar_fit <- cfa(
    lscalar_mod,
    data = PoliticalDemocracy,
    std.lv = TRUE,
    missing = "fiml",
    estimator = "mlr"
)
parameterEstimates(lscalar_fit)
#>      lhs op   rhs label    est    se      z pvalue ci.lower ci.upper
#> 1  dem60 =~    y1   .l1  2.085 0.211  9.884  0.000    1.672    2.498
#> 2  dem60 =~    y2   .l2  2.896 0.284 10.189  0.000    2.339    3.453
#> 3  dem60 =~    y3   .l3  2.552 0.246 10.366  0.000    2.069    3.034
#> 4  dem60 =~    y4   .l4  2.889 0.225 12.842  0.000    2.448    3.330
#> 5     y1 ~1         .i1  5.508 0.289 19.037  0.000    4.941    6.075
#> 6     y2 ~1         .i2  3.796 0.429  8.856  0.000    2.956    4.636
#> 7     y3 ~1         .i3  6.670 0.353 18.915  0.000    5.979    7.361
#> 8     y4 ~1         .i4  4.554 0.366 12.458  0.000    3.838    5.270
#> 9     y1 ~~    y1 .u1_1  2.141 0.488  4.384  0.000    1.184    3.098
#> 10    y2 ~~    y2 .u2_1  6.840 1.444  4.738  0.000    4.010    9.669
#> 11    y3 ~~    y3 .u3_1  5.424 1.101  4.927  0.000    3.266    7.582
#> 12    y4 ~~    y4 .u4_1  2.623 0.641  4.091  0.000    1.366    3.879
#> 13 dem65 =~    y5   .l1  2.085 0.211  9.884  0.000    1.672    2.498
#> 14 dem65 =~    y6   .l2  2.896 0.284 10.189  0.000    2.339    3.453
#> 15 dem65 =~    y7   .l3  2.552 0.246 10.366  0.000    2.069    3.034
#> 16 dem65 =~    y8   .l4  2.889 0.225 12.842  0.000    2.448    3.330
#> 17    y5 ~1         .i1  5.508 0.289 19.037  0.000    4.941    6.075
#> 18    y6 ~1         .i2  3.796 0.429  8.856  0.000    2.956    4.636
#> 19    y7 ~1         .i3  6.670 0.353 18.915  0.000    5.979    7.361
#> 20    y8 ~1         .i4  4.554 0.366 12.458  0.000    3.838    5.270
#> 21    y5 ~~    y5 .u1_2  2.824 0.577  4.895  0.000    1.693    3.954
#> 22    y6 ~~    y6 .u2_2  4.032 0.818  4.931  0.000    2.429    5.634
#> 23    y7 ~~    y7 .u3_2  3.650 0.639  5.713  0.000    2.398    4.903
#> 24    y8 ~~    y8 .u4_2  2.482 0.705  3.522  0.000    1.101    3.863
#> 25 dem60 ~~ dem60        1.000 0.000     NA     NA    1.000    1.000
#> 26 dem65 ~~ dem65        0.947 0.097  9.745  0.000    0.756    1.137
#> 27 dem60 ~1              0.000 0.000     NA     NA    0.000    0.000
#> 28 dem65 ~1             -0.210 0.067 -3.115  0.002   -0.342   -0.078
#> 29    y1 ~~    y5        0.845 0.433  1.953  0.051   -0.003    1.694
#> 30    y2 ~~    y6        1.670 0.934  1.788  0.074   -0.161    3.502
#> 31    y3 ~~    y7        1.206 0.646  1.868  0.062   -0.060    2.472
#> 32    y4 ~~    y8        0.261 0.465  0.563  0.574   -0.649    1.172
#> 33 dem60 ~~ dem65        0.918 0.057 16.191  0.000    0.807    1.030
```
