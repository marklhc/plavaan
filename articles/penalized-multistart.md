# Multistart Penalized Estimation

``` r

library(plavaan)
library(lavaan)
#> This is lavaan 0.7-2
#> lavaan is FREE software! Please report any bugs.
data(PoliticalDemocracy)
```

## Why multistart?

Penalized objectives such as `l0a` and `alf` are non-convex near zero.
When combined with a generic optimizer like
[`nlminb()`](https://rdrr.io/r/stats/nlminb.html), the risk of settling
in a local optimum is real. Running the optimization from multiple
starting values helps avoid this problem: if different starts converge
to the same solution, you can be more confident it is (locally) optimal;
if they diverge, the multistart wrapper selects the best among them.

## Single-start baseline

Start with the longitudinal invariance example from \[penalized_est()\].
We fit a two-wave CFA and penalize differences in loadings and
intercepts across time:

``` r

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

fit_un <- cfa(model, data = PoliticalDemocracy, std.lv = TRUE,
              meanstructure = TRUE, do.fit = FALSE)

pt <- parTable(fit_un)
load_60 <- pt$free[pt$op == "=~" & pt$lhs == "dem60"]
load_65 <- pt$free[pt$op == "=~" & pt$lhs == "dem65"]
int_60  <- pt$free[pt$op == "~1" & pt$lhs %in% c("y1", "y2", "y3", "y4")]
int_65  <- pt$free[pt$op == "~1" & pt$lhs %in% c("y5", "y6", "y7", "y8")]

fit_single <- penalized_est(
    x = fit_un,
    w = 0.03,
    pen_diff_id = list(
        loadings   = rbind(load_60, load_65),
        intercepts = rbind(int_60, int_65)
    ),
    pen_fn = "l0a"
)

fit_single@optim$fx
#> [1] 0.2100911
```

## Multistart with random jittering

[`penalized_est_multistart()`](https://marklhc.github.io/plavaan/reference/penalized_est_multistart.md)
wraps
[`penalized_est()`](https://marklhc.github.io/plavaan/reference/penalized_est.md)
and tries several starting vectors. The first start is always lavaan’s
default (unperturbed), so multistart can never do worse than the
single-start baseline:

``` r

set.seed(3821)

fit_multi <- penalized_est_multistart(
    x = fit_un,
    w = 0.03,
    pen_diff_id = list(
        loadings   = rbind(load_60, load_65),
        intercepts = rbind(int_60, int_65)
    ),
    pen_fn = "l0a",
    n_starts = 10,
    verbose = TRUE,
    eps = .0001
)
#> Start 1 / 10...
#>   Start 1: converged, objective = 0.279207
#> Start 2 / 10...
#>   Start 2: converged, objective = 0.279207
#> Start 3 / 10...
#>   Start 3: converged, objective = 0.279207
#> Start 4 / 10...
#>   Start 4: converged, objective = 0.279207
#> Start 5 / 10...
#>   Start 5: converged, objective = 0.279207
#> Start 6 / 10...
#>   Start 6: converged, objective = 0.279207
#> Start 7 / 10...
#>   Start 7: converged, objective = 0.279207
#> Start 8 / 10...
#>   Start 8: converged, objective = 0.279207
#> Start 9 / 10...
#>   Start 9: converged, objective = 0.279207
#> Start 10 / 10...
#>   Start 10: converged, objective = 0.279207

# Inspect the spread of objective values across starts
attr(fit_multi, "multistart")
#>    start_id objective converged
#> 10       10 0.2792071      TRUE
#> 3         3 0.2792071      TRUE
#> 8         8 0.2792071      TRUE
#> 1         1 0.2792071      TRUE
#> 5         5 0.2792071      TRUE
#> 6         6 0.2792071      TRUE
#> 2         2 0.2792071      TRUE
#> 4         4 0.2792071      TRUE
#> 7         7 0.2792071      TRUE
#> 9         9 0.2792071      TRUE
```

The `multistart` attribute is a data frame with one row per start,
showing the final penalized objective (`objective`) and whether it
converged. Rows are sorted by ascending objective.

## Comparing local solutions

When starts lead to different parameter estimates, it can be useful to
inspect all solutions:

``` r

set.seed(3821)

fit_all <- penalized_est_multistart(
    x = fit_un,
    w = 0.03,
    pen_diff_id = list(
        loadings   = rbind(load_60, load_65),
        intercepts = rbind(int_60, int_65)
    ),
    pen_fn = "l0a",
    n_starts = 5,
    keep_all = TRUE
)

# Extract loadings from the top-3 solutions for comparison
all_fits <- attr(fit_all, "all_fits")
ms_table <- attr(fit_all, "multistart")[order(attr(fit_all, "multistart")$objective), ]

loadings_top3 <- lapply(head(all_fits[ms_table$start_id[1:3]], 3), function(f) {
    coefs <- coef(f)[grep("^dem60=~", names(coef(f)))]
    round(coefs, 4)
})
names(loadings_top3) <- paste0("Start ", ms_table$start_id[1:3])
loadings_top3
#> $`Start 1`
#> dem60=~y1 dem60=~y2 dem60=~y3 dem60=~y4 
#>    2.0980    2.8300    2.5609    2.9006 
#> 
#> $`Start 2`
#> dem60=~y1 dem60=~y2 dem60=~y3 dem60=~y4 
#>    2.0980    2.8300    2.5609    2.9006 
#> 
#> $`Start 5`
#> dem60=~y1 dem60=~y2 dem60=~y3 dem60=~y4 
#>    2.0981    2.8300    2.5609    2.9006
```

## Using custom starting values

Instead of random jittering, you can supply your own starting vectors
via `starts`:

``` r

base <- lavaan::lav_export_estimation(fit_un)$starting_values
my_starts <- rbind(
    base,                              # unperturbed
    base + runif(length(base), -0.1, 0.1)  # perturbed
)

fit_custom <- penalized_est_multistart(
    x = fit_un,
    w = 0.03,
    pen_diff_id = list(
        loadings   = rbind(load_60, load_65),
        intercepts = rbind(int_60, int_65)
    ),
    starts = my_starts
)
#> Custom starting values provided (2 starts). Ignoring n_starts.

attr(fit_custom, "multistart")
#>   start_id objective converged
#> 1        1 0.2100911      TRUE
#> 2        2 0.2100911      TRUE
```

## Reproducibility

The function does not call
[`set.seed()`](https://rdrr.io/r/base/Random.html) internally. Set your
own seed before calling for reproducibility:

``` r

set.seed(7293)
fit_repro <- penalized_est_multistart(
    x = fit_un, w = 0.03,
    pen_diff_id = list(loadings = rbind(load_60, load_65)),
    n_starts = 5
)
```

## Parallel execution

Multistart runs sequentially for simplicity. For parallel execution,
call `penalized_est(start = ...)` directly with your preferred backend
(`future.apply`, `parallel`, etc.):

``` r

library(future.apply)
plan(multicore)

# Note: random_start() is an internal helper without API stability guarantees.
starts <- plavaan:::random_start(fit_un, n = 10)[-1, , drop = FALSE]  # skip base
fits <- future_lapply(seq_len(nrow(starts)), function(i) {
    penalized_est(
        x = fit_un, w = 0.03,
        pen_diff_id = list(loadings = rbind(load_60, load_65)),
        start = starts[i, ]
    )
})
```
