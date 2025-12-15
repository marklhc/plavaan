# Loss functions

For small eps this provides a smooth, numerically stable approximation
of \|x\|^(1/2) (i.e. the square root of the absolute value). The
function is vectorized over x.

## Usage

``` r
alf(x, eps = 0.001)

l0a(x, eps = 0.01)
```

## Arguments

- x:

  Numeric vector. Input values to transform.

- eps:

  Positive numeric scalar (default .001 for `alf()` and .01 for
  `l0a()`). Small regularization constant to avoid non-differentiability
  and division-by-zero issues.

## Value

Numeric vector of the same length as x.

## Details

The ALF, (x^2 + eps)^(1/4), is useful when a smooth surrogate for
sqrt(\|x\|) is required (for optimization or regularization) while
maintaining numerical stability near x = 0.

L0a, x^2/(x^2 + eps), is an approximation of the L0 penalty.

## Examples

``` r
alf(0)
#> [1] 0.1778279
alf(c(-4, -1, 0, 1, 4))
#> [1] 2.0000312 1.0002499 0.1778279 1.0002499 2.0000312
alf(0.5, eps = 1e-6)
#> [1] 0.7071075
l0a(0)
#> [1] 0
l0a(c(0, 1e-3, 0.1, 1))
#> [1] 0.00000000 0.00009999 0.50000000 0.99009901
l0a(c(-2, 0, 2), eps = 1e-4)
#> [1] 0.999975 0.000000 0.999975
```
