# Composite Pairwise Loss Function

Computes the total loss across all pairwise combinations of rows in a
matrix.

## Usage

``` r
composite_pair_loss(x, fun, trans = identity, rescale = "df", ...)
```

## Arguments

- x:

  A numeric vector, matrix, or data frame. If not a matrix, it will be
  coerced to one after applying the transformation function.

- fun:

  A function to compute the loss for each pairwise difference. The
  package supports the alignment loss (`alf`) and the approximate L0
  penalty (`l0a`), but users can provide custom functions as well.

- trans:

  A transformation function to apply to `x` before computing pairwise
  differences. Default is `identity` (no transformation).

- rescale:

  Either `"df"` (default) to rescale the total loss by
  `(nrow - 1) / ncombn(nrow, 2)`, where `nrow` is the number of rows, or
  a numeric value (likely between 0 and 1) to multiply the total loss
  by. With `rescale = "df"` and a 0/1-valued `fun` such as `l0a`, the
  returned value approximates the number of degrees of freedom of
  non-invariance in the block: `0` when all rows are equal (full
  invariance) and `nrow - 1` per column when all rows differ (fully
  free). It soft-counts how many of the `nrow - 1` pairwise-invariance
  constraints per column are violated.

- ...:

  Additional arguments passed to the loss function `fun`.

## Value

A numeric scalar representing the sum of losses across all pairwise
combinations of rows.

## Details

The function works by:

1.  Applying the transformation function `trans` to the input `x`

2.  Converting the result to a matrix

3.  Generating all possible pairwise combinations of row indices

4.  Computing the difference between each pair of rows

5.  Applying the loss function `fun` to each difference

6.  Summing all the individual losses

[`effective_df()`](https://marklhc.github.io/plavaan/reference/effective_df.md)
builds on this degrees-of-freedom interpretation (`rescale = "df"` with
a 0/1-valued `fun` such as `l0a`) to report the effective number of
parameters, and hence the effective model degrees of freedom, of
penalized fits returned by
[`penalized_est()`](https://marklhc.github.io/plavaan/reference/penalized_est.md).

## Examples

``` r
# Example with a simple matrix
x <- matrix(runif(12), nrow = 4)
composite_pair_loss(x, fun = alf)
#> [1] 5.075207

# Example with log transformation and L2 loss
composite_pair_loss(x, fun = function(x) x^2, trans = log)
#> [1] 35.18578
```
