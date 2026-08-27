# Print Method for Effective Degrees of Freedom Tables

Prints the component table from
[`effective_df()`](https://marklhc.github.io/plavaan/reference/effective_df.md)
followed by summary lines with the number of sample statistics, the
nominal and effective model degrees of freedom, and the penalty
settings.

## Usage

``` r
# S3 method for class 'plavaan_efdf'
print(x, ...)
```

## Arguments

- x:

  An object of class `plavaan_efdf`, typically the return value of
  [`effective_df()`](https://marklhc.github.io/plavaan/reference/effective_df.md).

- ...:

  Passed to
  [`print.data.frame()`](https://rdrr.io/r/base/print.dataframe.html).
