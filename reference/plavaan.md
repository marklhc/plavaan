# Penalized lavaan Fit Objects (‘plavaan’ Class)

An S4 subclass of the ‘lavaan’ class, defined by
`setClass("plavaan", contains = "lavaan")`. Objects returned by
[`penalized_est()`](https://marklhc.github.io/plavaan/reference/penalized_est.md)
are ‘plavaan’ objects: they carry a ‘penalized’ attribute recording the
penalty specification, a shared ‘plavaan.cache’ environment attribute
(used to cache the frozen refit), a
[`fitmeasures()`](https://rdrr.io/pkg/lavaan/man/fitMeasures.html)
method that reports fit indices at the effective degrees of freedom (see
[`effective_df()`](https://marklhc.github.io/plavaan/reference/effective_df.md)),
and a [`summary()`](https://rdrr.io/r/base/summary.html) method that
prints the frozen fit's summary together with the effective parameter
and degrees-of-freedom counts. Fit evaluation
([`fitmeasures()`](https://rdrr.io/pkg/lavaan/man/fitMeasures.html) and
the chi-square test in
[`summary()`](https://rdrr.io/r/base/summary.html)) is experimental and
only available when the fit was created with a non-`"none"` `test` (see
[`penalized_est()`](https://marklhc.github.io/plavaan/reference/penalized_est.md)).
