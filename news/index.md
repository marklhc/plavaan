# Changelog

## plavaan (development version)

- [`penalized_est()`](https://marklhc.github.io/plavaan/reference/penalized_est.md)
  gains `eps` control for built-in penalties, including
  `eps = "telescoping"` over a decreasing epsilon sequence. Telescoping
  reuses the original starting values by default; set
  `telescoping_control = list(warm_start = TRUE)` to warm-start later
  stages from prior estimates.
- New
  [`effective_df()`](https://marklhc.github.io/plavaan/reference/effective_df.md)
  for penalized fits: reports the effective number of parameters per
  penalty component and the effective model degrees of freedom, using
  the soft-count property of the penalties.
- [`effective_df()`](https://marklhc.github.io/plavaan/reference/effective_df.md)
  (and the fit-evaluation degrees of freedom) now count sample
  statistics correctly for models with ordinal/categorical indicators
  (e.g. WLSMV fits with thresholds). The count now uses lavaan’s own
  [`lav_pt_ndat()`](https://rdrr.io/pkg/lavaan/man/lav_partable.html)
  where available (with a structural fallback), so it accounts for
  ordinal thresholds, composites, and the correlation metric, and always
  matches lavaan’s own model df. The previous count used a
  continuous-variables formula that miscounted such models. Results for
  purely continuous models are unchanged.
- Fits returned by
  [`penalized_est()`](https://marklhc.github.io/plavaan/reference/penalized_est.md)
  are now `plavaan` objects.
  [`summary()`](https://rdrr.io/r/base/summary.html) reports the
  effective number of parameters and degrees of freedom, and — when fit
  evaluation is enabled with `test = "Chisq"` (or `"SatorraBentler"`) —
  [`fitmeasures()`](https://rdrr.io/pkg/lavaan/man/fitMeasures.html) and
  the chi-square test in
  [`summary()`](https://rdrr.io/r/base/summary.html) additionally report
  fit indices at the effective degrees of freedom via a frozen refit at
  the penalized estimates. This fit evaluation is experimental and
  disabled by default (`test = "none"`, the new default for the `test`
  argument); an experimental notice is shown whenever it is used.

## plavaan 0.0.2

CRAN release: 2026-07-28

- Added
  [`penalized_est_multistart()`](https://marklhc.github.io/plavaan/reference/penalized_est_multistart.md)
  and support for custom optimizer starting values via `start` in
  [`penalized_est()`](https://marklhc.github.io/plavaan/reference/penalized_est.md).
- Minor changes to speed up the optimization process.

### Breaking Changes

- **Penalty metric for factor loadings:** The penalty for cross-group
  differences in factor loadings is now computed on the **original
  scale** ($`\lambda_{g1} - \lambda_{g2}`$) rather than the log scale
  ($`\log(\lambda_{g1}) - \log(\lambda_{g2})`$).

## plavaan 0.0.1

CRAN release: 2025-12-30

- Initial CRAN submission.
