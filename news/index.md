# Changelog

## plavaan (development version)

- [`penalized_est()`](https://marklhc.github.io/plavaan/reference/penalized_est.md)
  gains `eps` control for built-in penalties, including
  `eps = "telescoping"` over a decreasing epsilon sequence. Telescoping
  reuses the original starting values by default; set
  `telescoping_control = list(warm_start = TRUE)` to warm-start later
  stages from prior estimates.

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
