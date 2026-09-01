# plavaan (development version)

* `penalized_est()` gains `eps` control for built-in penalties, including `eps = "telescoping"` over a decreasing epsilon sequence. Telescoping reuses the original starting values by default; set `telescoping_control = list(warm_start = TRUE)` to warm-start later stages from prior estimates.
* New `effective_df()` for penalized fits: reports the effective number of parameters per penalty component and the effective model degrees of freedom, using the soft-count property of the penalties.
* `effective_df()` (and the fit-evaluation degrees of freedom) now count sample statistics correctly for models with ordinal/categorical indicators (e.g. WLSMV fits with thresholds). The count now uses lavaan's own `lav_pt_ndat()` where available (with a structural fallback), so it accounts for ordinal thresholds, composites, and the correlation metric, and always matches lavaan's own model df. The previous count used a continuous-variables formula that miscounted such models. Results for purely continuous models are unchanged.
* Fits returned by `penalized_est()` are now `plavaan` objects. `summary()` reports the effective number of parameters and degrees of freedom, and — when fit evaluation is enabled with `test = "Chisq"` (or `"SatorraBentler"`) — `fitmeasures()` and the chi-square test in `summary()` additionally report fit indices at the effective degrees of freedom via a frozen refit at the penalized estimates. This fit evaluation is experimental and disabled by default (`test = "none"`, the new default for the `test` argument); an experimental notice is shown whenever it is used.
* `summary()` on a penalized fit prints standard errors (and the resulting z-values and p-values) when `penalized_est()` was called with `se = "robust.huber.white"`. A development-version change that computed the `summary()` from a frozen refit had suppressed the standard-error column even when standard errors were requested; this is now corrected. `parameterEstimates()` was unaffected.
* `fitmeasures()` on a penalized fit now honours its arguments for selecting and formatting fit indices (e.g. `fitmeasures(fit, "rmsea")`, `output = "matrix"`, `fm_args`). A development-version change had declared the method as `(object, ...)`, so lavaan's `fitmeasures()` generic's named arguments were silently dropped and every call returned the full set of fit measures; the method now declares and forwards them.

# plavaan 0.0.2

* Added `penalized_est_multistart()` and support for custom optimizer starting values via `start` in `penalized_est()`.
* Minor changes to speed up the optimization process.

## Breaking Changes

* **Penalty metric for factor loadings:** The penalty for cross-group differences in factor loadings is now computed on the **original scale** ($\lambda_{g1} - \lambda_{g2}$) rather than the log scale ($\log(\lambda_{g1}) - \log(\lambda_{g2})$).

# plavaan 0.0.1

* Initial CRAN submission.
