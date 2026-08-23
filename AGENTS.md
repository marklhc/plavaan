# AGENTS.md

plavaan is a CRAN R package that adds penalized estimation to `lavaan`
(Robitzsch 2023 style): it re-optimizes a lavaan model with `nlminb` on
a penalized objective.

## Commands

- Run all tests: `devtools::test()` (testthat edition 3). Full suite
  takes ~3s.
- Run one suite: `devtools::test(filter = "multistart")`.
- Do NOT verify test-multistart.R with
  [`testthat::test_file()`](https://testthat.r-lib.org/reference/test_file.html)
  directly: it has a top-level `skip_on_cran()`, and `test_file()` does
  not set `NOT_CRAN`, so the whole suite silently skips (shows
  `SKIP 1`). Use `devtools::test()` or `Sys.setenv(NOT_CRAN = "true")`.
- Regenerate docs: `devtools::document()` (roxygen2 8.1.0 pinned in
  DESCRIPTION) — updates `man/` and `NAMESPACE`. Never hand-edit those.
- `README.md` is generated from `README.Rmd` — edit the `.Rmd`.
- Before a CRAN submission: update `NEWS.md`, run `R CMD check`, update
  `cran-comments.md`.
- Dependencies: `lavaan (>= 0.6-15)`, `numDeriv` (Imports);
  `testthat (>= 3.1.7)` (Suggests — tests use
  `with_mocked_bindings`/`capture_warnings`). Tests use the built-in
  [`lavaan::PoliticalDemocracy`](https://rdrr.io/pkg/lavaan/man/PoliticalDemocracy.html)
  data; no external fixtures or services.

## Architecture

- Three source files: `R/penalized.R` (core
  [`penalized_est()`](https://marklhc.github.io/plavaan/reference/penalized_est.md),
  objective/gradient assembly, telescoping), `R/multistart.R`
  (`random_start()`,
  [`penalized_est_multistart()`](https://marklhc.github.io/plavaan/reference/penalized_est_multistart.md)),
  `R/alf.R` (penalties `l0a`/`alf` + gradients,
  [`composite_pair_loss()`](https://marklhc.github.io/plavaan/reference/composite_pair_loss.md)).
- Canonical usage pattern: fit a dry-run lavaan model with
  `do.fit = FALSE` → `parTable()` → take the **`free` column** values as
  `pen_par_id` / `pen_diff_id` entries →
  [`penalized_est()`](https://marklhc.github.io/plavaan/reference/penalized_est.md).
  These IDs are positions in the free-parameter vector (matching
  `lavaan::coef()` order), not parTable row numbers — a common
  off-by-design mistake.
- `pen_diff_id` is a named list of matrices: one row per group/time
  point, one column per matched parameter. The penalty is the sum over
  pairwise row differences, rescaled by (nrow - 1)/ncombn. Structural
  `NA`s in the matrices mean “parameter absent in that row” and must be
  preserved.
- The returned lavaan object is rebuilt with `do.fit = FALSE` (see
  `make_penalized_fit()`, R/penalized.R:352):
  [`summary()`](https://rdrr.io/r/base/summary.html) works but SEs are
  absent unless `se = "robust.huber.white"` (experimental sandwich
  approximation — interpret with caution).
- `eps = "telescoping"` runs a log-spaced sequence of epsilons, reusing
  the original start by default (`warm_start = TRUE` chains starts);
  stops when max \|Δparam\| ≤ 5e-4.
- [`penalized_est_multistart()`](https://marklhc.github.io/plavaan/reference/penalized_est_multistart.md)
  is sequential by design; for parallelism call
  `penalized_est(start = ...)` directly per start.

## Conventions

- Difference penalties on factor loadings are computed on the **original
  scale**, not log scale (breaking change in 0.0.2 — see `NEWS.md`).
- Vignettes in `vignettes/` are the primary usage docs; keep them
  consistent with API changes.
- pkgdown site deploys automatically via
  `.github/workflows/pkgdown.yaml`; site config is `_pkgdown.yml`.
- Gradient tests validate analytic gradients against
  [`numDeriv::grad()`](https://rdrr.io/pkg/numDeriv/man/grad.html) — new
  penalty/gradient code should follow this pattern.
