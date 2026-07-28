# plavaan 0.0.2

* Added `penalized_est_multistart()` and support for custom optimizer starting values via `start` in `penalized_est()`.
* Minor changes to speed up the optimization process.

## Breaking Changes

* **Penalty metric for factor loadings:** The penalty for cross-group differences in factor loadings is now computed on the **original scale** ($\lambda_{g1} - \lambda_{g2}$) rather than the log scale ($\log(\lambda_{g1}) - \log(\lambda_{g2})$).

# plavaan 0.0.1

* Initial CRAN submission.
