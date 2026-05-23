# betaARMA 1.2.0 (development)

## New features

* Add S3 method `plot.barma()` to generate diagnostic plots for fitted `barma`
  objects. It displays a four-panel grid by default: observed vs. fitted values,
  residuals over time, residual ACF, and a residual histogram with a density 
  overlay.
* Add `penalty` argument to `barma()`, enabling ridge-penalized conditional
  maximum likelihood estimation as proposed by Cribari-Neto, Costa and
  Fonseca (2025, BJPS). When `penalty = TRUE`, the penalized log-likelihood
  is optimized and the penalized FIM is returned.
* Add internal helper `utils_internal.R` with `.barma_ridge_lambda()` to
  centralize the ridge penalty scalar `lambda_n = 1/(n-a)^0.9`, shared
  across `loglik_barma()`, `score_vector_barma()`, and `fim_barma()`.
* Add `utils_optimizer.R` with internal helpers `make_par_indices()` and
  `run_optimizer()`, centralizing parameter indexing and optimization logic
  previously embedded in `barma()`.

## Improvements

* Forward `penalty` argument through `run_optimizer()` to `loglik_barma()`
  and `score_vector_barma()`, ensuring the penalized likelihood is used
  consistently during optimization when requested.
* Forward `penalty` argument to `fim_barma()` in `barma()`, so the returned
  Fisher Information Matrix reflects penalization when `penalty = TRUE`.

## Documentation

* Add two new vignettes: one detailing the numerical stability and ridge 
  penalization approach (Cribari-Neto, Costa, and Fonseca, 2025), and 
  another demonstrating applied hydro-environmental time series modeling 
  (Costa, Cribari-Neto, and Scher, 2024).
* Add Cribari-Neto, Costa and Fonseca (2025, BJPS, DOI: 10.1214/25-BJPS645)
  to `@references` in `barma()`, reflecting the new `penalty` argument.
* Update `@param optimization` in `barma()` to use `\describe` instead of
  `\itemize` for correct rendering of named list fields in HTML help pages.

# betaARMA 1.1.0 (2026-04-14)

## Bug fixes

* Fixed incorrect computation of `ystar` in `score_vector_barma()`; previously the estimates were only correct for the logit link function.

## New features

* Export `loglik_barma()`, `score_vector_barma()`, and `fim_barma()` to allow advanced user access to model internals.

## Breaking changes

* Deprecate `ar = NA` and `ma = NA` for absent AR/MA components. Use `ar = integer(0)` and `ma = integer(0)` or omit the arguments. Backward compatibility is preserved in `barma()` with a deprecation warning.

## Improvements

* Simplify internal handling of absent AR/MA components using vectorized
  operations (`drop(crossprod())`) in `loglik_barma()`, `score_vector_barma()`,
  and `fim_barma()`, replacing conditional `if/else` branching.
* Standardize parameter order across core functions to
  `(alpha, varphi, theta, beta, phi)`.
* Set default link to `"logit"` in `make_link_structure()` and all exported
  functions; improve error messages for unsupported links.
* Rename internal helper `._get_phi_start` to `.get_phi_start` for consistency
  with R conventions.
* Improve documentation of precision parameter initialization in
  `start_values()`.

# betaARMA 1.0.1 (2026-03-29)

## Documentation

* Update README for CRAN release.
* Add CRAN status badge.
* Use `install.packages()` as primary installation method.
* Update citation to version 1.0.1 with proper BibTeX key.
* Exclude generated `README.html` from version control.

# betaARMA 1.0.0 (2026-03-21)

* Initial CRAN submission.

