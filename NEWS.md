# betaARMA 1.1.0 (2026-04-14)

## Bug fixes

* Fixed incorrect computation of `ystar` in `score_vector_barma()`; previously the estimates were only correct for the logit link function.

## New features

* Export `loglik_barma()`, `score_vector_barma()`, and `fim_barma()` to allow advanced user access to model internals.

## Breaking changes

* Deprecate `ar = NA` and `ma = NA` for absent AR/MA components. Use `ar = integer(0)` and `ma = integer(0)` or omit the arguments. Backward compatibility is preserved in `barma()` with a deprecation warning.

## Improvements

* Simplify internal handling of absent AR/MA components using vectorized operations.
* Standardize parameter order across core functions to `(alpha, varphi, theta, beta, phi)`.
* Set default link to `"logit"` in `make_link_structure()` and improve error messages for unsupported links.
* Rename internal helper `._get_phi_start` to `.get_phi_start` for consistency with R conventions.
* Improve documentation of precision parameter initialization in `start_values()`.

# betaARMA 1.0.1 (2026-03-29)

## Documentation

* Update README for CRAN release.
* Add CRAN status badge.
* Use `install.packages()` as primary installation method.
* Update citation to version 1.0.1 with proper BibTeX key.
* Exclude generated `README.html` from version control.

# betaARMA 1.0.0 (2026-03-21)

* Initial CRAN submission.

