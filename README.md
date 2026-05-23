
# betaARMA: Beta Autoregressive Moving Average Models

[![CRAN
status](https://www.r-pkg.org/badges/version/betaARMA)](https://CRAN.R-project.org/package=betaARMA)
[![DOI](https://img.shields.io/badge/DOI-10.32614%2FCRAN.package.betaARMA-blue)](https://doi.org/10.32614/CRAN.package.betaARMA)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/betaARMA)](https://CRAN.R-project.org/package=betaARMA)
[![CRAN all-time
downloads](https://cranlogs.r-pkg.org/badges/grand-total/betaARMA)](https://CRAN.R-project.org/package=betaARMA)
[![R-CMD-check](https://github.com/Everton-da-Costa/betaARMA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Everton-da-Costa/betaARMA/actions/workflows/R-CMD-check.yaml)
[![R-hub](https://github.com/Everton-da-Costa/betaARMA/actions/workflows/rhub.yaml/badge.svg)](https://github.com/Everton-da-Costa/betaARMA/actions/workflows/rhub.yaml)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

The `betaARMA` package provides a comprehensive suite of tools for
fitting, diagnosing, and forecasting Beta Autoregressive Moving Average
models. These models are specifically tailored for time series data that
assume values in the standard unit interval (0, 1), such as rates,
proportions, and relative indices.

## Features

- **Estimation:** Supports both standard Conditional Maximum Likelihood
  Estimation (CMLE) and Penalized Conditional Maximum Likelihood
  Estimation (PCMLE) via ridge penalization to handle numerical
  instabilities and flat likelihood surfaces (Cribari-Neto, Costa, and
  Fonseca, 2025).
- **Diagnostics:** Highly optimized, fully vectorized implementations of
  the Ljung-Box and Monti portmanteau tests, correctly adjusted for the
  effective sample size and the number of estimated parameters.
- **Forecasting:** Built-in methods for out-of-sample predictions. True
  to the theoretical properties of the betaARMA model, all forecasts
  naturally remain strictly within the (0, 1) bounds.

------------------------------------------------------------------------

## What’s New in v1.2.0

- **Ridge Penalization:** New `penalty` argument in `barma()` enables
  Penalized CMLE (PCMLE) as proposed by Cribari-Neto, Costa, and Fonseca
  (2025). Addresses numerical instability caused by flat log-likelihood
  surfaces in challenging data structures.
- **Advanced Optimization Backends:** Users can now specify the
  optimization algorithm via the `optimization` argument in `barma()`.
  Supported methods include `stats::optim` (BFGS and L-BFGS-B with
  bounds) and `lbfgs::lbfgs`.
- **Enhanced Model Summaries:** The `summary()` method has been
  significantly upgraded. It now outputs significance levels (stars) for
  coefficients, explicitly details the optimization procedure and
  convergence status used, and reports the $\lambda$ penalty parameter
  when ridge penalization is applied.
- **Diagnostic Plots:** New `plot()` method for `barma` objects displays
  a comprehensive diagnostic grid.
- **New Vignette:** A new, comprehensive applied vignette is now
  included. It demonstrates the handling of numerical instability, the
  application of ridge penalization (PCMLE), and a full modeling
  workflow using monthly relative humidity data from Brasília, Brazil.

------------------------------------------------------------------------

## Table of Contents

- [Project Motivation](#-project-motivation)
- [Foundational Literature](#-foundational-literature)
- [Key Features](#-key-features)
- [Repository Structure](#-repository-structure)
- [Installation](#️-installation)
- [Getting Started](#-getting-started)
- [Citation](#-citation)
- [Contributing](#-contributing)
- [License](#-license)
- [Contact](#-contact)

------------------------------------------------------------------------

## Project Motivation

Modeling time series data bounded within the unit interval $(0, 1)$
presents unique challenges. Standard Gaussian methods (like ARIMA) are
often inappropriate because they do not respect the natural boundaries
of the data, potentially leading to fitted values or forecasts outside
the admissible range.

The **$\beta\text{ARMA}$ model** addresses this by assuming the
conditional distribution of the variable follows a Beta law. While the
theoretical foundations exist, there has been a need for a modern,
robust R package that:

1.  **Unifies the workflow** for $\beta$AR, $\beta$MA, and $\beta$ARMA
    specifications.
2.  **Ensures stability** through numerical optimization with analytic
    gradients, including ridge penalization for challenging data
    structures.
3.  **Provides a standard API** consistent with popular time series
    tools (like `forecast`).

This project aims to fill that gap, serving as a go-to resource for
hydrologists, economists, and data scientists working with bounded data.

------------------------------------------------------------------------

## Foundational Literature

The foundational methodology implemented in this package was published
in **TEST**, journal in the fields of Statistics and Probability.

[![SCImago Journal & Country
Rank](https://www.scimagojr.com/journal_img.php?id=14882)](https://www.scimagojr.com/journalsearch.php?q=14882&tip=sid)

**Key Publications:**

- **Rocha, A. V., & Cribari-Neto, F. (2009).** “Beta autoregressive
  moving average models.” *TEST*, 18(3), 529–545.
  [doi:10.1007/s11749-008-0112-z](https://doi.org/10.1007/s11749-008-0112-z)
- **Rocha, A. V., & Cribari-Neto, F. (2017).** “Erratum to: Beta
  autoregressive moving average models.” *TEST*, 26(2), 451–459.
  [doi:10.1007/s11749-017-0528-4](https://doi.org/10.1007/s11749-017-0528-4)
- **Cribari-Neto, F., Costa, E., & Fonseca, R. V. (2025).** “Numerical
  stability enhancements in beta autoregressive moving average model
  estimation.” *Brazilian Journal of Probability and Statistics*, 39(4),
  410–437. [doi:10.1214/25-BJPS645](https://doi.org/10.1214/25-BJPS645)
  — *Basis for the `penalty` argument in `barma()`.*
- **Costa, E., Cribari-Neto, F., & Scher, V. T. (2024).** “Test
  inferences and link function selection in dynamic beta modeling of
  seasonal hydro-environmental time series with temporary abnormal
  regimes.” *Journal of Hydrology*, 638, 131489.

------------------------------------------------------------------------

## Key Features

This package utilizes modern R development standards (S3 classes,
roxygen2 documentation) to provide a seamless user experience.

- **Unified Model Fitting:** A single core function, `barma()`, handles
  any combination of AR and MA lags, exogenous regressors (`xreg`), and
  ridge penalization (`penalty = TRUE`).
- **Estimation:** Implements Conditional Maximum Likelihood Estimation
  (CMLE) and Penalized CMLE (PCMLE). Users can select the numerical
  optimization backend (`BFGS`, `L-BFGS-B`, or `lbfgs`) utilizing
  **analytical gradients** (score vectors) for improved convergence and
  speed.
- **Correctness:** Fully implements the corrections to the Information
  Matrix and Score Vector detailed in the 2017 Erratum, for all
  supported link functions (`logit`, `cloglog`, `loglog`, `probit`).
- **Standard S3 Methods:** Works with standard R functions:
  - `summary()`: Detailed coefficients, significance levels (stars),
    optimization tracking, and penalization metrics.
  - `plot()`: Six-panel diagnostic grid (observed vs. fitted, residuals
    over time, ACF, PACF, Ljung-Box p-values, Monti p-values) via which
    = ‘all’, or a four-panel default view.
  - `forecast()`: Multi-step-ahead point forecasts.
  - `residuals()`: Extraction of Pearson, quantile, raw, or link-scale
    residuals.
  - `fitted()`: Extraction of fitted values.
- **Advanced Access:** Core mathematical functions (`loglik_barma()`,
  `score_vector_barma()`, `fim_barma()`) are exported for researchers
  who need direct access to gradients or model internals.

------------------------------------------------------------------------

## Repository Structure

The repository is structured as a standard R package for clarity and
reproducibility.

``` plaintext
.
├── R/                  # Source code for all R functions.
├── data/               # Processed data objects (e.g., brasilia_ts).
├── data-raw/           # Scripts and raw CSV files used to generate the datasets.
├── doc/                # Compiled HTML vignettes and extracted R code.
├── inst/               # Additional package files (CITATION, REFERENCES.bib).
├── man/                # R package documentation files (generated by roxygen2).
├── vignettes/          # R Markdown source files for vignettes.
├── DESCRIPTION         # Package metadata and dependencies.
├── NAMESPACE           # Manages the package's namespace (generated by roxygen2).
├── NEWS.md             # Package changelog.
├── LICENSE             # MIT License file.
├── betaARMA.Rproj      # RStudio project file.
└── README.md           # This file.
```

------------------------------------------------------------------------

## Installation

### CRAN (Recommended)

You can install the stable version of `betaARMA` directly from CRAN:

``` r
install.packages("betaARMA")
```

### Development Version

You can install the development version from GitHub using the remotes
package:

``` r
# install.packages("remotes")
remotes::install_github("Everton-da-Costa/betaARMA", dependencies = TRUE)
```

------------------------------------------------------------------------

## Usage and Vignettes

The package includes a detailed vignette demonstrating a complete
modeling workflow, including how to handle convergence failures using
the penalized estimator.

To explore the application of betaARMA modeling on empirical data, check
out the main vignette:

``` r
vignette("humidity_brasilia", package = "betaARMA")
```

### Application: Modeling Relative Humidity in Brasília

This vignette walks through a full analysis of monthly relative humidity
in Brasília, Brazil. It highlights a specific lag specification where
standard CMLE fails to converge. It demonstrates how applying the PCMLE
(ridge-penalized) estimator achieves convergence, allows for valid model
reduction, and ultimately produces a stable, well-specified model that
passes all residual diagnostic tests.

------------------------------------------------------------------------

## Getting Started

Here is a quick example of how to simulate data, fit a model, and
generate forecasts.

``` r
library(betaARMA)

# Example 1: Fit a BAR(1) model

# 1. Simulate data from a BAR(1) process
set.seed(2025)
y_sim_bar <- simu_barma(
  n = 250,
  alpha = 0.0,
  varphi = 0.6,
  phi = 25.0,
  link = "logit",
  freq = 12
)

# 2. Fit the model
fit_bar <- barma(y_sim_bar, ar = 1, link = "logit")

# 3. Inspect results
summary(fit_bar)
coef(fit_bar)

# 4. Diagnostic plots
plot(fit_bar)

# 5. Forecast the next 6 steps
pred_bar <- forecast(fit_bar, h = 6)
print(pred_bar)

# Example 2: Fit a BARMA(1,1) model using L-BFGS-B optimization

# 1. Simulate data from a BARMA(1,1) process
set.seed(2025)
y_sim_barma <- simu_barma(
    n = 250,
    alpha = 0.0,
    varphi = 0.6, 
    theta = 0.3, 
    phi = 25,
    link = "logit",
    freq = 12
)

# 2. Fit the model explicitly using L-BFGS-B
fit_barma <- barma(
    y_sim_barma, 
    ar = 1, 
    ma = 1, 
    link = "logit",
    optimization = list(method = "L-BFGS-B")
)

# 3. Inspect results
summary(fit_barma)
coef(fit_barma)

# 4. Forecast the next 6 steps
pred_barma <- forecast(fit_barma, h = 6)
print(pred_barma)

# Example 3: Ridge-penalized estimation (PCMLE)
# Useful when standard CMLE produces implausible estimates.

fit_penalized <- barma(y_sim_barma, ar = 1, ma = 1, link = "logit", penalty = TRUE)
summary(fit_penalized)
```

------------------------------------------------------------------------

## Citation

If you use this package in your research, please cite it as follows:

``` bibtex
@Manual{,
  title = {betaARMA: Beta Autoregressive Moving Average Models},
  author = {Everton da Costa, Francisco Cribari-Neto and Vinícius T. Scher},
  year = {2026},
  note = {R package version 1.2.0},
  doi = {10.32614/CRAN.package.betaARMA},
  url = {https://CRAN.R-project.org/package=betaARMA}
}
```

If you use the ridge penalization (`penalty = TRUE`), please also cite:

``` bibtex
@Article{Cribari_Costa_Fonseca2025,
  title   = {Numerical stability enhancements in beta autoregressive moving
             average model estimation},
  author  = {Francisco Cribari-Neto and Everton Costa and Rodney V. Fonseca},
  journal = {Brazilian Journal of Probability and Statistics},
  year    = {2025},
  volume  = {39},
  number  = {4},
  pages   = {410--437},
  doi     = {10.1214/25-BJPS645}
}
```

If your application involves hydro-environmental modeling, please
consider citing:

``` bibtex
@Article{Costa_Cribari_Scher2024,
  title     = {Test inferences and link function selection in dynamic beta 
               modeling of seasonal hydro-environmental time series with 
               temporary abnormal regimes},
  author    = {Costa, E. and Cribari-Neto, F. and Scher, V. T.},
  journal   = {Journal of Hydrology},
  year      = {2024},
  note      = {forthcoming}
}
```

------------------------------------------------------------------------

## Contributing

Contributions are welcome! If you find any issues or have suggestions
for improvements, please open an issue or submit a pull request.

## License

This project is licensed under the **MIT License**. See the `LICENSE`
file for details.

## Contact

For questions, suggestions, or issues related to the code, please
contact:

**Everton da Costa** 📧 <everto.cost@gmail.com>

------------------------------------------------------------------------

## Code of Conduct

Please note that the betaARMA project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
