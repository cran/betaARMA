# betaARMA: Beta Autoregressive Moving Average Models

[![CRAN status](https://www.r-pkg.org/badges/version/betaARMA)](https://CRAN.R-project.org/package=betaARMA)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/betaARMA)](https://CRAN.R-project.org/package=betaARMA)
[![DOI](https://img.shields.io/badge/DOI-10.32614%2FCRAN.package.betaARMA-blue)](https://doi.org/10.32614/CRAN.package.betaARMA)
[![R-hub](https://github.com/Everton-da-Costa/betaARMA/actions/workflows/rhub.yaml/badge.svg)](https://github.com/Everton-da-Costa/betaARMA/actions/workflows/rhub.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the R package **`betaARMA`**, a comprehensive toolkit for fitting, forecasting, and simulating Beta Autoregressive Moving Average models. It provides a unified workflow for modeling time series data bounded on the (0, 1) interval, such as rates, proportions, and indices.

---

## What's New in v1.1.0
* **Advanced Internals Exported:** The core mathematical functions (`loglik_barma()`, `score_vector_barma()`, and `fim_barma()`) are now exported, allowing researchers and advanced users to build custom extensions or extract gradients directly.
* **Unified Parameter Architecture:** The underlying codebase has been heavily refactored for performance and stability. All functions now strictly adhere to the unified parameter vector order: `(alpha, varphi, theta, beta, phi)`.
* **Link Function Fixes:** Fixed a bug in the score vector derivation to ensure mathematically accurate inference for all non-logit models (`cloglog`, `loglog`, and `probit`).

---

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

---

## Project Motivation

Modeling time series data bounded within the unit interval $(0, 1)$ presents unique challenges. Standard Gaussian methods (like ARIMA) are often inappropriate because they do not respect the natural boundaries of the data, potentially leading to fitted values or forecasts outside the admissible range.

The **$\beta\text{ARMA}$ model** addresses this by assuming the conditional distribution of the variable follows a Beta law. While the theoretical foundations exist, there has been a need for a modern, robust R package that:

1.  **Unifies the workflow** for $\beta$AR, $\beta$MA, and $\beta$ARMA specifications.
2.  **Ensures stability** through numerical optimization with analytic gradients.
3.  **Provides a standard API** consistent with popular time series tools (like `forecast`).

This project aims to fill that gap, serving as a go-to resource for hydrologists, economists, and data scientists working with bounded data.

---

## Foundational Literature

This package implements the methodology established in the following key publications. The original code foundation was developed by Fabio M. Bayer and has been substantially optimized and refactored for this package.

* **Rocha, A. V., & Cribari-Neto, F. (2009).** "Beta autoregressive moving average models." *TEST*, 18(3), 529-545. [doi:10.1007/s11749-008-0112-z](https://doi.org/10.1007/s11749-008-0112-z)
* **Rocha, A. V., & Cribari-Neto, F. (2017).** "Erratum to: Beta autoregressive moving average models." *TEST*, 26(2), 451-459. [doi:10.1007/s11749-017-0528-4](https://doi.org/10.1007/s11749-017-0528-4)

### Journal Quality Metrics (TEST)

[![SCImago Journal & Country Rank](https://www.scimagojr.com/journal_img.php?id=14882)](https://www.scimagojr.com/journalsearch.php?q=14882&tip=sid)

* **SJR (2024):** 0.505 (Q2)
* **H-Index:** 52

---

## Key Features

This package utilizes modern R development standards (S3 classes, roxygen2 documentation) to provide a seamless user experience.

* **Unified Model Fitting:** A single core function, `barma()`, handles any combination of AR and MA lags, as well as exogenous regressors (`xreg`).
* **Estimation:** Implements Conditional Maximum Likelihood Estimation (CMLE) using the BFGS algorithm with **analytical gradients** (score vectors) for improved convergence and speed.
* **Correctness:** Fully implements the corrections to the Information Matrix and Score Vector detailed in the 2017 Erratum.
* **Standard S3 Methods:** Works with standard R functions:
    * `summary()`: Detailed coefficients, standard errors, and significance tests.
    * `forecast()`: Multi-step-ahead point forecasts.
    * `residuals()`: Extraction of standardized residuals.
    * `fitted()`: Extraction of fitted values.

---

## Repository Structure

The repository is structured as a standard R package for clarity and reproducibility.

```plaintext
.
├── R/                  # Source code for all R functions.
├── man/                # R package documentation files (generated by roxygen2).
├── tests/              # Unit tests (using testthat).
├── validation/         # Validation scripts comparing results to literature.
├── DESCRIPTION         # Package metadata and dependencies.
├── NAMESPACE           # Manages the package's namespace (generated by roxygen2).
├── LICENSE             # MIT License file.
└── README.md           # This file.
```

---

## Installation

### CRAN (Recommended)
You can install the stable version of `betaARMA` directly from CRAN:

```R
install.packages("betaARMA")
```

### Development Version
You can install the development version from GitHub using the remotes package:

```R
# install.packages("remotes")
remotes::install_github("Everton-da-Costa/betaARMA", dependencies = TRUE)
```

---

## Getting Started

Here is a quick example of how to simulate data, fit a model, and generate forecasts.

```R
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

# 4. Forecast the next 6 steps
pred_bar <- forecast(fit_bar, h = 6)
print(pred_bar)


# Example 2: Fit a BARMA(1,1) model

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


# 2. Fit the model
fit_barma <- barma(y_sim_barma, ar = 1, ma = 1, link = "logit")

# 3. Inspect results
summary(fit_barma)
coef(fit_barma)

# 4. Forecast the next 6 steps
pred_barma <- forecast(fit_barma, h = 6)
print(pred_barma)
```

---

## Citation

If you use this package in your research, please cite it as follows:

```bibtex
@Manual{,
  title = {betaARMA: Beta Autoregressive Moving Average Models},
  author = {Everton da Costa, Francisco Cribari-Neto and Vinícius T. Scher},
  year = {2026},
  note = {R package version 1.0.1},
  doi = {10.32614/CRAN.package.betaARMA},
  url = {https://CRAN.R-project.org/package=betaARMA}
}
```
## Contributing

Contributions are welcome! If you find any issues or have suggestions for improvements, please open an issue or submit a pull request.

## License

This project is licensed under the **MIT License**. See the `LICENSE` file for details.

## Contact

For questions, suggestions, or issues related to the code, please contact:

**Everton da Costa**
📧 <everto.cost@gmail.com>

---

## Code of Conduct

Please note that the betaARMA project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.