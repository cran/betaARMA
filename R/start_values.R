#' Generate Initial Values for BARMA Model Estimation
#'
#' @description
#' This function calculates reasonable starting values for the parameters of
#' various Beta Autoregressive Moving Average (BARMA) models. The method is
#' based on the approach proposed by Ferrari & Cribari-Neto (2004) for
#' beta regression, adapted here for the time series context.
#'
#' @details
#' The function computes initial values by fitting a linear model
#' (`lm.fit`) to the link-transformed response variable `g(y)`.
#' This provides a computationally cheap and stable way to initialize the
#' main optimization algorithm.
#'
#' The specific procedure depends on the model's structure:
#' \itemize{
#'   \item For models with autoregressive (AR) terms, the lagged values of
#'     the transformed series `g(y)` are used as predictors to estimate
#'     the initial `alpha` (intercept) and `varphi` (AR) coefficients.
#'   \item If exogenous variables `X` are included, they are added as
#'     predictors in the linear model to obtain initial `beta` values.
#'   \item Moving average (MA) coefficients (`theta`) are initialized to zero,
#'     a standard practice in ARMA model estimation.
#'   \item The initial value for the precision parameter `phi` is derived
#'     from the variance of the residuals of the initial linear fit,
#'     following the methodology from Ferrari & Cribari-Neto (2004).
#'   \item For a pure BMA model (no AR or X components), a simpler method is
#'     used where `alpha` is the mean of `g(y)` and `phi` is based on the
#'     unconditional variance of `y`.
#' }
#'
#' This function is called internally by the main model fitting function
#' but is exported for standalone use and inspection.
#'
#' @author
#' Original code by Fabio M. Bayer (bayer@ufsm.br).
#' Substantially modified and improved by Everton da Costa
#' (everto.cost@gmail.com).
#'
#' @references
#' Ferrari, S. L. P., & Cribari-Neto, F. (2004). Beta regression for
#' modelling rates and proportions. *Journal of Applied Statistics*,
#' 31(7), 799-818. <doi:10.1080/0266476042000214501>
#'
#' @param y A numeric time series with values in the open interval (0, 1).
#' @param link A string specifying the link function for the mean, such as
#'   `"logit"`, `"probit"`, or `"cloglog"`.
#' @param ar A numeric vector of autoregressive (AR) lags. Defaults to `NA`
#'   for models without an AR component.
#' @param ma A numeric vector of moving average (MA) lags. Defaults to `NA`
#'   for models without an MA component.
#' @param X An optional numeric matrix or data frame of exogenous
#'   variables (regressors).
#'
#' @importFrom stats lm.fit fitted residuals var
#'
#' @return
#' A named numeric vector containing the initial values for the model
#' parameters (`alpha`, `varphi`, `theta`, `phi`, `beta`), ready to be
#' used by an optimization routine. Returns `NULL` if the model
#' specification is not recognized.
#'
#' @export
start_values <- function(y, link,
                         ar = NA, ma = NA, X = NA) {
  # from:
  #       Beta Regression for Modelling Rates and Proportions
  #       Silvia Ferrari & Francisco Cribari-Neto
  #       p. 805
  
  # Link functions
  # ----------------------------------------------------------------------- #
  link_structure <- make_link_structure(link)
  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv
  mu.eta  <- link_structure$mu.eta
  
  ynew <- linkfun(y)
  n_obs <- length(y)
  
  # Determine model components presence
  has_ar <- !is.null(ar) && !any(is.na(ar)) && length(ar) > 0
  has_ma <- !is.null(ma) && !any(is.na(ma)) && length(ma) > 0
  has_X <- !is.null(X) && !all(is.na(X)) &&
    (is.matrix(X) || is.data.frame(X))
  
  # Use consistent naming convention
  ar_lags <- if (has_ar) ar else integer(0)
  ma_lags <- if (has_ma) ma else integer(0)
  
  ar_order <- ifelse(has_ar, max(ar_lags, 0L), 0L)
  ma_order <- ifelse(has_ma, max(ma_lags, 0L), 0L)
  
  n_ar_params <- length(ar_lags)
  n_ma_params <- length(ma_lags)
  
  max_lag <- max(ar_order, ma_order)
  
  if (has_X) {
    X <- as.matrix(X) # Ensure X is a matrix if it's used
  }
  
  names_varphi <- if (has_ar) paste0("varphi", ar_lags) else character(0)
  names_theta <- if (has_ma) paste0("theta", ma_lags) else character(0)
  names_beta <- if (has_X) colnames(X) else character(0)
  
  n_eff <- n_obs - max_lag # Effective number of observations
  
  # ========================================================================= #
  # BARMA initial values (has_ar, has_ma, !has_X)
  # ========================================================================= #
  if (has_ar && has_ma && !has_X) {
    
    P <- matrix(NA, nrow = n_eff, ncol = n_ar_params)
    for (i in 1:n_eff) P[i, ] <- ynew[i + max_lag - ar_lags]
    
    x_inter <- matrix(1, nrow = n_eff, ncol = 1)
    x_start <- cbind(x_inter, P)
    y_start <- ynew[(max_lag + 1):n_obs]
    
    fit_start  <- lm.fit(x = x_start, y = y_start)
    mqo <- fit_start$coef
    
    alpha_start <- mqo[1]
    varphi_start <- mqo[2:(n_ar_params + 1)] 
    
    phi_start <- ._get_phi_start(
      fit = fit_start,
      n_eff = n_eff,
      linkinv = linkinv,
      mu.eta = mu.eta
    )
    
    theta_start <- rep(0, n_ma_params)
    
    start_value <- c(alpha_start, varphi_start, theta_start, phi_start)
    names(start_value) <- c("alpha", names_varphi, names_theta, "phi")
    
    return(start_value)
  }
  
  # ============================================================================
  # BAR initial values (has_ar, !has_ma, !has_X)
  # ============================================================================
  if (has_ar && !has_ma && !has_X) {
    
    P <- matrix(NA, nrow = n_eff, ncol = n_ar_params)
    for (i in 1:n_eff) P[i, ] <- ynew[i + max_lag - ar_lags]
    
    x_inter <- matrix(1, nrow = n_eff, ncol = 1)
    x_start <- cbind(x_inter, P)
    y_start <- ynew[(max_lag + 1):n_obs]
    
    fit_start  <- lm.fit(x = x_start, y = y_start)
    mqo <- fit_start$coef
    
    alpha_start <- mqo[1]
    varphi_start <- mqo[2:(n_ar_params + 1)]
    
    phi_start <- ._get_phi_start(
      fit = fit_start,
      n_eff = n_eff,
      linkinv = linkinv,
      mu.eta = mu.eta
    )
    
    start_value <- c(alpha_start, varphi_start, phi_start)
    names(start_value) <- c("alpha", names_varphi, "phi")
    
    return(start_value)
  }
  
  # ============================================================================
  # BMA initial values (!has_ar, has_ma, !has_X)
  # ============================================================================
  if (!has_ar && has_ma && !has_X) {
    
    mean_y <- mean(y)
    alpha_start <- mean(ynew)
    theta_start <- rep(0, n_ma_params)
    
    phi_start <- (mean_y * (1 - mean_y)) / var(y)
    
    start_value <- c(alpha_start, theta_start, phi_start)
    names(start_value) <- c("alpha", names_theta, "phi")
    
    return(start_value)
  }
  
  # ========================================================================= #
  # BARMAX initial values (has_ar, has_ma, has_X)
  # ========================================================================= #
  if (has_ar && has_ma && has_X) {
    
    P <- matrix(NA, nrow = n_eff, ncol = n_ar_params)
    for (i in 1:n_eff) P[i, ] <- ynew[i + max_lag - ar_lags]
    
    x_inter <- matrix(1, nrow = n_eff, ncol = 1)
    x_reg <- X[(max_lag + 1):n_obs, , drop = FALSE]
    x_start <- cbind(x_inter, P, x_reg)
    y_start <- ynew[(max_lag + 1):n_obs]
    
    fit_start <- lm.fit(x = x_start, y = y_start)
    mqo <- c(fit_start$coef) 
    
    phi_start <- ._get_phi_start(
      fit = fit_start,
      n_eff = n_eff,
      linkinv = linkinv,
      mu.eta = mu.eta
    )
    
    alpha_start <- mqo[1]
    varphi_start <- mqo[2:(n_ar_params + 1)]
    theta_start <- rep(0, n_ma_params)
    beta_start <- mqo[(n_ar_params + 2):length(mqo)]
    
    start_value <-
      c(alpha_start, varphi_start, theta_start, phi_start, beta_start)
    
    names(start_value) <-
      c("alpha", names_varphi, names_theta, "phi", names_beta)
    
    return(start_value)
  }
  
  # ========================================================================= #
  # BARX initial values (has_ar, !has_ma, has_X)
  # ========================================================================= #
  if (has_ar && !has_ma && has_X) {
    
    P <- matrix(NA, nrow = n_eff, ncol = n_ar_params)
    for (i in 1:n_eff) P[i, ] <- ynew[i + max_lag - ar_lags]
    
    x_inter <- matrix(1, nrow = n_eff, ncol = 1)
    x_reg <- X[(max_lag + 1):n_obs, , drop = FALSE]
    x_start <- cbind(x_inter, P, x_reg)
    y_start <- ynew[(max_lag + 1):n_obs]
    
    fit_start <- lm.fit(x = x_start, y = y_start)
    mqo <- c(fit_start$coef) 
    
    phi_start <- ._get_phi_start(
      fit = fit_start,
      n_eff = n_eff,
      linkinv = linkinv,
      mu.eta = mu.eta
    )
    
    alpha_start <- mqo[1]
    varphi_start <- mqo[2:(n_ar_params + 1)]
    beta_start <- mqo[(n_ar_params + 2):length(mqo)]
    
    start_value <- c(alpha_start, varphi_start, phi_start, beta_start)
    names(start_value) <- c("alpha", names_varphi, "phi", names_beta)
    
    return(start_value)
  }
  
  # ========================================================================= #
  # BMAX initial values (!has_ar, has_ma, has_X)
  # ========================================================================= #
  if (!has_ar && has_ma && has_X) {
    
    x_inter <- matrix(1, nrow = n_eff, ncol = 1)
    x_reg <- X[(max_lag + 1):n_obs, , drop = FALSE]
    x_start <- cbind(x_inter, x_reg)
    y_start <- ynew[(max_lag + 1):n_obs]
    
    fit_start <- lm.fit(x = x_start, y = y_start)
    mqo <- fit_start$coef
    
    phi_start <- ._get_phi_start(
      fit = fit_start,
      n_eff = n_eff,
      linkinv = linkinv,
      mu.eta = mu.eta
    )
    
    alpha_start <- mqo[1]
    theta_start <- rep(0, n_ma_params)
    beta_start <- mqo[2:length(mqo)]
    
    start_value <- c(alpha_start, theta_start, phi_start, beta_start)
    names(start_value) <- c("alpha", names_theta, "phi", names_beta)
    
    return(start_value)
  }
  
  # ========================================================================= #
  # Beta Regression (!has_ar, !has_ma, has_X)
  # ========================================================================= #
  if (!has_ar && !has_ma && has_X) {
    
    x_inter <- matrix(1, nrow = n_obs, ncol = 1) # max_lag is 0
    x_start <- cbind(x_inter, X)
    y_start <- ynew 
    
    fit_start <- lm.fit(x = x_start, y = y_start)
    mqo <- fit_start$coef
    
    phi_start <- ._get_phi_start(
      fit = fit_start,
      n_eff = n_obs,
      linkinv = linkinv,
      mu.eta = mu.eta
    )
    
    alpha_start <- mqo[1]
    beta_start <- mqo[2:length(mqo)]
    
    start_value <- c(alpha_start, phi_start, beta_start)
    names(start_value) <- c("alpha", "phi", names_beta)
    
    return(start_value)
  }
  
  # ========================================================================= #
  # Intercept-only Model (!has_ar, !has_ma, !has_X)
  # ========================================================================= #
  if (!has_ar && !has_ma && !has_X) {
    
    mean_y <- mean(y)
    alpha_start <- mean(ynew)
    phi_start <- (mean_y * (1 - mean_y)) / var(y) 
    
    start_value <- c(alpha_start, phi_start)
    names(start_value) <- c("alpha", "phi")
    
    return(start_value)
  }
  
  warning("No matching model configuration found for initial values.")
  return(NULL)
}

#' Internal Helper to Calculate Initial Phi
#'
#' @description
#' Calculates the initial value for the precision parameter phi based on the
#' residuals of an initial 'lm.fit' object, following Ferrari & 
#' Cribari-Neto (2004).
#'
#' This helper uses the *exact* mathematical logic from the original
#' 'start_values.R' file to ensure identical numerical output.
#'
#' @param fit An 'lm.fit' object.
#' @param n_eff The effective number of observations used in the fit.
#' @param linkinv The inverse link function.
#' @param mu.eta The derivative of the mean w.r.t. eta (d(mu)/d(eta)).
#'
#' @return A single numeric value for phi_start.
#' @keywords internal
._get_phi_start <- function(fit, n_eff, linkinv, mu.eta) {
  
  mqo <- fit$coef
  k <- length(mqo)
  n1 <- n_eff
  
  y_hat_fit_start <- fitted(fit)
  mean_fit_start <- linkinv(y_hat_fit_start)
  
  
  linkfun_deriv_aux <- mu.eta(eta = y_hat_fit_start)
  linkfun_deriv <- 1 / linkfun_deriv_aux
  
  er <- residuals(fit)
  
  # ---
  sigma2 <- sum(er^2) / ((n1 - k) * linkfun_deriv^2)
  phi_start_aux <- sum(mean_fit_start * (1 - mean_fit_start) / sigma2)
  
  # NOTE: The formula in Ferrari & Cribari-Neto (2004) implies a 
  # slightly different estimator (`phi_start_aux - 1` / n1).
  # However, extensive simulation studies showed that omitting the
  # '-1' yields better initial values with lower bias and RMSE
  # for the final parameter estimates.
  phi_start <- phi_start_aux / n1
  # ---
  
  return(phi_start)
}