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
#' @author
#' Original code by Fabio M. Bayer (bayer@ufsm.br).
#' Substantially modified and improved by Everton da Costa
#' (everto.cost@gmail.com).
#'
#' @param y A numeric time series with values in the open interval (0, 1).
#' @param link A string specifying the link function for the mean.
#' @param ar A numeric vector of autoregressive (AR) lags. Defaults to
#'   \code{integer(0)}, which omits the AR component entirely.
#' @param ma A numeric vector of moving average (MA) lags. Defaults to
#'   \code{integer(0)}, which omits the MA component entirely.
#' @param xreg An optional numeric matrix or data frame of exogenous variables.
#'
#' @importFrom stats lm.fit fitted residuals var
#'
#' @return
#' A named numeric vector containing the initial values for the model
#' parameters (`alpha`, `varphi`, `theta`, `beta`, `phi`), ready to be
#' used by an optimization routine.
#'
#' @export
start_values <- function(y, link, ar = integer(0), ma = integer(0), xreg = NA) {
  
  # Link functions
  # ------------------------------------------------------------------------- #
  link_structure <- make_link_structure(link)
  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv
  mu.eta  <- link_structure$mu.eta
  
  y_transformed <- linkfun(y)
  n_obs <- length(y)
  
  # Resolve lag vectors to integer(0) if absent
  ar_lags <- if (length(ar) > 0) ar else integer(0)
  ma_lags <- if (length(ma) > 0) ma else integer(0)
  
  # Named flags for readability throughout the function
  has_ar   <- length(ar_lags) > 0
  has_ma   <- length(ma_lags) > 0
  has_xreg <- !is.null(xreg) && !all(is.na(xreg))
  
  ar_order <- if (has_ar) max(ar_lags, 0L) else 0L
  ma_order <- if (has_ma) max(ma_lags, 0L) else 0L
  
  n_ar_params <- length(ar_lags)
  n_ma_params <- length(ma_lags)
  
  max_lag <- max(ar_order, ma_order)
  
  if (has_xreg) {
    xreg <- as.matrix(xreg)
  }
  
  names_varphi <- if (has_ar) paste0("varphi", ar_lags) else character(0)
  names_theta  <- if (has_ma) paste0("theta", ma_lags) else character(0)
  names_beta   <- if (has_xreg) colnames(xreg) else character(0)
  
  n_eff <- n_obs - max_lag
  
  # ======================================================================== #
  # Initialize empty containers (Single Point of Truth prep)
  # ======================================================================== #
  alpha_start  <- numeric(0)
  varphi_start <- numeric(0)
  theta_start  <- numeric(0)
  beta_start   <- numeric(0)
  phi_start    <- numeric(0)
  
  # ======================================================================== #
  # 1. BARMA initial values (has_ar, has_ma, !has_xreg)
  # ======================================================================== #
  if (has_ar && has_ma && !has_xreg) {
    P <- matrix(NA, nrow = n_eff, ncol = n_ar_params)
    for (i in 1:n_eff) {
      P[i, ] <- y_transformed[i + max_lag - ar_lags]
    }
    x_inter <- matrix(1, nrow = n_eff, ncol = 1)
    
    fit_start <- lm.fit(
      x = cbind(x_inter, P),
      y = y_transformed[(max_lag + 1):n_obs]
    )
    mqo <- fit_start$coef
    
    alpha_start  <- mqo[1]
    varphi_start <- mqo[2:(n_ar_params + 1)]
    theta_start  <- rep(0, n_ma_params)
    phi_start    <- .get_phi_start(fit_start, n_eff, linkinv, mu.eta)
    
  }
  
  # ======================================================================== #
  # 2. BAR initial values (has_ar, !has_ma, !has_xreg)
  # ======================================================================== #
  if (has_ar && !has_ma && !has_xreg) {
    P <- matrix(NA, nrow = n_eff, ncol = n_ar_params)
    for (i in 1:n_eff) {
      P[i, ] <- y_transformed[i + max_lag - ar_lags]
    }
    x_inter <- matrix(1, nrow = n_eff, ncol = 1)
    
    fit_start <- lm.fit(
      x = cbind(x_inter, P),
      y = y_transformed[(max_lag + 1):n_obs]
    )
    mqo <- fit_start$coef
    
    alpha_start  <- mqo[1]
    varphi_start <- mqo[2:(n_ar_params + 1)]
    phi_start    <- .get_phi_start(fit_start, n_eff, linkinv, mu.eta)
    
  }
  # ======================================================================== #
  # 3. BMA initial values (!has_ar, has_ma, !has_xreg)
  # ======================================================================== #
  if (!has_ar && has_ma && !has_xreg) {
    mean_y <- mean(y)
    alpha_start <- mean(y_transformed)
    theta_start <- rep(0, n_ma_params)
    phi_start   <- (mean_y * (1 - mean_y)) / var(y)
    
  }
  
  # ======================================================================== #
  # 4. BARMAX initial values (has_ar, has_ma, has_xreg)
  # ======================================================================== #
  if (has_ar && has_ma && has_xreg) {
    P <- matrix(NA, nrow = n_eff, ncol = n_ar_params)
    for (i in 1:n_eff) {
      P[i, ] <- y_transformed[i + max_lag - ar_lags]
    }
    x_inter <- matrix(1, nrow = n_eff, ncol = 1)
    x_reg_mat <- xreg[(max_lag + 1):n_obs, , drop = FALSE]
    
    fit_start <- lm.fit(
      x = cbind(x_inter, P, x_reg_mat),
      y = y_transformed[(max_lag + 1):n_obs]
    )
    mqo <- fit_start$coef
    
    alpha_start  <- mqo[1]
    varphi_start <- mqo[2:(n_ar_params + 1)]
    theta_start  <- rep(0, n_ma_params)
    beta_start   <- mqo[(n_ar_params + 2):length(mqo)]
    phi_start    <- .get_phi_start(fit_start, n_eff, linkinv, mu.eta)
    
  }
  
  # ======================================================================== #
  # 5. BARX initial values (has_ar, !has_ma, has_xreg)
  # ======================================================================== #
  if (has_ar && !has_ma && has_xreg) {
    P <- matrix(NA, nrow = n_eff, ncol = n_ar_params)
    for (i in 1:n_eff) {
      P[i, ] <- y_transformed[i + max_lag - ar_lags]
    }
    x_inter <- matrix(1, nrow = n_eff, ncol = 1)
    x_reg_mat <- xreg[(max_lag + 1):n_obs, , drop = FALSE]
    
    fit_start <- lm.fit(
      x = cbind(x_inter, P, x_reg_mat),
      y = y_transformed[(max_lag + 1):n_obs]
    )
    mqo <- fit_start$coef
    
    alpha_start  <- mqo[1]
    varphi_start <- mqo[2:(n_ar_params + 1)]
    beta_start   <- mqo[(n_ar_params + 2):length(mqo)]
    phi_start    <- .get_phi_start(fit_start, n_eff, linkinv, mu.eta)
    
  }
  
  # ======================================================================== #
  # 6. BMAX initial values (!has_ar, has_ma, has_xreg)
  # ======================================================================== #
  if (!has_ar && has_ma && has_xreg) {
    x_inter <- matrix(1, nrow = n_eff, ncol = 1)
    x_reg_mat <- xreg[(max_lag + 1):n_obs, , drop = FALSE]
    
    fit_start <- lm.fit(
      x = cbind(x_inter, x_reg_mat),
      y = y_transformed[(max_lag + 1):n_obs]
    )
    mqo <- fit_start$coef
    
    alpha_start <- mqo[1]
    theta_start <- rep(0, n_ma_params)
    beta_start  <- mqo[2:length(mqo)]
    phi_start   <- .get_phi_start(fit_start, n_eff, linkinv, mu.eta)
    
  }
  
  # ======================================================================== #
  # 7. Beta Regression (!has_ar, !has_ma, has_xreg)
  # ======================================================================== #
  if (!has_ar && !has_ma && has_xreg) {
    x_inter <- matrix(1, nrow = n_obs, ncol = 1)
    
    fit_start <- lm.fit(
      x = cbind(x_inter, xreg),
      y = y_transformed
    )
    mqo <- fit_start$coef
    
    alpha_start <- mqo[1]
    beta_start  <- mqo[2:length(mqo)]
    phi_start   <- .get_phi_start(fit_start, n_obs, linkinv, mu.eta)
    
  }
  
  # ======================================================================== #
  # 8. Intercept-only Model (!has_ar, !has_ma, !has_xreg)
  # ======================================================================== #
  if (!has_ar && !has_ma && !has_xreg) {
    mean_y <- mean(y)
    alpha_start <- mean(y_transformed)
    phi_start <- (mean_y * (1 - mean_y)) / var(y)
    
  }
  
  # ======================================================================== #
  # FINAL ASSEMBLY
  # ======================================================================== #
  start_value <- c(
    alpha_start,
    varphi_start,
    theta_start,
    beta_start,
    phi_start
  )
  
  names(start_value) <- c(
    "alpha",
    if (length(varphi_start) > 0) names_varphi else character(0),
    if (length(theta_start) > 0)  names_theta  else character(0),
    if (length(beta_start) > 0)   names_beta   else character(0),
    "phi"
  )
  
  return(start_value)
}

#' Internal Helper to Calculate Initial Phi
#'
#' @description
#' Calculates the initial value for the precision parameter phi based on the
#' residuals of an initial 'lm.fit' object. The formula is adapted and
#' empirically improved from Ferrari & Cribari-Neto (2004) — specifically,
#' the "-1" term in the original precision estimator is omitted, as Monte
#' Carlo simulations showed better optimization performance without it.
#'
#' @param fit An 'lm.fit' object.
#' @param n_eff The effective number of observations used in the fit.
#' @param linkinv The inverse link function.
#' @param mu.eta The derivative of the mean w.r.t. eta (d(mu)/d(eta)).
#'
#' @return A single numeric value for phi_start.
#' @keywords internal
.get_phi_start <- function(fit, n_eff, linkinv, mu.eta) {
  
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
  
  phi_start <- phi_start_aux / n1
  # ---
  
  return(phi_start)
}