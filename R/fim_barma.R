#' Fisher Information Matrix for BARMA Models
#'
#' @title
#' Compute Fisher Information Matrix for Beta Autoregressive
#' Moving Average Models
#'
#' @description
#' Computes the observed Fisher Information Matrix (FIM) of a Beta
#' Autoregressive Moving Average (BARMA) model. This function also
#' efficiently returns auxiliary values like fitted values and residuals.
#' 
#' This function is designed for users who:
#' \itemize{
#'   \item Compute standard errors of parameter estimates
#'   \item Construct confidence intervals
#'   \item Perform hypothesis tests
#'   \item Verify theoretical properties
#'   \item Conduct simulation studies
#' }
#'
#' @details
#' The Fisher Information Matrix is computed from the outer product of
#' score vectors at the MLE, which provides the observed FIM. The FIM
#' is used to obtain standard errors via \code{sqrt(diag(solve(FIM)))}.
#'
#' **Non-Diagonal Structure**: Unlike generalized linear models, the FIM
#' is not block-diagonal due to the coupling introduced by the ARMA
#' dynamics. This is documented in the Rocha & Cribari-Neto (2009) paper.
#'
#' **Important**: This function implements the corrections from the 2017
#' Erratum (Rocha & Cribari-Neto, 2017) for moving average components.
#' See References section for details.
#'
#' **Parameter Order**: Parameters should be supplied in the order:
#' alpha, varphi (AR), theta (MA), phi, beta (regressors).
#' This matches the parameter order used by \code{\link{barma}}.
#'
#' @param y
#' A numeric vector representing the time series data, in (0, 1).
#'
#' @param ar
#' A numeric vector specifying the autoregressive (AR) lags
#' (e.g., \code{c(1, 2)}). Can be \code{NA} or \code{NULL} if no AR component.
#'
#' @param ma
#' A numeric vector specifying the moving average (MA) lags
#' (e.g., \code{1}). Can be \code{NA} or \code{NULL} if no MA component.
#'
#' @param alpha
#' The intercept term (numeric scalar).
#'
#' @param varphi
#' A numeric vector of autoregressive (AR) parameters.
#' Use \code{numeric(0)} or empty vector if no AR component.
#' Length must match \code{ar} specification.
#'
#' @param theta
#' A numeric vector of moving average (MA) parameters.
#' Use \code{numeric(0)} or empty vector if no MA component.
#' Length must match \code{ma} specification.
#'
#' @param phi
#' The precision parameter of the BARMA model (must be positive and finite).
#' Larger values indicate less variance for a given mean.
#'
#' @param link
#' A character string specifying the link function:
#' \code{"logit"} (default), \code{"probit"}, or \code{"cloglog"}.
#'
#' @param xreg
#' A matrix or data frame of static regressors (optional).
#' Must have the same number of rows as length of \code{y}.
#'
#' @param beta
#' A numeric vector of regression coefficients for \code{xreg} (optional).
#' Length must match number of columns in \code{xreg}.
#'
#' @return
#' A list containing:
#' \item{fisher_info_mat}{The Fisher Information Matrix (numeric matrix).
#'   Dimensions: \code{(k+p+q+2) x (k+p+q+2)} where k is number of regressors,
#'   p is AR order, q is MA order. The matrix is symmetric and should be
#'   positive definite at the MLE.}
#' \item{fitted_ts}{The fitted values (conditional mean) as a \code{ts} object,
#'   with the same time index as the input \code{y}. Values for the
#'   burn-in period (first max_lag observations) are \code{NA}.}
#' \item{muhat_effective}{The fitted conditional means (numeric vector)
#'   excluding the burn-in period.}
#' \item{etahat_full}{The estimated linear predictor values (numeric vector,
#'   length = length(y)), with \code{NA} for burn-in period.}
#' \item{errorhat_full}{The estimated errors on the predictor scale
#'   (numeric vector, length = length(y)), zero-padded for burn-in period.}
#'
#' @references
#' Rocha, A.V., & Cribari-Neto, F. (2009). Beta autoregressive moving
#' average models. \emph{TEST}, 18(3), 529-545.
#' \doi{10.1007/s11749-008-0112-z}
#'
#' Rocha, A.V., & Cribari-Neto, F. (2017). Erratum to: Beta autoregressive
#' moving average models. \emph{TEST}, 26, 451-459.
#' \doi{10.1007/s11749-017-0528-4}
#'
#' @seealso
#' \code{\link{barma}} for model fitting,
#' \code{\link{loglik_barma}} for log-likelihood computation,
#' \code{\link{score_vector_barma}} for score vector (gradient)
#'
#' @examples
#' \donttest{
#'   # Example 1: Fisher Information Matrix for a BAR(1) model
#'   set.seed(2025)
#'   y_sim_bar <- simu_barma(
#'     n = 250,
#'     alpha = 0.0,
#'     varphi = 0.6,
#'     phi = 25.0,
#'     link = "logit",
#'     freq = 12
#'   )
#'
#'   result_bar <- fim_barma(
#'     y = y_sim_bar,
#'     ar = 1,
#'     ma = NA,
#'     alpha = 0.0,
#'     varphi = 0.6,
#'     theta = numeric(0),
#'     phi = 25.0,
#'     link = "logit"
#'   )
#'
#'   # Check positive definiteness
#'   fim <- result_bar$fisher_info_mat
#'   all(eigen(fim)$values > 0)
#'
#'   # Standard errors from inverse of FIM
#'   sqrt(diag(solve(fim)))
#'
#'   # Example 2: Fisher Information Matrix for a BARMA(1,1) model
#'   set.seed(2025)
#'   y_sim_barma <- simu_barma(
#'     n = 250,
#'     alpha = 0.0,
#'     varphi = 0.6,
#'     theta = 0.3,
#'     phi = 25.0,
#'     link = "logit",
#'     freq = 12
#'   )
#'
#'   result_barma <- fim_barma(
#'     y = y_sim_barma,
#'     ar = 1,
#'     ma = 1,
#'     alpha = 0.0,
#'     varphi = 0.6,
#'     theta = 0.3,
#'     phi = 25.0,
#'     link = "logit"
#'   )
#'
#'   # Standard errors
#'   sqrt(diag(solve(result_barma$fisher_info_mat)))
#' }
#'
#' @export
fim_barma <- function(
    y,
    ar,
    ma,
    alpha,
    varphi,
    theta,
    phi,
    link,
    xreg = NULL,
    beta = NULL
) {
  
  # --------------------------------------------------------------------------
  # 1. SETUP LINK FUNCTIONS AND MODEL STRUCTURE
  # --------------------------------------------------------------------------
  
  # Get link function structures
  link_structure <- make_link_structure(link)
  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv
  mu_eta_fun <- link_structure$mu.eta
  
  # Transform response using link function
  y_transformed <- linkfun(y)
  
  # Check for presence of AR, MA, and external regressors
  has_ar <- !is.null(ar) && !any(is.na(ar)) && length(ar) > 0
  has_ma <- !is.null(ma) && !any(is.na(ma)) && length(ma) > 0
  has_xreg <- !is.null(xreg)
  
  # Handle empty lag cases
  ar_lags <- if (has_ar) ar else integer(0)
  ma_lags <- if (has_ma) ma else integer(0)
  
  n_ar_params <- length(ar_lags)
  n_ma_params <- length(ma_lags)
  
  # Setup Regressors ---
  if (has_xreg) {
    if (!is.matrix(xreg)) xreg <- as.matrix(xreg)
    n_beta_params <- ncol(xreg)
    # Pre-compute X * beta for efficiency (vectorized operation)
    xb <- as.vector(xreg %*% beta)
  } else {
    n_beta_params <- 0
    xb <- numeric(length(y))
  }
  
  # Determine maximum lag for burn-in period
  ar_order <- if (has_ar) max(ar_lags) else 0L
  ma_order <- if (has_ma) max(ma_lags) else 0L
  max_lag  <- max(ar_order, ma_order)
  
  n_obs <- length(y)
  
  # --------------------------------------------------------------------------
  # 2. RECURSIVE CALCULATION OF PREDICTOR, ERROR, AND DERIVATIVES
  # --------------------------------------------------------------------------
  
  # Initialize containers for linear predictor and errors
  error   <- rep(0, n_obs)
  eta     <- rep(NA_real_, n_obs)
  
  # Initialize derivative matrices for chain rule calculations
  d_eta_d_alpha  <- rep(0, n_obs)
  
  d_eta_d_varphi <- if (has_ar) {
    matrix(0, nrow = n_obs, ncol = n_ar_params)
  } else {
    matrix(0, nrow = n_obs, ncol = 0)
  }
  
  d_eta_d_theta  <- if (has_ma) {
    matrix(0, nrow = n_obs, ncol = n_ma_params)
  } else {
    matrix(0, nrow = n_obs, ncol = 0)
  }
  
  d_eta_d_beta <- if (has_xreg) {
    matrix(0, nrow = n_obs, ncol = n_beta_params)
  } else {
    matrix(0, nrow = n_obs, ncol = 0)
  }
  
  # Recursively compute predictor, errors, and their derivatives
  for (t in (max_lag + 1):n_obs) {
    
    # Part A: Compute Linear Predictor and Score Error ---
    eta[t] <- alpha + xb[t]
    
    if (has_ar) {
      prev_terms <- y_transformed[t - ar_lags] - xb[t - ar_lags]
      eta[t] <- eta[t] + drop(crossprod(varphi, prev_terms))
    }
    if (has_ma) {
      eta[t] <- eta[t] + drop(crossprod(theta, error[t - ma_lags]))
    }
    
    # Compute prediction error (on transformed scale)
    error[t] <- y_transformed[t] - eta[t]
    
    # Part B: Compute Recursive Derivatives ---
    
    # 1. Derivative w.r.t. alpha
    d_eta_d_alpha[t] <- 1
    if (has_ma) {
      ma_adj <- drop(crossprod(theta, d_eta_d_alpha[t - ma_lags]))
      d_eta_d_alpha[t] <- 1 - ma_adj
    }
    
    # 2. Derivative w.r.t. AR (varphi)
    if (has_ar) {
      d_eta_d_varphi[t, ] <- y_transformed[t - ar_lags] - xb[t - ar_lags]
      if (has_ma) {
        ma_adj <- drop(
          crossprod(theta, d_eta_d_varphi[t - ma_lags, , drop = FALSE])
        )
        d_eta_d_varphi[t, ] <- d_eta_d_varphi[t, ] - ma_adj
      }
    }
    
    # 3. Derivative w.r.t. MA (theta)
    if (has_ma) {
      d_eta_d_theta[t, ] <- error[t - ma_lags]
      ma_adj <- drop(
        crossprod(theta, d_eta_d_theta[t - ma_lags, , drop = FALSE])
      )
      d_eta_d_theta[t, ] <- d_eta_d_theta[t, ] - ma_adj
    }
    
    # 4. Derivative w.r.t. beta (Regressors)
    if (has_xreg) {
      base_grad <- xreg[t, ]
      if (has_ar) {
        ar_adj <- drop(
          crossprod(varphi, xreg[t - ar_lags, , drop = FALSE])
        )
        base_grad <- base_grad - ar_adj
      }
      
      d_eta_d_beta[t, ] <- base_grad
      
      if (has_ma) {
        ma_adj <- drop(
          crossprod(theta, d_eta_d_beta[t - ma_lags, , drop = FALSE])
        )
        d_eta_d_beta[t, ] <- d_eta_d_beta[t, ] - ma_adj
      }
    }
  }
  
  # --------------------------------------------------------------------------
  # 3. GET EFFECTIVE (NON-NA) OBSERVATIONS
  # --------------------------------------------------------------------------
  
  # Extract observations used in likelihood (excluding burn-in period)
  idx_effective <- (max_lag + 1):n_obs
  eta_eff <- eta[idx_effective]
  mu_eff  <- linkinv(eta = eta_eff)
  
  # Numerical stability check
  if (any(mu_eff <= 0 | mu_eff >= 1 | !is.finite(mu_eff))) {
    warning("mu_effective out of bounds during FIM calculation.")
    n_params <- 1 + n_ar_params + n_ma_params + n_beta_params + 1
    return(list(
      fisher_info_mat = matrix(NA_real_, nrow = n_params, ncol = n_params),
      fitted_ts = ts(rep(NA_real_, n_obs), start = start(y), 
                     frequency = frequency(y)),
      muhat_effective = rep(NA_real_, length(idx_effective)),
      etahat_full = eta,
      errorhat_full = error
    ))
  }
  
  # Extract effective derivative vectors/matrices
  s_vec <- d_eta_d_alpha[idx_effective]
  P_mat <- d_eta_d_varphi[idx_effective, , drop = FALSE]
  R_mat <- d_eta_d_theta[idx_effective, , drop = FALSE]
  M_mat <- d_eta_d_beta[idx_effective, , drop = FALSE]
  
  # --------------------------------------------------------------------------
  # 4. CALCULATE FIM COMPONENT VECTORS (Efficiently)
  # --------------------------------------------------------------------------
  
  # Trigamma functions: derivatives of digamma
  trigamma_p <- trigamma(mu_eff * phi)
  trigamma_q <- trigamma((1 - mu_eff) * phi)
  
  # d(mu)/d(eta) - derivative of inverse link function
  mu_eta_val <- mu_eta_fun(eta = eta_eff)
  
  # Weight Vector ---
  # W = phi * {psi'(mu*phi) + psi'((1-mu)*phi)} * (mu_eta)^2
  w_t_vec <- phi * (trigamma_p + trigamma_q) * mu_eta_val^2
  
  # Vector c ---
  # c = phi * {psi'(mu*phi)*mu - psi'((1-mu)*phi)*(1-mu)}
  c_t_vec <- phi * (trigamma_p * mu_eff - trigamma_q * (1 - mu_eff))
  
  # Vector d ---
  # d = psi'(mu*phi)*mu^2 + psi'((1-mu)*phi)*(1-mu)^2 - psi'(phi)
  d_t_vec <- trigamma_p * mu_eff^2 +
    trigamma_q * (1 - mu_eff)^2 -
    trigamma(phi)
  
  # --------------------------------------------------------------------------
  # 5. CALCULATE FIM BLOCKS
  # --------------------------------------------------------------------------
  
  # Helper to compute weighted cross-product
  calc_block <- function(X, Y, w, scale = 1) {
    if (length(X) == 0 || length(Y) == 0) {
      return(matrix(numeric(0), 0, 0))
    }
    
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    
    if (length(w) != nrow(Y)) {
      stop("Weight vector length must match number of rows in Y.")
    }
    
    result <- crossprod(X, Y * w) * scale
    return(result)
  }
  
  # Blocks scaling with phi ---
  K_aa <- calc_block(s_vec, s_vec, w_t_vec, scale = phi)
  
  K_ap <- if (has_ar) {
    calc_block(s_vec, P_mat, w_t_vec, scale = phi)
  } else {
    numeric(0)
  }
  
  K_at <- if (has_ma) {
    calc_block(s_vec, R_mat, w_t_vec, scale = phi)
  } else {
    numeric(0)
  }
  
  K_ab <- if (has_xreg) {
    calc_block(s_vec, M_mat, w_t_vec, scale = phi)
  } else {
    numeric(0)
  }
  
  K_pp <- if (has_ar) {
    calc_block(P_mat, P_mat, w_t_vec, scale = phi)
  } else {
    matrix(0, 0, 0)
  }
  
  K_pt <- if (has_ar && has_ma) {
    calc_block(P_mat, R_mat, w_t_vec, scale = phi)
  } else {
    matrix(0, 0, 0)
  }
  
  K_pb <- if (has_ar && has_xreg) {
    calc_block(P_mat, M_mat, w_t_vec, scale = phi)
  } else {
    matrix(0, 0, 0)
  }
  
  K_tt <- if (has_ma) {
    calc_block(R_mat, R_mat, w_t_vec, scale = phi)
  } else {
    matrix(0, 0, 0)
  }
  
  K_tb <- if (has_ma && has_xreg) {
    calc_block(R_mat, M_mat, w_t_vec, scale = phi)
  } else {
    matrix(0, 0, 0)
  }
  
  K_bb <- if (has_xreg) {
    calc_block(M_mat, M_mat, w_t_vec, scale = phi)
  } else {
    matrix(0, 0, 0)
  }
  
  # Blocks involving Phi ---
  K_aphi <- as.numeric(crossprod(s_vec, mu_eta_val * c_t_vec))
  
  K_pphi <- if (has_ar) {
    crossprod(P_mat, mu_eta_val * c_t_vec)
  } else {
    numeric(0)
  }
  
  K_tphi <- if (has_ma) {
    crossprod(R_mat, mu_eta_val * c_t_vec)
  } else {
    numeric(0)
  }
  
  K_bphi <- if (has_xreg) {
    crossprod(M_mat, mu_eta_val * c_t_vec)
  } else {
    numeric(0)
  }
  
  # Phi vs Phi Block ---
  K_phiphi <- sum(d_t_vec)
  
  # --------------------------------------------------------------------------
  # 6. ASSEMBLE AND RETURN THE FINAL FIM
  # --------------------------------------------------------------------------
  
  # Parameter order: Alpha, AR, MA, Beta, Phi
  n_params <- 1 + n_ar_params + n_ma_params + n_beta_params + 1
  fim <- matrix(NA_real_, nrow = n_params, ncol = n_params)
  
  # Define parameter indices
  idx_a <- 1
  
  idx_p <- if (has_ar) 2:(1 + n_ar_params) else integer(0)
  last_idx <- if (has_ar) max(idx_p) else idx_a
  
  idx_t <- if (has_ma) {
    (last_idx + 1):(last_idx + n_ma_params)
  } else {
    integer(0)
  }
  if (has_ma) last_idx <- max(idx_t)
  
  idx_b <- if (has_xreg) {
    (last_idx + 1):(last_idx + n_beta_params)
  } else {
    integer(0)
  }
  if (has_xreg) last_idx <- max(idx_b)
  
  idx_phi <- n_params
  
  # Fill Matrix (Symmetric) ---
  fim[idx_a, idx_a]   <- K_aa
  if (has_ar) {
    fim[idx_a, idx_p] <- K_ap
    fim[idx_p, idx_a] <- t(K_ap)
  }
  if (has_ma) {
    fim[idx_a, idx_t] <- K_at
    fim[idx_t, idx_a] <- t(K_at)
  }
  if (has_xreg) {
    fim[idx_a, idx_b] <- K_ab
    fim[idx_b, idx_a] <- t(K_ab)
  }
  fim[idx_a, idx_phi] <- K_aphi
  fim[idx_phi, idx_a] <- t(K_aphi)
  
  if (has_ar) {
    fim[idx_p, idx_p] <- K_pp
    if (has_ma) {
      fim[idx_p, idx_t] <- K_pt
      fim[idx_t, idx_p] <- t(K_pt)
    }
    if (has_xreg) {
      fim[idx_p, idx_b] <- K_pb
      fim[idx_b, idx_p] <- t(K_pb)
    }
    fim[idx_p, idx_phi] <- K_pphi
    fim[idx_phi, idx_p] <- t(K_pphi)
  }
  
  if (has_ma) {
    fim[idx_t, idx_t] <- K_tt
    if (has_xreg) {
      fim[idx_t, idx_b] <- K_tb
      fim[idx_b, idx_t] <- t(K_tb)
    }
    fim[idx_t, idx_phi] <- K_tphi
    fim[idx_phi, idx_t] <- t(K_tphi)
  }
  
  if (has_xreg) {
    fim[idx_b, idx_b] <- K_bb
    fim[idx_b, idx_phi] <- K_bphi
    fim[idx_phi, idx_b] <- t(K_bphi)
  }
  
  fim[idx_phi, idx_phi] <- K_phiphi
  
  # Name the matrix rows and columns ---
  names_varphi <- if (has_ar) paste0("varphi", ar_lags) else character(0)
  names_theta  <- if (has_ma) paste0("theta", ma_lags) else character(0)
  names_beta   <- if (has_xreg) colnames(xreg) else character(0)
  
  names_fim <- c("alpha", names_varphi, names_theta, names_beta, "phi")
  colnames(fim) <- names_fim
  rownames(fim) <- names_fim
  
  # --------------------------------------------------------------------------
  # 7. PREPARE FITTED VALUES AND OUTPUT
  # --------------------------------------------------------------------------
  
  # Create fitted values as a time series object
  fitted_ts <- ts(
    c(rep(NA, max_lag), mu_eff),
    start = start(y),
    frequency = frequency(y)
  )
  
  # Prepare output list with all required components
  output_list <- list(
    fisher_info_mat = fim,           # Fisher Information Matrix
    fitted_ts = fitted_ts,           # Fitted conditional means (ts object)
    muhat_effective = mu_eff,        # Effective fitted means (vector)
    etahat_full = eta,               # Full linear predictor (NA-padded)
    errorhat_full = error            # Full error/residuals (0-padded)
  )
  
  return(output_list)
}
