#' Build parameter index list for BARMA models
#'
#' Constructs the named list of positional indices mapping each parameter
#' group to its position in the flat parameter vector. Called by both
#' `barma()` and `run_optimizer()` to guarantee a single source of truth
#' for the parameter order.
#'
#' Parameter order:
#'   1. alpha   — intercept
#'   2. varphi  — AR coefficients (if present)
#'   3. theta   — MA coefficients (if present)
#'   4. beta    — regression coefficients (if present)
#'   5. phi     — precision parameter
#'
#' @param n_ar_params   Integer. Number of AR parameters.
#' @param n_ma_params   Integer. Number of MA parameters.
#' @param n_beta_params Integer. Number of regression parameters.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{`idx_alpha`}{Integer scalar. Index of the intercept.}
#'     \item{`idx_varphi`}{Integer vector. Indices of AR coefficients,
#'       or `integer(0)` if absent.}
#'     \item{`idx_theta`}{Integer vector. Indices of MA coefficients,
#'       or `integer(0)` if absent.}
#'     \item{`idx_beta`}{Integer vector. Indices of regression coefficients,
#'       or `integer(0)` if absent.}
#'     \item{`idx_phi`}{Integer scalar. Index of the precision parameter.}
#'     \item{`n_params`}{Integer scalar. Total number of parameters.}
#'   }
#' @keywords internal
#' @noRd
make_par_indices <- function(
    n_ar_params, 
    n_ma_params, 
    n_beta_params
) {
  
  idx_alpha  <- 1L
  last_idx   <- 1L
  
  idx_varphi <- seq_len(n_ar_params)   + last_idx
  last_idx   <- last_idx + n_ar_params
  
  idx_theta  <- seq_len(n_ma_params)   + last_idx
  last_idx   <- last_idx + n_ma_params
  
  idx_beta   <- seq_len(n_beta_params) + last_idx
  last_idx   <- last_idx + n_beta_params
  
  idx_phi    <- last_idx + 1L
  n_params   <- idx_phi
  
  list(
    idx_alpha  = idx_alpha,
    idx_varphi = idx_varphi,
    idx_theta  = idx_theta,
    idx_beta   = idx_beta,
    idx_phi    = idx_phi,
    n_params   = n_params
  )
}

# -------------------------------------------------------------------------- #
# OPTIMIZATION (CMLE / PCMLE) — internal helper
# Not exported. Called by barma().
# -------------------------------------------------------------------------- #

#' Run optimization for BARMA models (CMLE / PCMLE)
#'
#' Internal function that wraps `stats::optim` (BFGS or L-BFGS-B) and
#' `lbfgs::lbfgs` into a unified interface for standard Conditional Maximum 
#' Likelihood Estimation (CMLE) or Penalized CMLE (PCMLE) of BARMA/BARMAX models.
#'
#' @param y Numeric vector of observations on the unit interval (0, 1).
#' @param init_pars Numeric vector of initial parameter values. Order must
#'   match: `alpha`, `varphi` (AR), `theta` (MA), `beta` (regression), `phi`.
#' @param ar_lags Integer vector of AR lag indices, or `integer(0)` if absent.
#' @param ma_lags Integer vector of MA lag indices, or `integer(0)` if absent.
#' @param xreg Numeric matrix of exogenous regressors, or `NULL` if absent.
#' @param link Character string specifying the link function (e.g. `"logit"`).
#' @param par_idx Named list of parameter indices as returned by 
#'   \code{make_par_indices()}.
#' @param optimization A named list controlling the optimizer. Recognized
#'   fields:
#'   \describe{
#'     \item{`method`}{Character. One of `"BFGS"` (default), `"L-BFGS-B"`,
#'       or `"lbfgs"`.}
#'     \item{`lower`}{Numeric vector or scalar. Lower bounds for `"L-BFGS-B"`.
#'       Default `-Inf`.}
#'     \item{`upper`}{Numeric vector or scalar. Upper bounds for `"L-BFGS-B"`.
#'       Default `Inf`.}
#'   }
#' @param penalty Logical. If \code{TRUE}, applies ridge penalization (PCMLE). 
#'   If \code{FALSE}, performs standard CMLE.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{`opt`}{A list containing the estimated `par`, maximized `loglik`, 
#'       `counts` of function evaluations, the `convergence` code, and any 
#'       diagnostic `message`.}
#'     \item{`raw`}{The unmodified object returned by the underlying optimizer.}
#'   }
#'
#' @seealso [loglik_barma()], [score_vector_barma()], [barma()]
#' @keywords internal
#' @noRd
run_optimizer <- function(
    y,
    init_pars,
    ar_lags,
    ma_lags,
    xreg,
    link,
    par_idx,
    optimization,
    penalty
) {
  
  # 1. Setup =================================================================
  method <- optimization$method
  
  # Fallback to -Inf/Inf if the user omitted them from the list
  lower <- if (!is.null(optimization$lower)) optimization$lower else -Inf
  upper <- if (!is.null(optimization$upper)) optimization$upper else Inf
  
  has_xreg <- !is.null(xreg)
  
  n_ar_params   <- length(ar_lags)
  n_ma_params   <- length(ma_lags)
  n_beta_params <- if (has_xreg) ncol(xreg) else 0L
  
  # 2. Parameter indices =====================================================
  
  # Parameter order (matching init_pars order):
  #   1. alpha   — intercept
  #   2. varphi  — AR coefficients (if present)
  #   3. theta   — MA coefficients (if present)
  #   4. beta    — regression coefficients (if present)
  #   5. phi     — precision parameter
  
  idx_alpha  <- par_idx$idx_alpha
  idx_varphi <- par_idx$idx_varphi
  idx_theta  <- par_idx$idx_theta
  idx_beta   <- par_idx$idx_beta
  idx_phi    <- par_idx$idx_phi
  n_params   <- par_idx$n_params
  
  if (length(init_pars) != n_params)
    stop(
      "'init_pars' length (", length(init_pars), ") does not match ",
      "the expected number of parameters (", n_params, ")."
    )
  
  # 3. Helper closures =======================================================
  
  call_barma <- function(fn, x) {
    fn(
      y      = y,
      ar     = ar_lags,
      ma     = ma_lags,
      alpha  = x[idx_alpha],
      varphi = x[idx_varphi],
      theta  = x[idx_theta],
      phi    = x[idx_phi],
      link   = link,
      xreg   = xreg,
      beta   = x[idx_beta],
      penalty = penalty
    )
  }
  
  neg_loglik       <- function(x) -call_barma(loglik_barma,       x)
  neg_score_vector <- function(x) -call_barma(score_vector_barma, x)
  
  # 4. Execute Optimization ================================================
  
  opt_final <- tryCatch({
    
    if (method == "BFGS") {
      
      opt_aux <- stats::optim(
        par    = init_pars,
        fn     = neg_loglik,
        gr     = neg_score_vector,
        method = "BFGS"
      )
      
    } else if (method == "L-BFGS-B") {
      
      opt_aux <- stats::optim(
        par    = init_pars,
        fn     = neg_loglik,
        gr     = neg_score_vector,
        method = "L-BFGS-B",
        lower  = lower,
        upper  = upper
      )
      
    } else if (method == "lbfgs") {
      
      opt_aux <- lbfgs::lbfgs(
        vars      = init_pars,
        call_eval = neg_loglik,
        call_grad = neg_score_vector,
        invisible = 1
      )
      
    } else {
      
      stop("Unrecognized optimization method: '", method, 
           "'. Use 'BFGS', 'L-BFGS-B', or 'lbfgs'.")
      
    }
  
  # 5. Output ==============================================================
  
  beta_names <- if (has_xreg) {
    if (!is.null(colnames(xreg))) colnames(xreg) else 
      paste0("beta", seq_len(n_beta_params))
  } else {
    character(0)
  }
  par_names <- c(
    "alpha",
    if (n_ar_params   > 0) paste0("varphi", ar_lags),
    if (n_ma_params   > 0) paste0("theta",  ma_lags),
    beta_names,
    "phi"
  )
  
  counts <- if (method %in% c("BFGS", "L-BFGS-B")) {
    opt_aux$counts
  } else {
    NULL
  }
  
  opt_aux_par <- opt_aux$par
  names(opt_aux_par) <- par_names
  
  output <- list(
    
    opt = list(
      par         = opt_aux_par,
      loglik      = -opt_aux$value,
      convergence = opt_aux$convergence,
      message     = opt_aux$message,
      counts      = counts
    ),
    
    opt_raw = opt_aux
    
  )
  
  # Error message ==========================================================
  
  }, error = function(e) {
    stop("Optimization failed (method = '", method, "'): ", conditionMessage(e))
  })
  
  return(output)
  
}