#' Residuals for a Fitted betaARMA Model
#'
#' @description
#' Computes residuals for a fitted \eqn{\beta}ARMA model object of class
#' \code{"barma"}. Four types are available, controlled by the \code{type}
#' argument.
#'
#' @details
#' **Pearson residuals** (default) are defined on the response scale as:
#' \deqn{r_t = \frac{y_t - \hat{\mu}_t}{\sqrt{\hat{\mu}_t(1-\hat{\mu}_t)/(1+\hat{\phi})}}}
#' where \eqn{\hat{\mu}_t} is the fitted conditional mean and \eqn{\hat{\phi}}
#' is the estimated precision parameter. This residual was defined in Ferrari
#' and Cribari-Neto (2004) for beta regression models and is adopted here in
#' the \eqn{\beta}ARMA setting. These residuals are appropriate for ACF/PACF
#' analysis and portmanteau tests (Ljung-Box and Monti).
#'
#' **Quantile residuals** (Dunn and Smyth, 1996) are defined as:
#' \deqn{r_t^Q = \Phi^{-1}\!\left(F(y_t;\,\hat{\mu}_t\hat{\phi},\,(1-\hat{\mu}_t)\hat{\phi})\right)}
#' where \eqn{F(\cdot\,;\,a,b)} is the Beta CDF with shape parameters
#' \eqn{a = \hat{\mu}_t\hat{\phi}} and \eqn{b = (1-\hat{\mu}_t)\hat{\phi}},
#' and \eqn{\Phi^{-1}} is the standard normal quantile function. Under a
#' correctly specified model, quantile residuals are approximately i.i.d.
#' \eqn{N(0,1)}, making them suitable for normality assessment via
#' histograms and Q-Q plots.
#' 
#' **Link-scale residuals** are defined on the predictor scale as:
#' \deqn{r_t^L = \frac{g(y_t) - \hat{\eta}_t}{\sqrt{[g'(\hat{\mu}_t)]^2
#'   \cdot \hat{\mu}_t(1-\hat{\mu}_t)/(1+\hat{\phi})}}}
#' where \eqn{g(\cdot)} is the link function, \eqn{\hat{\eta}_t} is the
#' fitted linear predictor, and \eqn{g'(\hat{\mu}_t) = d\eta/d\mu}. This
#' residual was used in earlier versions of the package and is retained for
#' compatibility.
#' 
#' **Raw residuals** are the simple difference on the response scale:
#' \deqn{r_t^{\mathrm{raw}} = y_t - \hat{\mu}_t}
#'
#'
#' @param object
#'   A fitted model object of class \code{"barma"}, as returned by
#'   \code{\link{barma}}.
#' @param type
#'   Character string specifying the type of residuals to compute. One of:
#'   \describe{
#'     \item{\code{"pearson"} (default)}{Pearson residuals on the response
#'       scale, as defined in Ferrari and Cribari-Neto (2004). Recommended
#'       for ACF/PACF analysis and portmanteau tests.}
#'     \item{\code{"quantile"}}{Quantile residuals (Dunn and Smyth, 1996),
#'       approximately i.i.d. \eqn{N(0,1)} under a correctly specified model.
#'       Recommended for normality assessment.}
#'     \item{\code{"link"}}{Residuals on the predictor (link) scale. Retained
#'       for compatibility with earlier versions of the package.}
#'     \item{\code{"raw"}}{Raw residuals \eqn{y_t - \hat{\mu}_t} on the
#'       response scale.}
#'   }
#' @param ...
#'   Additional arguments (currently unused, included for S3 consistency).
#'
#' @return
#'   A numeric \code{ts} object of the same length as the response series,
#'   with the first \code{max_lag} values set to \code{NA}. The \code{ts}
#'   attributes (\code{start} and \code{frequency}) are inherited from the
#'   original response series stored in \code{object$y}.
#'
#' @references
#'   Dunn, P. K., & Smyth, G. K. (1996). Randomized quantile residuals.
#'   \emph{Journal of Computational and Graphical Statistics}, 5(3), 236--244.
#'   \doi{10.1080/10618600.1996.10474708}
#'
#'   Ferrari, S. L. P., & Cribari-Neto, F. (2004). Beta regression for
#'   modelling rates and proportions. \emph{Journal of Applied Statistics},
#'   31(7), 799--815. \doi{10.1080/0266476042000214501}
#'
#'   Rocha, A. V., & Cribari-Neto, F. (2009). Beta autoregressive moving
#'   average models. \emph{TEST}, 18(3), 529--545.
#'   \doi{10.1007/s11749-008-0112-z}
#'
#'   Rocha, A. V., & Cribari-Neto, F. (2017). Erratum to: Beta autoregressive
#'   moving average models. \emph{TEST}, 26, 451--459.
#'   \doi{10.1007/s11749-017-0528-4}
#'
#' @seealso
#'   \code{\link{barma}} for model fitting;
#'   \code{\link{fitted.barma}} for fitted values;
#'   \code{\link{plot.barma}} for diagnostic plots.
#'
#' @examples
#' data(brasilia_ts)
#' fit <- barma(brasilia_ts, ar = 1, ma = 1)
#'
#' # Pearson residuals (default) -- Ferrari & Cribari-Neto (2004)
#' res <- residuals(fit)
#'
#' # Raw residuals
#' res_raw <- residuals(fit, type = "raw")
#'
#' # Quantile residuals -- recommended for normality assessment
#' res_q <- residuals(fit, type = "quantile")
#'
#' # Link-scale residuals
#' res_link <- residuals(fit, type = "link")
#'
#' @importFrom stats qnorm pbeta
#' @export
#' @method residuals barma
residuals.barma <- function(object,
                            type = c("pearson", "quantile", "link", "raw"),
                            ...) {
  
  type <- match.arg(type)
  
  # ---------------------------------------------------------------------------
  # 1. Extract components from the barma object
  # ---------------------------------------------------------------------------
  y <- object$y
  if (is.null(y)) {
    stop("Original time series 'y' not found in the 'barma' object.")
  }
  
  fitted_mus <- object$fitted  # Fitted means mu_hat (full length, NA-padded ts)
  phi        <- object$phi     # Estimated precision parameter
  max_lag    <- object$max_lag # Number of initial NAs (= max(ar, ma) lags)
  
  if (is.null(max_lag)) {
    stop("'max_lag' not found in the 'barma' object.")
  }
  
  n_obs <- length(y)
  
  if (n_obs <= max_lag) {
    warning("Insufficient observations relative to max_lag.")
    return(rep(NA_real_, n_obs))
  }
  
  # ---------------------------------------------------------------------------
  # 2. Restrict to effective observations (drop the first max_lag NAs)
  # ---------------------------------------------------------------------------
  idx_effective   <- (max_lag + 1):n_obs
  y_effective     <- as.numeric(y[idx_effective])
  muhat_effective <- as.numeric(fitted_mus[idx_effective])
  
  # ---------------------------------------------------------------------------
  # 3. Compute residuals by type
  # ---------------------------------------------------------------------------
  if (type == "pearson") {
    
    # Pearson residuals on the response scale.
    # Reference: Ferrari & Cribari-Neto (2004), adapted to the betaARMA setting.
    # r_t = (y_t - mu_hat_t) / sqrt( mu_hat_t*(1 - mu_hat_t) / (1 + phi) )
    Vmu_hat <- muhat_effective * (1 - muhat_effective)
    resids  <- (y_effective - muhat_effective) / sqrt(Vmu_hat / (1 + phi))
    
  } else if (type == "raw") {
    
    # Raw residuals on the response scale: y_t - mu_hat_t
    resids <- y_effective - muhat_effective
    
  } else if (type == "quantile") {
    
    # Quantile residuals: Phi^{-1}( F(y_t; mu_hat_t*phi, (1-mu_hat_t)*phi) )
    # Reference: Dunn & Smyth (1996).
    # Under a correctly specified model, these are approximately i.i.d. N(0,1).
    shape1 <- muhat_effective * phi
    shape2 <- (1 - muhat_effective) * phi
    resids  <- stats::qnorm(
      stats::pbeta(y_effective, shape1 = shape1, shape2 = shape2)
    )
    
  } else {
    
    # Link-scale residuals: (g(y_t) - eta_hat_t) / sd(r_t) on predictor scale.
    # Retained for compatibility with earlier versions of the package.
    etahat <- object$etahat  # Fitted linear predictor (full length, NA-padded)
    link   <- object$link
    
    link_structure <- make_link_structure(link)
    linkfun    <- link_structure$linkfun   # g(mu): maps (0,1) -> R
    mu.eta_fun <- link_structure$mu.eta    # dmu/deta
    
    etahat_effective <- etahat[idx_effective]
    
    # Numerator: g(y_t) - eta_hat_t
    num <- linkfun(y_effective) - etahat_effective
    
    # Denominator: sqrt( [g'(mu_hat_t)]^2 * mu_hat_t*(1-mu_hat_t) / (1+phi) )
    # where g'(mu) = deta/dmu = 1 / (dmu/deta)
    dmu_deta_sq <- mu.eta_fun(eta = etahat_effective)^2
    Vmu_hat     <- muhat_effective * (1 - muhat_effective)
    den         <- sqrt((1 / dmu_deta_sq) * Vmu_hat / (1 + phi))
    
    resids <- num / den
    
  }
  
  # ---------------------------------------------------------------------------
  # 4. Pad the first max_lag positions with NA and return as ts object
  # ---------------------------------------------------------------------------
  final_resids <- c(rep(NA_real_, max_lag), resids)
  
  if (stats::is.ts(y)) {
    final_resids <- stats::ts(
      final_resids,
      start     = stats::start(y),
      frequency = stats::frequency(y)
    )
  }
  
  return(final_resids)
}
