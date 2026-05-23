#' Diagnostic Plots for a Fitted betaARMA Model
#'
#' @description
#' Produces diagnostic plots for a \eqn{\beta}ARMA model fitted by
#' \code{\link{barma}}. By default (\code{which = "default"}), four panels are
#' displayed in a 2\eqn{\times}2 grid: observed vs. fitted values, residual
#' ACF, Ljung-Box p-values, and residual PACF.
#'
#' @param x
#'   An object of class \code{"barma"}, as returned by \code{\link{barma}}.
#' @param which
#'   A character string controlling which panel(s) to display. One of:
#'   \describe{
#'     \item{\code{"default"} (default)}{Four-panel 2\eqn{\times}2 grid:
#'       observed vs. fitted (top left), residual ACF (top right),
#'       Ljung-Box p-values (bottom left), and residual PACF (bottom right).}
#'     \item{\code{"all"}}{Six-panel 3\eqn{\times}2 grid:
#'       observed vs. fitted (top left), residuals over time (top right),
#'       residual ACF (middle left), residual PACF (middle right),
#'       Ljung-Box p-values (bottom left), Monti p-values (bottom right).}
#'     \item{\code{"fitted"}}{Single panel: observed vs. fitted values.}
#'     \item{\code{"tsplot"}}{Single panel: residuals plotted as a time series,
#'       with a horizontal reference line at zero and dashed lines at
#'       \eqn{\pm 3} (ad hoc threshold; not shown for \code{"raw"} residuals).}
#'     \item{\code{"acf"}}{Single panel: residual ACF.}
#'     \item{\code{"pacf"}}{Single panel: residual PACF.}
#'     \item{\code{"ljungbox"}}{Single panel: Ljung-Box p-values.}
#'     \item{\code{"monti"}}{Single panel: Monti test p-values.}
#'     \item{\code{"hist"}}{Single panel: residual distribution histogram with
#'       kernel density overlay.}
#'     \item{\code{"qq"}}{Single panel: Normal Q-Q plot. Always uses quantile
#'       residuals regardless of \code{residual_type}; a message is issued if
#'       the type is overridden.}
#'   }
#' @param residual_type
#'   Character string passed to \code{\link{residuals.barma}} controlling
#'   which residuals are computed and displayed. One of \code{"pearson"}
#'   (default), \code{"quantile"}, \code{"link"} or \code{"raw"}. Pearson
#'   residuals are recommended for ACF/PACF analysis and portmanteau
#'   tests (Scher, Cribari-Neto, Pumi, and Bayer, 2020); quantile residuals
#'   are recommended for normality assessment (Dunn and Smyth, 1996). See also
#'   Ferrari and Cribari-Neto (2004).
#' @param lag_max
#'   Maximum lag to display in the ACF, PACF, Ljung-Box, and Monti panels.
#'   Defaults to \code{24}, appropriate for monthly series.
#' @param title
#'   An optional character string used as the overall title of the plot.
#'   Defaults to \code{NULL} (no title).
#' @param colour_observed
#'   Colour for observed values in the fitted-values panel.
#'   Defaults to \code{"#00BFC4"} (teal).
#' @param colour_fitted
#'   Colour for fitted values in the fitted-values panel.
#'   Defaults to \code{"#F8766D"} (red).
#' @param colour_residual
#'   Colour for residual lines, ACF/PACF bars, portmanteau test points, and
#'   histogram fill. Defaults to \code{"#7CAE00"} (green).
#' @param ...
#'   Additional arguments (currently unused, included for S3 consistency).
#'
#' @details
#' **Residual type.** All panels (except \emph{Observed vs. Fitted}) use the
#' residuals selected by \code{residual_type}. The residual type is displayed
#' in each panel subtitle for transparency.
#'
#' **Reference lines in the residuals-over-time panel.** The dashed horizontal 
#' lines are drawn at \eqn{\pm 3} as an ad hoc threshold for identifying 
#' atypical observations. No distributional assumption is made; observations 
#' outside these bounds may warrant individual inspection. Pearson and 
#' link-scale residuals do not follow a standard normal distribution, so 
#' normal-distribution-based quantiles (e.g., \eqn{\pm 1.96} or 
#' \eqn{\pm 2.576}) are not appropriate reference values for such residuals. 
#' Quantile residuals are approximately distributed as \eqn{N(0,1)} under a 
#' correctly specified model, so \eqn{\pm 3} is a conservative but reasonable 
#' threshold for them.
#'
#' **Q-Q plot residuals.** The \code{"qq"} panel always uses quantile
#' residuals (Dunn and Smyth, 1996), regardless of \code{residual_type}.
#' Quantile residuals are the only type theoretically expected to follow
#' \eqn{N(0,1)} under a correctly specified model; using Pearson, raw, or
#' link-scale residuals in a normal Q-Q plot would produce misleading results.
#' A message is issued when \code{residual_type} is overridden.
#'
#' **Effective sample size.** Both the Ljung-Box and Monti test statistics
#' rely on the effective sample size \eqn{N_{eff} = n - \max(p, q)} (where
#' \eqn{\max(p, q)} corresponds to the maximum lag in the model).
#'
#' **Ljung-Box test.** The p-values are computed from the Ljung-Box test
#' statistic.
#' \eqn{Q_{LB}(k) = N_{eff}(N_{eff}+2)\sum_{j=1}^{k}\hat{\rho}_j^2/(N_{eff}-j)},
#' where \eqn{\hat{\rho}_j} is the \eqn{j}-th sample autocorrelation of the
#' residuals and \eqn{N_{eff} = n - \max(p, q)} is the effective sample size.
#' The degrees of freedom are adjusted by subtracting the total
#' number of estimated  AR and MA parameters,
#' \eqn{n_{ar} + n_{ma}} (Scher, Cribari-Neto, Pumi, and Bayer (2020)). The
#'  first \eqn{n_{ar} + n_{ma}} lags are set to \code{NA}
#' since the test statistic is not defined for \eqn{k \leq n_{ar} + n_{ma}}.
#'
#' **Monti test.** The p-values are computed from the test statistic
#' \eqn{M(k) = N_{eff}(N_{eff}+2)\sum_{j=1}^{k}\hat{\rho}_{jj}^2/(N_{eff}-j)},
#' where \eqn{\hat{\rho}_{jj}} is the \eqn{j}-th sample partial autocorrelation
#' of the residuals, and \eqn{N_{eff} = n - \max\_lag} is the effective sample
#' size after stripping initial NAs. The degrees of freedom are adjusted to
#' \eqn{k - (n_{ar} + n_{ma})}
#' (Monti, 1994; Scher, Cribari-Neto, Pumi, and Bayer (2020)). The first
#' \eqn{n_{ar} + n_{ma}} lags are set to \code{NA} as the statistic is not
#' defined for \eqn{k \leq n_{ar} + n_{ma}}.
#'
#' @return
#'   For grid requests (\code{"default"} or \code{"all"}), invisibly returns a
#'   \code{gtable} object produced by \code{gridExtra::arrangeGrob}, which can
#'   be saved with \code{ggplot2::ggsave}. For single-panel requests, invisibly
#'   returns the individual \code{ggplot} object.
#'
#' @seealso
#'   \code{\link{barma}} for model fitting;
#'   \code{\link{residuals.barma}} for the residual extraction method;
#'   \code{\link{fitted.barma}} for the fitted-value extraction method;
#'   \code{\link{forecast.barma}} for out-of-sample forecasting.
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
#'   Monti, A. C. (1994). A proposal for a residual autocorrelation test in
#'   linear models. \emph{Biometrika}, 81(4), 776--780.
#'   \doi{10.1093/biomet/81.4.776}
#'
#'   Rocha, A. V., & Cribari-Neto, F. (2009). Beta autoregressive moving
#'   average models. \emph{TEST}, 18(3), 529--545.
#'   \doi{10.1007/s11749-008-0112-z}
#'
#'   Rocha, A. V., & Cribari-Neto, F. (2017). Erratum to: Beta autoregressive
#'   moving average models. \emph{TEST}, 26, 451--459.
#'   \doi{10.1007/s11749-017-0528-4}
#'
#' @examples
#' data(brasilia_ts)
#' fit <- barma(brasilia_ts, ar = 1, ma = 1)
#'
#' # Default four-panel diagnostic grid (Pearson residuals)
#' plot(fit)
#'
#' # Six-panel diagnostic grid
#' plot(fit, which = "all")
#'
#' # Single panel: residuals as time series
#' plot(fit, which = "tsplot")
#'
#' # Single panel: Ljung-Box p-values only
#' plot(fit, which = "ljungbox")
#'
#' # Single panel: Monti p-values only
#' plot(fit, which = "monti")
#'
#' # Normal Q-Q plot (always uses quantile residuals)
#' plot(fit, which = "qq")
#'
#' # Histogram with quantile residuals
#' plot(fit, which = "hist", residual_type = "quantile")
#'
#' @importFrom forecast ggAcf ggPacf
#' @importFrom gridExtra arrangeGrob
#' @importFrom grid grid.draw
#' @importFrom ggplot2 ggplot aes geom_line geom_hline geom_point
#'   geom_histogram geom_density stat_qq stat_qq_line scale_colour_manual
#'   scale_y_continuous scale_x_continuous labs theme_minimal theme
#'   element_text element_blank after_stat
#' @importFrom stats acf pacf time fitted residuals density pchisq
#' @importFrom rlang .data
#' @export
plot.barma <- function(x,
                       which = c(
                         "default", "all", "fitted",
                         "tsplot", "acf", "pacf",
                         "ljungbox", "monti", "hist", "qq"
                       ),
                       residual_type = c("pearson", "quantile", "link", "raw"),
                       lag_max = 24,
                       title = NULL,
                       colour_observed = "#00BFC4",
                       colour_fitted = "#F8766D",
                       colour_residual = "#7CAE00",
                       ...) {
  # --------------------------------------------------------------------------
  # 0. Input validation
  # --------------------------------------------------------------------------
  if (!inherits(x, "barma")) {
    stop("'x' must be an object of class \"barma\".")
  }

  which <- match.arg(which)
  residual_type <- match.arg(residual_type)

  # Q-Q plot requires quantile residuals — override silently with a message
  if (which == "qq" && residual_type != "quantile") {
    message(
      "Note: the Q-Q plot uses quantile residuals regardless of ",
      "'residual_type'. Switching to residual_type = \"quantile\"."
    )
    residual_type <- "quantile"
  }

  # Label used in panel subtitles to make residual type explicit
  resid_label <- switch(residual_type,
    pearson  = "Pearson residuals",
    raw      = "Raw residuals",
    quantile = "Quantile residuals",
    link     = "Link-scale residuals"
  )

  # --------------------------------------------------------------------------
  # 1. Extract components from the barma object
  # --------------------------------------------------------------------------
  y_obs <- as.numeric(x$y)
  y_fit <- as.numeric(fitted(x))
  y_res <- as.numeric(residuals(x, type = residual_type))
  t_axis <- as.numeric(time(x$y))
  n <- length(y_obs)

  # Pad fitted values with leading NAs
  n_pad_fit <- n - length(y_fit)
  if (n_pad_fit > 0) {
    y_fit <- c(rep(NA_real_, n_pad_fit), y_fit)
  }

  # Guard against edge cases where residual lengths differ
  n_pad_res <- n - length(y_res)
  if (n_pad_res > 0) {
    y_res_padded <- c(rep(NA_real_, n_pad_res), y_res)
  } else {
    y_res_padded <- y_res
  }

  # Strip NAs for panels that require a clean residual vector
  y_res_clean <- y_res_padded[!is.na(y_res_padded)]
  n_clean <- length(y_res_clean)

  # Degrees of freedom for portmanteau tests: Number of estimated parameters
  n_ar <- length(x$ar_lags)
  n_ma <- length(x$ma_lags)
  fitdf <- n_ar + n_ma

  lags <- seq_len(lag_max)

  # --------------------------------------------------------------------------
  # 2. Compute Ljung-Box p-values (vectorized)
  #
  # Original Ljung-Box formula:
  # Q(k) = N_eff * (N_eff + 2) * sum_{j=1}^k r_j^2 / (N_eff - j)
  #
  # Variable Mapping in this code:
  # N_eff = n_clean (effective sample size after stripping initial NAs).
  # j     = lags    (the sequence 1, 2, ..., lag_max).
  #
  # Degrees of freedom for pchisq: k - fitdf
  # (Scher, Cribari-Neto, Pumi, and Bayer (2020)).
  # The first 'fitdf' lags are set to NA because the chi-square distribution
  # requires strictly positive degrees of freedom (k > fitdf).
  # --------------------------------------------------------------------------
  valid_lags <- lags > fitdf

  # acf()$acf includes lag 0 (always 1.0); drop it to get lags 1..lag_max
  acf_vals <- as.numeric(
    stats::acf(y_res_clean, lag.max = lag_max, plot = FALSE)$acf
  )[-1]

  lb_inner_terms <- (acf_vals^2) / (n_clean - lags)
  Q_stats <- n_clean * (n_clean + 2) * cumsum(lb_inner_terms)

  lb_pvals <- rep(NA_real_, lag_max)
  lb_pvals[valid_lags] <- stats::pchisq(
    Q_stats[valid_lags],
    df         = lags[valid_lags] - fitdf,
    lower.tail = FALSE
  )
  lb_df <- data.frame(lag = lags, pvalue = lb_pvals)

  # --------------------------------------------------------------------------
  # 3. Compute Monti test p-values (vectorized)
  #
  # Original Monti (1994) formula:
  # M(k) = N_eff * (N_eff + 2) * sum_{j=1}^k rho_jj^2 / (N_eff - j)
  #
  # Variable Mapping in this code:
  # N_eff = n_clean (effective sample size after stripping initial NAs).
  # j     = lags    (the sequence 1, 2, ..., lag_max).
  #
  # Degrees of freedom for pchisq: k - (n_ar + n_ma)
  # (Scher, Cribari-Neto, Pumi, and Bayer (2020)).
  # --------------------------------------------------------------------------
  pacf_vals <- as.numeric(
    stats::pacf(y_res_clean, lag.max = lag_max, plot = FALSE)$acf
  )

  inner_terms <- (pacf_vals^2) / (n_clean - lags)
  M_stats <- n_clean * (n_clean + 2) * cumsum(inner_terms)

  monti_pvals <- rep(NA_real_, lag_max)
  valid_lags <- lags > fitdf
  monti_pvals[valid_lags] <- stats::pchisq(
    M_stats[valid_lags],
    df         = lags[valid_lags] - fitdf,
    lower.tail = FALSE
  )
  monti_df <- data.frame(lag = lags, pvalue = monti_pvals)

  # --------------------------------------------------------------------------
  # 4. Shared data frame and ggplot theme
  # --------------------------------------------------------------------------
  df <- data.frame(
    time     = t_axis,
    observed = y_obs,
    fitted   = y_fit,
    residual = y_res_padded
  )

  barma_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text       = ggplot2::element_text(size = 9),
      axis.title      = ggplot2::element_text(size = 9),
      legend.text     = ggplot2::element_text(size = 9),
      legend.title    = ggplot2::element_blank(),
      legend.position = "bottom"
    )

  # Shared x-axis breaks aligned to seasonal lags (6, 12, 18, 24, ...)
  x_breaks <- ggplot2::scale_x_continuous(breaks = seq(6, lag_max, by = 6))

  # --------------------------------------------------------------------------
  # 5. Build all individual panels
  # --------------------------------------------------------------------------

  # --- Panel 1: Observed vs. Fitted ---
  p_fitted <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time)) +
    ggplot2::geom_line(ggplot2::aes(y = .data$observed, colour = "Observed"),
      linewidth = 0.6, na.rm = TRUE
    ) +
    ggplot2::geom_line(ggplot2::aes(y = .data$fitted, colour = "Fitted"),
      linewidth = 0.6, na.rm = TRUE
    ) +
    ggplot2::scale_colour_manual(
      values = c("Observed" = colour_observed, "Fitted" = colour_fitted),
      breaks = c("Observed", "Fitted")
    ) +
    ggplot2::labs(x = "Time", y = "Value", title = "Observed vs. Fitted") +
    barma_theme

  # --- Panel 2: Residuals over time ---
  # Reference lines at ±3 as an ad hoc threshold. No distributional
  # assumption is made — see @details for the rationale.

  # 2.1 Dynamically set the subtitle based on residual type
  resid_subtitle <- if (residual_type == "raw") {
    resid_label
  } else {
    paste0(resid_label, "  |  dashed lines: \u00b13 (ad hoc)")
  }

  # 2.2. Build the base plot (without the +/- 3 lines)
  p_resid <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$time,
      y = .data$residual
    )
  ) +
    ggplot2::geom_line(
      colour = colour_residual, linewidth = 0.6,
      na.rm = TRUE
    ) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "solid",
      colour = "grey20", linewidth = 0.3
    ) +
    ggplot2::labs(
      x        = "Time",
      y        = resid_label,
      title    = "Residuals over Time",
      subtitle = resid_subtitle
    ) +
    barma_theme

  # 2.3. Conditionally add the dashed reference lines
  if (residual_type != "raw") {
    p_resid <- p_resid +
      ggplot2::geom_hline(
        yintercept = c(-3, 3), linetype = "dashed",
        colour = "grey40", linewidth = 0.4
      )
  }

  # --- Panel 3: Residual ACF ---
  # suppressMessages() silences "Scale for x is already present" warning
  # that arises when overriding the default x scale set by ggAcf internally.
  p_acf <- suppressMessages(
    forecast::ggAcf(y_res_clean, lag.max = lag_max) +
      x_breaks +
      ggplot2::labs(
        title = "Residual ACF",
        subtitle = resid_label,
        x = "Lag", y = "ACF"
      ) +
      barma_theme
  )

  # --- Panel 4: Residual PACF ---
  # suppressMessages() silences "Scale for x is already present" warning
  # that arises when overriding the default x scale set by ggPacf internally.
  p_pacf <- suppressMessages(
    forecast::ggPacf(y_res_clean, lag.max = lag_max) +
      x_breaks +
      ggplot2::labs(
        title = "Residual PACF",
        subtitle = resid_label,
        x = "Lag", y = "PACF"
      ) +
      barma_theme
  )

  # --- Panel 5: Ljung-Box p-values ---
  p_lb <- ggplot2::ggplot(
    lb_df,
    ggplot2::aes(x = .data$lag, y = .data$pvalue)
  ) +
    ggplot2::geom_point(colour = colour_residual, size = 2, na.rm = TRUE) +
    ggplot2::geom_hline(
      yintercept = 0.05, linetype = "dashed",
      colour = "grey40", linewidth = 0.4
    ) +
    x_breaks +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      x = "Lag", y = "p-value",
      title = "Ljung-Box p-values",
      subtitle = resid_label
    ) +
    barma_theme

  # --- Panel 6: Monti test p-values ---
  p_monti <- ggplot2::ggplot(
    monti_df,
    ggplot2::aes(x = .data$lag, y = .data$pvalue)
  ) +
    ggplot2::geom_point(colour = colour_residual, size = 2, na.rm = TRUE) +
    ggplot2::geom_hline(
      yintercept = 0.05, linetype = "dashed",
      colour = "grey40", linewidth = 0.4
    ) +
    x_breaks +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      x = "Lag", y = "p-value",
      title = "Monti p-values",
      subtitle = resid_label
    ) +
    barma_theme

  # --- Panel 7: Residual distribution histogram (built on demand) ---
  needs_hist <- which %in% c("hist")
  p_hist <- if (needs_hist) {
    ggplot2::ggplot(
      data.frame(residual = y_res_clean),
      ggplot2::aes(x = .data$residual)
    ) +
      ggplot2::geom_histogram(
        ggplot2::aes(y = ggplot2::after_stat(density)),
        bins = 20,
        fill = colour_residual,
        colour = "white",
        alpha = 0.8
      ) +
      ggplot2::geom_density(linewidth = 0.7, colour = "grey20") +
      ggplot2::labs(
        x = resid_label,
        y = "Density",
        title = "Residual Distribution",
        subtitle = resid_label
      ) +
      barma_theme
  } else {
    NULL
  }

  # --- Panel 8: Normal Q-Q plot (built on demand; always quantile residuals) ---
  # Quantile residuals are the only type with a theoretical N(0,1) guarantee
  # under a correctly specified model (Dunn & Smyth, 1996). Residuals for the
  # Q-Q plot are therefore always computed with type = "quantile", independent
  # of the global residual_type chosen by the user.
  needs_qq <- which %in% c("qq", "all")
  p_qq <- if (needs_qq) {
    qq_res <- as.numeric(residuals(x, type = "quantile"))
    qq_res <- qq_res[!is.na(qq_res)]
    n_qq <- length(qq_res)

    conf_level <- 0.95
    p_vals <- (seq_len(n_qq) - 0.5) / n_qq
    q_theor <- stats::qnorm(p_vals)
    se <- sqrt(p_vals * (1 - p_vals) / n_qq) / stats::dnorm(q_theor)
    z_val <- stats::qnorm((1 + conf_level) / 2)
    mean_qq <- mean(qq_res, na.rm = TRUE)
    sd_qq <- stats::sd(qq_res, na.rm = TRUE)

    qq_band_df <- data.frame(
      theoretical = q_theor,
      ymin = (q_theor - z_val * se) * sd_qq + mean_qq,
      ymax = (q_theor + z_val * se) * sd_qq + mean_qq
    )

    ggplot2::ggplot() +
      ggplot2::geom_ribbon(
        data = qq_band_df,
        ggplot2::aes(
          x = .data$theoretical,
          ymin = .data$ymin,
          ymax = .data$ymax
        ),
        fill = colour_residual, alpha = 0.2
      ) +
      ggplot2::stat_qq(
        data = data.frame(residual = qq_res),
        ggplot2::aes(sample = .data$residual),
        colour = colour_residual, size = 2, alpha = 0.6
      ) +
      ggplot2::stat_qq_line(
        data = data.frame(residual = qq_res),
        ggplot2::aes(sample = .data$residual),
        colour = "grey20", linewidth = 0.7
      ) +
      ggplot2::labs(
        x        = "Theoretical Quantiles",
        y        = "Sample Quantiles",
        title    = "Normal Q-Q Plot",
        subtitle = "Quantile residuals"
      ) +
      barma_theme
  } else {
    NULL
  }

  # --------------------------------------------------------------------------
  # 6. Return single panel if requested
  # --------------------------------------------------------------------------
  if (!(which %in% c("all", "default"))) {
    p_single <- switch(which,
      fitted   = p_fitted,
      tsplot   = p_resid,
      acf      = p_acf,
      pacf     = p_pacf,
      ljungbox = p_lb,
      monti    = p_monti,
      hist     = p_hist,
      qq       = p_qq
    )
    print(p_single)
    return(invisible(p_single))
  }

  # --------------------------------------------------------------------------
  # 7. Return "default" (4-panel grid):
  #    Observed vs. Fitted  |  Ljung-Box p-values
  #    Residual ACF         |  Residual PACF
  # --------------------------------------------------------------------------
  if (which == "default") {
    g <- gridExtra::arrangeGrob(
      grobs = list(
        p_fitted, p_lb,
        p_acf, p_pacf
      ),
      ncol = 2,
      top = title
    )
    grid::grid.draw(g)
    return(invisible(g))
  }

  # --------------------------------------------------------------------------
  # 8. Return "all" (6-panel grid):
  #    Observed vs. Fitted  |  Residuals over Time
  #    Residual ACF         |  Residual PACF
  #    Ljung-Box p-values   |  Monti p-values
  # --------------------------------------------------------------------------
  if (which == "all") {
    g <- gridExtra::arrangeGrob(
      grobs = list(
        p_fitted, p_resid,
        p_acf,    p_pacf,
        p_lb,     p_monti
      ),
      ncol = 2,
      top = title
    )
    grid::grid.draw(g)
    return(invisible(g))
  }
}