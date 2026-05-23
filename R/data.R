#' Monthly relative humidity in Brasília
#'
#' @description
#' Monthly relative humidity in Brasília, Brazil, expressed as a proportion.
#' The data were obtained from the NASA POWER project and cover
#' January 1999 to June 2024.
#'
#' @format A \code{ts} object with monthly observations, frequency 12,
#'   starting in January 1999.
#'
#' @source NASA Prediction of Worldwide Energy Resources (POWER) project.
#'   \url{https://power.larc.nasa.gov/}
#'
#' @references
#' Cribari-Neto, F., Costa, E., & Fonseca, R. V. (2025). Numerical stability
#' enhancements in beta autoregressive moving average model estimation.
#' \emph{Brazilian Journal of Probability and Statistics}, 39(4), 410--437.
#' \doi{10.1214/25-BJPS645}
#'
#' @seealso \code{\link{brasilia_df}} for the same data as a data frame.
#'
#' @examples
#' data(brasilia_ts)
#' plot(brasilia_ts, ylab = "Relative humidity (proportion)", xlab = "Year")
"brasilia_ts"


#' Monthly relative humidity in Brasília (data frame)
#'
#' @description
#' The same series as \code{\link{brasilia_ts}}, stored as a data frame
#' with a \code{yearmon} time column.
#'
#' @format A data frame with 2 columns:
#' \describe{
#'   \item{time}{Month of observation (\code{yearmon}).}
#'   \item{y}{Relative humidity as a proportion (numeric, in (0, 1)).}
#' }
#'
#' @source NASA Prediction of Worldwide Energy Resources (POWER) project.
#'   \url{https://power.larc.nasa.gov/}
#'
#' @seealso \code{\link{brasilia_ts}} for the same data as a \code{ts} object.
#'
#' @examples
#' data(brasilia_df)
#' head(brasilia_df)
"brasilia_df"