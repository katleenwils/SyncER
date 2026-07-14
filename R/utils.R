#' Check Distribution Shape via Skewness and Kurtosis Z-Tests
#'
#' Flags a numeric vector as approximately normal-shaped using the standard
#' errors of skewness and kurtosis (Tabachnick & Fidell rule of thumb):
#'   z_skewness = skewness / sqrt(6 / N)
#'   z_kurtosis = excess_kurtosis / sqrt(24 / N)
#' If either |z| exceeds the critical value, the distribution is flagged
#' non-normal on that characteristic. Unlike Shapiro-Wilk, this still uses N
#' explicitly in its standard error, but is the conventional approach for
#' assessing skew/peakedness rather than overall distributional fit.
#'
#' @param x Numeric vector to test.
#' @param critical_z Critical z value (default 1.96, the .05 significance level;
#'   use 2.58 for the .01 level).
#'
#' @return Logical: TRUE if both z_skewness and z_kurtosis are within bounds.
#'
#' @keywords internal
shape_ok <- function(x, critical_z = 1.96) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 3) return(NA)

  m <- mean(x)
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(NA)

  skewness <- (sum((x - m)^3) / n) / s^3
  excess_kurtosis <- (sum((x - m)^4) / n) / s^4 - 3

  z_skewness <- skewness / sqrt(6 / n)
  z_kurtosis <- excess_kurtosis / sqrt(24 / n)

  abs(z_skewness) <= critical_z && abs(z_kurtosis) <= critical_z
}

#' Current BP Datum (Year - 1950)
#'
#' Returns the number of years elapsed since the radiocarbon BP datum (1950),
#' calculated from the current system date. Used as the default \code{offset}/
#' \code{age_offset} argument throughout \emph{SyncER} instead of a hardcoded
#' value, so the default stays correct in future years.
#'
#' @return Integer: current year minus 1950.
#'
#' @export
bp_datum <- function() as.integer(format(Sys.Date(), "%Y")) - 1950L

#' Build Group Lookup from Horizon Groups
#'
#' Constructs a named character vector mapping each individual horizon name to
#' its group name, based on the \code{horizon_groups} list supplied by the user.
#'
#' @param horizon_groups Named list as passed to \code{compute_synchronicity()} or
#'   \code{compute_synchronized_ages()}, or \code{NULL}.
#'
#' @return Named character vector (\code{group_of}) where each element name is a
#'   horizon name and each value is the corresponding group name, or \code{NULL}
#'   when \code{horizon_groups} is \code{NULL}.
#'
#' @keywords internal
build_group_lookup <- function(horizon_groups) {
  if (is.null(horizon_groups)) return(NULL)
  group_of <- character(0)
  for (grp in names(horizon_groups)) {
    group_of[as.character(horizon_groups[[grp]])] <- grp
  }
  group_of
}
