#' Calculate Validation Thresholds from Synchronized Data
#'
#' Calculates empirically-based validation thresholds and confidence levels from synchronized
#' age data for use in subsequent synchronicity testing.
#'
#' @param adjusted_ages Output from \code{set_to_zero()} containing synchronized ages and errors.
#' @param sigma_multiplier Numeric value or named vector specifying standard deviation
#'   multiplier(s) for threshold calculation (default: 2, equivalent to ~95% confidence).
#' @param age_offset Numeric offset value that was added to ages to avoid negative values
#'   (default: \code{bp_datum()}, i.e. the current year minus 1950). This offset is subtracted
#'   before calculating the coefficient of variation to ensure CV is based on true ages, not
#'   offset ages.
#'
#' @return Named list containing:
#'   \itemize{
#'     \item \code{validation_thresholds}: Named numeric vector of age difference thresholds
#'           (as proportions)
#'     \item \code{confidence_levels}: Named numeric vector of corresponding confidence levels
#'           (as proportions)
#'     \item \code{sigma_values}: Named numeric vector of sigma multipliers used per horizon
#'     \item \code{empirical_cv}: Named numeric vector of empirical coefficients of variation
#'     \item \code{age_offset}: The age offset value used during threshold calculation
#'   }
#'
#' @details Calculates the coefficient of variation (CV = sigma/mu) from synchronized ages,
#'   multiplies by user-specified sigma values to generate age difference thresholds, and
#'   converts sigma values to confidence levels using the standard normal distribution.
#'   The age offset is removed before calculating CV to ensure it reflects true age variability.
#'   These thresholds can be directly used in \code{verify_synchronicity()} for validation testing.
#'
#' @examples
#' \dontrun{
#' # Single sigma for all horizons
#' thresholds <- calculate_validation_thresholds(adjusted_ages, sigma_multiplier = 2)
#'
#' # Different sigma per horizon
#' thresholds <- calculate_validation_thresholds(
#'   adjusted_ages,
#'   sigma_multiplier = c(horizon1 = 2, horizon2 = 1, horizon3 = 1.5)
#' )
#'
#' # With custom age offset
#' thresholds <- calculate_validation_thresholds(
#'   adjusted_ages,
#'   sigma_multiplier = 2,
#'   age_offset = 100
#' )
#'
#' # Use in verify_synchronicity
#' verify_synchronicity(
#'   event_stats,
#'   event_names = names(thresholds$validation_thresholds),
#'   confidence_level = thresholds$confidence_levels,
#'   age_difference = thresholds$validation_thresholds
#' )
#' }
#'
#' @export
calculate_validation_thresholds <- function(adjusted_ages,
                                            sigma_multiplier = 2,
                                            age_offset = bp_datum()) {

  # Validate inputs
  if (!is.list(adjusted_ages) || length(adjusted_ages) == 0) {
    stop("adjusted_ages must be a non-empty list")
  }

  if (!is.numeric(sigma_multiplier) || any(sigma_multiplier <= 0)) {
    stop("sigma_multiplier must be positive numeric value(s)")
  }

  if (!is.numeric(age_offset) || length(age_offset) != 1) {
    stop("age_offset must be a single numeric value")
  }

  # Calculate empirical coefficient of variation for each horizon
  # CV = sigma / |mu|
  empirical_cv <- vapply(adjusted_ages, function(x) {
    if (!"adjusted_error" %in% names(x) || !"adjusted_age" %in% names(x)) {
      stop("Each element in adjusted_ages must have 'adjusted_error' and 'adjusted_age' columns")
    }

    # adjusted_ages stores ages relative to the offset-subtracted baseline.
    # Add the offset back to recover the true age (in cal yrs BP) before computing CV.
    true_age <- x$adjusted_age[1] + age_offset

    # CV based on true age, not offset age
    x$adjusted_error[1] / abs(true_age)
  }, numeric(1))

  # Determine sigma multiplier for each horizon
  if (is.null(names(sigma_multiplier))) {
    # Single value for all horizons
    sigma_values <- setNames(rep(sigma_multiplier[1], length(empirical_cv)),
                             names(empirical_cv))
  } else {
    # Named vector - match horizon names
    sigma_values <- vapply(names(empirical_cv), function(horizon) {
      if (horizon %in% names(sigma_multiplier)) {
        sigma_multiplier[[horizon]]
      } else {
        # Default to first value if horizon not specified
        warning(sprintf(
          "Horizon '%s' not found in sigma_multiplier, using default value %.1f",
          horizon, sigma_multiplier[1]
        ))
        sigma_multiplier[1]
      }
    }, numeric(1))
    names(sigma_values) <- names(empirical_cv)
  }

  # Calculate validation thresholds
  # Threshold = 2 * sigma * CV (the factor of 2 accounts for +/- range)
  validation_thresholds <- 2 * sigma_values * empirical_cv

  # Convert sigma to confidence level
  # Using standard normal distribution: P(-sigma < Z < sigma)
  sigma_to_confidence <- function(sigma) {
    2 * pnorm(sigma) - 1
  }

  confidence_levels <- vapply(sigma_values, sigma_to_confidence, numeric(1))

  return(list(
    validation_thresholds = validation_thresholds,
    confidence_levels = confidence_levels,
    sigma_values = sigma_values,
    empirical_cv = empirical_cv,
    age_offset = age_offset
  ))
}

#' Print Validation Thresholds Summary
#'
#' Prints a formatted summary table of the thresholds returned by
#' \code{calculate_validation_thresholds()}.
#'
#' @param thresholds Output list from \code{calculate_validation_thresholds()}.
#'
#' @return Invisibly returns \code{thresholds} (unchanged).
#'
#' @export
print_validation_summary <- function(thresholds) {
  cat("\n=== Validation Thresholds Summary ===\n")
  cat(sprintf("Age offset used: %.0f years\n", thresholds$age_offset))
  cat(sprintf("%-20s %10s %10s %12s %12s\n",
              "Horizon", "Sigma", "Conf. %", "CV %", "Threshold %"))
  cat(strrep("-", 66), "\n")

  for (horizon in names(thresholds$validation_thresholds)) {
    cat(sprintf("%-20s %10.1f %10.1f %12.2f %12.2f\n",
                horizon,
                thresholds$sigma_values[[horizon]],
                thresholds$confidence_levels[[horizon]] * 100,
                thresholds$empirical_cv[[horizon]] * 100,
                thresholds$validation_thresholds[[horizon]] * 100))
  }
  cat(strrep("=", 66), "\n\n")
  invisible(thresholds)
}
