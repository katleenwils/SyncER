#' Assess Log-Ratio Vector for Distributional Quality
#'
#' Tests whether a log-ratio vector has a normal shape (via \code{shape_ok()}) and
#' is sufficiently centred (mean within half a standard deviation of zero). Both
#' conditions must hold for the precision estimate derived from the distribution to
#' be trustworthy.
#'
#' @param ABlr Numeric vector of log-ratio values.
#'
#' @return Named list with elements:
#'   \itemize{
#'     \item \code{ok}: Logical; \code{TRUE} when both shape and centring criteria pass.
#'     \item \code{lr_mean}: Mean of \code{ABlr} (NAs removed).
#'     \item \code{lr_sd}: Standard deviation of \code{ABlr} (NAs removed).
#'   }
#'
#' @keywords internal
assess_lr <- function(ABlr) {
  lr_mean <- mean(ABlr, na.rm = TRUE)
  lr_sd   <- sd(ABlr,   na.rm = TRUE)
  list(
    ok      = isTRUE(shape_ok(ABlr)) && abs(lr_mean) <= 0.5 * lr_sd,
    lr_mean = lr_mean,
    lr_sd   = lr_sd
  )
}

#' Sample Pairwise Log-Ratio
#'
#' Draws \code{n_samples} values with replacement from each of two age-sample
#' vectors and returns the vector of log-ratios log(A / B).
#'
#' @param samples_a Numeric vector of posterior age samples for record A.
#' @param samples_b Numeric vector of posterior age samples for record B.
#' @param n_samples Integer number of Monte Carlo draws (default: 10000).
#'
#' @return Numeric vector of length \code{n_samples}.
#'
#' @keywords internal
compute_pairwise_lr <- function(samples_a, samples_b, n_samples = 10000) {
  as.vector(log(sample(samples_a, size = n_samples, replace = TRUE) /
                sample(samples_b, size = n_samples, replace = TRUE)))
}

#' Find Precision Threshold for a Log-Ratio Distribution
#'
#' Returns the smallest absolute log-ratio value \eqn{t} such that at least
#' \code{conf_level} of \code{ABlr} values fall within \eqn{[-t, t]}.
#'
#' @param ABlr Numeric vector of log-ratio values.
#' @param conf_level Numeric confidence level (e.g., 0.95).
#'
#' @return Numeric threshold, or \code{NA} if \code{ABlr} is empty or the
#'   desired coverage cannot be achieved.
#'
#' @keywords internal
find_precision_threshold <- function(ABlr, conf_level) {
  ABlr <- ABlr[!is.na(ABlr)]
  if (length(ABlr) == 0) return(NA)
  abs_vals <- sort(abs(ABlr))
  props <- vapply(abs_vals, function(t) mean(abs(ABlr) <= t), numeric(1))
  idx <- which(props >= conf_level)[1]
  if (is.na(idx)) NA else abs_vals[idx]
}

#' Calculate Overall Synchronicity Score Using ALR Approach
#'
#' Calculates overall synchronicity scores across multiple records using additive
#' log-ratio (ALR) transformation.
#'
#' @param samples_list Named list of numeric vectors containing posterior age samples,
#'   one per record.
#' @param conf_level Numeric value specifying the confidence level for the test (e.g., 0.95).
#' @param age_diff_log_bounds Two-element numeric vector with log-space age difference bounds.
#'
#' @return Named list containing:
#'   \itemize{
#'     \item \code{overall_score}: Numeric value representing the proportion of ALR samples
#'           within age difference bounds
#'     \item \code{overall_precision}: Numeric value representing the synchronicity precision
#'           (as proportion) or NA if distribution is unsuitable
#'     \item \code{n_records}: Integer count of records included in the analysis
#'     \item \code{ref_record}: Character string identifying the reference record used as
#'           denominator in ALR
#'   }
#'
#' @details Uses the first record alphabetically as reference. Calculates ALR for all other
#'   records, tests for multivariate normality using Mahalanobis distances, and computes
#'   overall synchronicity metrics. Issues warnings if multivariate normality assumptions
#'   are violated.
#'
#' @export
calculate_overall_synchronicity <- function(samples_list,
                                            conf_level,
                                            age_diff_log_bounds,
                                            n_samples = 10000,
                                            seed = 5128) {

  set.seed(seed)

  # Get valid record names
  valid_records <- names(samples_list)
  n_records <- length(valid_records)

  # Need at least 2 records to compare
  if (n_records < 2) {
    return(list(overall_score = NA, overall_precision = NA,
                n_records = n_records, ref_record = NA))
  }

  # Choose reference record (first alphabetically for consistency)
  ref_record <- sort(valid_records)[1]
  ref_samples <- samples_list[[ref_record]]

  # Calculate ALR (Additive Log-Ratio) for each non-reference record
  alr_matrix <- list()

  for (rec in valid_records) {
    if (rec == ref_record) next  # Skip reference itself

    rec_samples <- samples_list[[rec]]

    # Sample with replacement to standardize sample size
    n_samples <- min(length(ref_samples), length(rec_samples), n_samples)
    ref_samp <- sample(ref_samples, n_samples, replace = TRUE)
    rec_samp <- sample(rec_samples, n_samples, replace = TRUE)

    # Calculate log-ratio (reference as denominator)
    alr_vec <- log(rec_samp / ref_samp)
    alr_matrix[[rec]] <- alr_vec
  }

  # Check if we have any comparisons
  if (length(alr_matrix) == 0) {
    return(list(overall_score = NA, overall_precision = NA,
                n_records = n_records, ref_record = ref_record))
  }

  # Combine ALR vectors into matrix (each column = one comparison)
  alr_mat <- do.call(cbind, alr_matrix)

  # Test for multivariate normality using Mahalanobis distances
  alr_mean <- colMeans(alr_mat)
  alr_cov_mat <- cov(alr_mat)

  # Calculate Mahalanobis distance for each observation
  mahal_dist <- mahalanobis(alr_mat, center = alr_mean, cov = alr_cov_mat)
  df <- ncol(alr_mat)

  # Test if distances follow chi-squared distribution (expected under MVN)
  test_data <- if (length(mahal_dist) > 5000) sample(mahal_dist, 5000) else mahal_dist
  ks_test <- suppressWarnings(ks.test(test_data, "pchisq", df = df))

  # Calculate overall synchronicity score
  # Score = proportion of samples within age difference bounds
  all_alr <- as.vector(alr_mat)
  overall_score <- mean(all_alr >= age_diff_log_bounds[1] &
                          all_alr <= age_diff_log_bounds[2], na.rm = TRUE)

  # Calculate overall synchronicity precision
  qa <- assess_lr(all_alr)

  # Calculate precision only if distribution is suitable
  if (qa$ok && any(all_alr <= 0, na.rm = TRUE) && any(all_alr >= 0, na.rm = TRUE)) {
    # Find threshold where conf_level proportion of samples fall within ±t
    abs_vals <- sort(abs(all_alr[!is.na(all_alr)]))
    props <- vapply(abs_vals, function(t) mean(abs(all_alr) <= t, na.rm = TRUE), numeric(1))
    idx <- which(props >= conf_level)[1]

    if (!is.na(idx)) {
      t_needed <- abs_vals[idx]
      overall_precision <- exp(t_needed) - 1
    } else {
      overall_precision <- NA
    }
  } else {
    overall_precision <- NA
  }

  return(list(
    overall_score = overall_score,
    overall_precision = overall_precision,
    n_records = n_records,
    ref_record = ref_record
  ))
}
