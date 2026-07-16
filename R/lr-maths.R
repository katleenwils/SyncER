#' Asses log-ratio value distribution for normality and centering to allow for synchronicity precision calculation
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

#' MC sampling of age input to calculate log-ratio values
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

#' Calculate minimal synchronicity precision value in log-space
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
compute_minimal_precision_threshold <- function(ABlr, conf_level) {
  ABlr <- ABlr[!is.na(ABlr)]
  if (length(ABlr) == 0) return(NA)
  abs_vals <- sort(abs(ABlr))
  props <- vapply(abs_vals, function(t) mean(abs(ABlr) <= t), numeric(1))
  idx <- which(props >= conf_level)[1]
  if (is.na(idx)) NA else abs_vals[idx]
}

#' Calculate overall synchronicity score using additive log-ratios
#'
#' Calculates overall synchronicity scores (i.e. the probability that all horizons are simultaneously synchronous) 
#' across multiple records using additive log-ratio (ALR) transformation. 
#' Uses the first record alphabetically as reference. The reference record's Monte Carlo resample is drawn once and shared across every non-reference 
#' column, since all comparisons are against the same single (uncertain) reference record; each non-reference
#'   record is independently resampled, reflecting that their age models are unrelated.
#'
#' @param samples_list Named list of numeric vectors containing posterior age samples,
#'   one per record.
#' @param conf_level Numeric value specifying the confidence level for the test (e.g., 0.95).
#' @param age_diff_log_bounds Two-element numeric vector with log-space age difference bounds.
#' @param n_samples Integer number of Monte Carlo draws per record (default: 10000).
#' @param seed Integer random seed for reproducibility (default: 5128).
#'
#' @return Named list containing:
#'   \itemize{
#'     \item \code{overall_score}: Numeric value representing the proportion of joint Monte
#'           Carlo draws (rows of the (k-1)-dimensional ALR distribution) for which
#'           \emph{every} non-reference record is simultaneously within bounds of the shared
#'           reference draw for that row. This is the joint probability that all considered
#'           records are synchronous.
#'     \item \code{overall_precision}: Numeric value representing the smallest age-difference
#'           tolerance (as a proportion) for which \code{overall_score}'s joint
#'           (all-records-simultaneously) requirement would reach \code{conf_level}, or NA
#'           if the distribution is unsuitable
#'     \item \code{n_records}: Integer count of records included in the analysis
#'     \item \code{ref_record}: Character string identifying the reference record used as
#'           denominator in ALR
#'   }
#'
#' @export
compute_overall_synchronicity <- function(samples_list,
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

  # Calculate ALR (Additive Log-Ratio) for each non-reference record.
  # The reference record's resample is drawn ONCE and reused across every
  # non-reference column, since all comparisons share the same single
  # (uncertain) reference age -- resampling it independently per column would
  # treat each comparison as if it had its own separate reference record,
  # erasing the correlation this (k-1)-dimensional distribution is meant to
  # capture.
  n_samples <- min(vapply(samples_list[valid_records], length, integer(1)), n_samples)
  ref_samp  <- sample(ref_samples, n_samples, replace = TRUE)

  alr_matrix <- list()

  for (rec in valid_records) {
    if (rec == ref_record) next  # Skip reference itself

    rec_samples <- samples_list[[rec]]
    rec_samp <- sample(rec_samples, n_samples, replace = TRUE)

    # Calculate log-ratio (shared reference draw as denominator)
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

  # Calculate overall synchronicity score: proportion of joint Monte Carlo draws
  # (rows of the (k-1)-dimensional ALR distribution) for which every
  # non-reference record is simultaneously within bounds of the shared
  # reference draw for that row -- the joint probability that all considered
  # records are synchronous, not merely the marginal average across records.
  all_alr <- as.vector(alr_mat)
  row_within_bounds <- rowSums(alr_mat >= age_diff_log_bounds[1] &
                                  alr_mat <= age_diff_log_bounds[2]) == ncol(alr_mat)
  overall_score <- mean(row_within_bounds, na.rm = TRUE)

  # Calculate overall synchronicity precision: the smallest tolerance t such that
  # at least conf_level proportion of joint draws have ALL k-1 columns
  # simultaneously within [-t, t] -- the same joint (row-wise) requirement used
  # for overall_score, rather than pooling all individual log-ratio values.
  qa <- assess_lr(all_alr)

  # Calculate precision only if distribution is suitable
  if (qa$ok && any(all_alr <= 0, na.rm = TRUE) && any(all_alr >= 0, na.rm = TRUE)) {
    # A row is jointly within +/-t iff every column's absolute value is <= t,
    # i.e. iff the row's largest absolute log-ratio is <= t.
    row_max_abs <- apply(alr_mat, 1, function(row) max(abs(row)))
    abs_vals <- sort(row_max_abs[!is.na(row_max_abs)])
    props <- vapply(abs_vals, function(t) mean(row_max_abs <= t, na.rm = TRUE), numeric(1))
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
