#' Retrieve confidence level, age difference, and log-ratio bounds per horizon
#'
#' Extracts \code{conf_level_h} and \code{age_diff_h} for a given horizon from
#' scalar or named-vector inputs, converts an absolute age difference (>= 1) to a
#' relative one using the mean age from \code{summaries}, and computes symmetric
#' log-space bounds.
#'
#' @param horizon Character string identifying the horizon (key for named lookups).
#' @param confidence_level Numeric scalar or named numeric vector.
#' @param age_difference Numeric scalar or named numeric vector.
#' @param summaries Named list of per-record summary data frames from
#'   \code{process_event_ages()$summaries}. Used only when \code{age_difference >= 1}
#'   to derive a relative proportion from the mean age.
#'
#' @return Named list with:
#'   \itemize{
#'     \item \code{conf_level_h}: Effective confidence level.
#'     \item \code{age_diff_h}: Effective (relative) age difference.
#'     \item \code{age_diff_log_bounds}: Two-element numeric \code{c(-t, t)} where
#'           \code{t = log(1 + age_diff_h)}.
#'   }
#'
#' @keywords internal
get_horizon_thresholds <- function(horizon, confidence_level, age_difference, summaries) {

  if (!is.null(names(confidence_level)) && horizon %in% names(confidence_level)) {
    conf_level_h <- confidence_level[[horizon]]
  } else {
    conf_level_h <- confidence_level[[1]]
  }

  if (!is.null(names(age_difference)) && horizon %in% names(age_difference)) {
    age_diff_h <- age_difference[[horizon]]
  } else {
    unnamed_idx <- which(!nzchar(names(age_difference)))
    age_diff_h <- if (length(unnamed_idx) > 0) {
      age_difference[[unnamed_idx[length(unnamed_idx)]]]
    } else {
      age_difference[[1]]
    }
  }

  if (age_diff_h >= 1) {
    mean_ages <- vapply(names(summaries), function(rec) {
      sum_cols <- names(summaries[[rec]])
      mean_col <- sum_cols[grepl(paste0("^", horizon, ".*_mean$"), sum_cols)]
      if (length(mean_col) == 0) return(NA_real_)
      mean(unlist(summaries[[rec]][mean_col]), na.rm = TRUE)
    }, numeric(1))
    horizon_mean <- mean(mean_ages, na.rm = TRUE)
    if (is.na(horizon_mean) || horizon_mean <= 0) {
      warning(sprintf(
        "Could not compute valid mean age for horizon '%s' to convert absolute age_difference. Using value as-is.",
        horizon
      ))
    } else {
      age_diff_h <- age_diff_h / horizon_mean
    }
  }

  t <- log(1 + age_diff_h)
  list(conf_level_h = conf_level_h, age_diff_h = age_diff_h, age_diff_log_bounds = c(-t, t))
}

#' Compute a log-ratio values for a single pair-wise comparison
#'
#' Core computation for one record pair within \code{compute_synchronicity_values()}.
#' Draws \code{n_samples} log-ratios, scores them against the age-difference bounds,
#' assesses distribution shape, estimates precision, and assembles the statistics and
#' visualization rows.
#'
#' @param col1 Numeric vector of posterior age samples for record \code{rec_i}.
#' @param col2 Numeric vector of posterior age samples for record \code{rec_j}.
#' @param rec_i Character string identifying record i.
#' @param rec_j Character string identifying record j.
#' @param var_i Column name used from record i. Pass \code{NULL} (default) for single-horizon
#'   comparisons; \code{variant1} in the stats row is then set to \code{NA}.
#' @param var_j Column name used from record j. Same semantics as \code{var_i}.
#' @param horizon Character string identifying the event/horizon being compared.
#' @param thresholds Named list as returned by \code{get_horizon_thresholds()}.
#' @param summaries Named list of per-record summary data frames from
#'   \code{process_event_ages()$summaries} (for min/max in the viz row).
#' @param n_samples Integer; number of Monte Carlo samples (default: 10000).
#'
#' @return Named list with elements \code{lr_name}, \code{variant_key}, \code{ABlr},
#'   \code{score}, \code{precision}, \code{stats_row}, \code{viz_row}.
#'
#' @keywords internal
compare_pair <- function(col1, col2,
                         rec_i, rec_j,
                         var_i = NULL, var_j = NULL,
                         horizon, thresholds, summaries,
                         n_samples = 10000) {

  conf_level_h        <- thresholds$conf_level_h
  age_diff_h          <- thresholds$age_diff_h
  age_diff_log_bounds <- thresholds$age_diff_log_bounds

  is_variant <- !is.null(var_i) && !is.null(var_j)
  col_i <- if (is_variant) var_i else horizon
  col_j <- if (is_variant) var_j else horizon

  ABlr  <- compute_pairwise_lr(col1, col2, n_samples)
  score <- mean(ABlr >= age_diff_log_bounds[1] & ABlr <= age_diff_log_bounds[2], na.rm = TRUE)

  qa        <- assess_lr(ABlr)
  precision <- if (qa$ok) {
    t_needed <- compute_minimal_precision_threshold(ABlr, conf_level_h)
    if (!is.na(t_needed)) exp(t_needed) - 1 else NA_real_
  } else {
    NA_real_
  }

  alpha     <- (1 - conf_level_h) / 2
  z         <- qnorm(1 - alpha)
  conf_low  <- qa$lr_mean - z * qa$lr_sd
  conf_high <- qa$lr_mean + z * qa$lr_sd

  lr_name     <- paste(rec_i, rec_j, "lr", sep = "_")
  variant_key <- if (is_variant) paste(col_i, col_j, sep = "_") else horizon

  stats_row <- data.frame(
    Event                   = horizon,
    record1                 = rec_i,
    record2                 = rec_j,
    variant1                = if (is_variant) col_i else NA_character_,
    variant2                = if (is_variant) col_j else NA_character_,
    mean_LR                 = qa$lr_mean,
    SD_LR                   = qa$lr_sd,
    confidence_level        = conf_level_h,
    age_difference          = age_diff_h,
    synchronicity_score     = score,
    synchronicity_precision = precision,
    stringsAsFactors        = FALSE
  )

  viz_row <- data.frame(
    Event    = horizon,
    record1  = rec_i,
    record2  = rec_j,
    variant1 = if (is_variant) col_i else NA_character_,
    variant2 = if (is_variant) col_j else NA_character_,
    Amin     = summaries[[rec_i]][[paste(col_i, "min", sep = "_")]],
    Amax     = summaries[[rec_i]][[paste(col_i, "max", sep = "_")]],
    Bmin     = summaries[[rec_j]][[paste(col_j, "min", sep = "_")]],
    Bmax     = summaries[[rec_j]][[paste(col_j, "max", sep = "_")]],
    ABlr     = I(list(ABlr)),
    col1     = I(list(col1)),
    col2     = I(list(col2)),
    age_diff_log_bounds_low  = age_diff_log_bounds[1],
    age_diff_log_bounds_high = age_diff_log_bounds[2],
    conf_log_bounds_low      = conf_low,
    conf_log_bounds_high     = conf_high,
    conf_level               = conf_level_h,
    age_diff                 = age_diff_h,
    stringsAsFactors         = FALSE
  )

  list(lr_name = lr_name, variant_key = variant_key, ABlr = ABlr,
       score = score, precision = precision, stats_row = stats_row, viz_row = viz_row)
}

#' Seperate per-record excluded horizons from globally excluded horizons
#'
#' Normalizes the \code{nonsynchro_horizons} argument into two tidy outputs used
#' downstream to filter age columns. Accepts a named list (CASE 1), a named
#' character vector (CASE 2), or an unnamed character vector (global exclusion).
#'
#' @param nonsynchro_horizons Named list, named character vector, or unnamed character
#'   vector. Named entries are matched against \code{all_records} by substring;
#'   unnamed entries (or elements keyed \code{"global"}) are global.
#' @param all_records Character vector of all record names.
#'
#' @return Named list with:
#'   \itemize{
#'     \item \code{per_record_excluded}: Named list of per-record horizon exclusions.
#'     \item \code{global_excluded}: Character vector of globally excluded horizons.
#'   }
#'
#' @keywords internal
parse_excluded_horizons <- function(nonsynchro_horizons, all_records) {

  per_record_excluded <- list()
  global_excluded <- character(0)

  if (is.null(nonsynchro_horizons))
    return(list(per_record_excluded = per_record_excluded, global_excluded = global_excluded))

  add_per_record <- function(key, vals) {
    rec_match <- all_records[grepl(key, all_records)]
    if (length(rec_match) == 1) {
      per_record_excluded[[rec_match]] <<-
        unique(c(per_record_excluded[[rec_match]], vals))
    } else if (length(rec_match) > 1) {
      stop(sprintf("Ambiguous exclusion key '%s' matched multiple records: %s",
                   key, paste(rec_match, collapse = ", ")))
    } else {
      warning(sprintf("Exclusion key '%s' did not match any record name.", key))
    }
  }

  if (is.list(nonsynchro_horizons)) {
    for (key in names(nonsynchro_horizons)) {
      vals <- as.character(nonsynchro_horizons[[key]])
      if (is.na(key) || !nzchar(key) || key == "global") {
        global_excluded <- c(global_excluded, vals)
      } else {
        add_per_record(key, vals)
      }
    }
  } else {
    named_keys <- names(nonsynchro_horizons)
    named_vals <- as.character(nonsynchro_horizons)
    if (!is.null(named_keys)) {
      for (i in seq_along(named_keys)) {
        if (nzchar(named_keys[i])) {
          add_per_record(named_keys[i], named_vals[i])
        } else {
          global_excluded <- c(global_excluded, named_vals[i])
        }
      }
    } else {
      global_excluded <- named_vals
    }
  }

  list(per_record_excluded = per_record_excluded, global_excluded = global_excluded)
}

#' Retreve posterior age Samples for a horizon group across all records
#'
#' Builds a named list of posterior age sample vectors for the specified horizon
#' column names, after applying per-record and global exclusions.
#'
#' @param records_with_horizon Character vector of record names to search.
#' @param relevant_horizons Character vector of column names belonging to the horizon group.
#' @param processed Named list of per-record data frames from
#'   \code{process_event_ages()$processed}.
#' @param per_record_excluded Named list of per-record exclusions from
#'   \code{parse_excluded_horizons()}.
#' @param global_excluded Character vector of globally excluded horizons from
#'   \code{parse_excluded_horizons()}.
#'
#' @return Named list of numeric vectors (posterior samples), keyed by record name.
#'   Records with no valid samples after filtering are omitted.
#'
#' @keywords internal
get_horizon_posteriors <- function(records_with_horizon, relevant_horizons,
                                    processed, per_record_excluded, global_excluded) {
  result <- list()
  for (rec in records_with_horizon) {
    df      <- processed[[rec]]
    present <- intersect(relevant_horizons, names(df))
    excl    <- per_record_excluded[[rec]]
    if (!is.null(excl)) present <- setdiff(present, excl)
    present <- setdiff(present, global_excluded)

    if (length(present) > 1)
      stop(sprintf("Record '%s' still has multiple horizons after exclusion filtering: %s",
                   rec, paste(present, collapse = ", ")))

    if (length(present) == 1) {
      samp <- df[[present]]
      if (!is.null(samp) && length(samp) > 1)
        result[[rec]] <- as.numeric(samp)
    }
  }
  result
}
