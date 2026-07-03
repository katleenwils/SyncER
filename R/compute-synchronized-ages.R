#' Compute Synchronized Ages
#'
#' Pure-compute core for age synchronization. Performs all method dispatch and age
#' calculation without any plotting or file I/O.
#'
#' @param event_stats List containing processed age data and summaries (output from
#'   \code{process_event_ages()}).
#' @param synchro_results Results from \code{verify_synchronicity()}.
#' @param method Character string or named vector specifying synchronization method(s):
#'   "mean", "mean_fixederror", "Bayesian", "ageofrecord", or "age" (default: NULL).
#' @param horizons Character vector of horizon names (default: NULL).
#' @param excluded_horizons Named or unnamed character vector of horizons to exclude (default: NULL).
#' @param excluded_records Character vector or named list of records to exclude from pooled
#'   calculations (default: NULL).
#' @param age_record Character string for \code{method = "ageofrecord"} (default: NULL).
#' @param age_value Numeric value or named vector for \code{method = "age"} (default: NULL).
#' @param age_error Numeric value or named vector of age errors (default: NULL).
#' @param age_cc Numeric value or named vector of calibration curve codes (default: NULL).
#' @param offset Numeric offset correction value (default: \code{bp_datum()}, i.e.
#'   the current year minus 1950).
#' @param seed Integer random seed for reproducibility (default: 5128).
#' @param bayes_plot_opts Named list of Bayesian plot layout options. Supported keys:
#'   \code{fig_width}, \code{fig_height}, \code{plot_range_sigma}, \code{posterior_lwd},
#'   \code{combined_lwd}, \code{legend_pos}.
#' @param horizon_groups Named list mapping group names to character vectors of horizon names
#'   that belong to each group (e.g., \code{list(tephra = c("tephra1", "tephra1a"))}). When
#'   \code{NULL} (default) each horizon in \code{horizons} is treated as standalone. Use
#'   this to explicitly declare which horizon names across records belong to the same
#'   depositional event.
#'
#' @return Named list with two elements:
#'   \itemize{
#'     \item \code{adjusted_ages}: Named list of data frames (one per horizon) with columns
#'           \code{record}, \code{adjusted_age}, \code{adjusted_error}, \code{cc}
#'     \item \code{bayesian_plot_data}: Named list (one entry per Bayesian horizon) with
#'           raw plot data — pass each element to \code{plot_bayesian_combination()} to
#'           draw to the console or save a PDF
#'   }
#'
#' @export
compute_synchronized_ages <- function(event_stats,
                                      synchro_results,
                                      method = NULL,
                                      horizons = NULL,
                                      excluded_horizons = NULL,
                                      excluded_records = NULL,
                                      age_record = NULL,
                                      age_value = NULL,
                                      age_error = NULL,
                                      age_cc = NULL,
                                      offset = bp_datum(),
                                      seed = 5128,
                                      n_samples = 10000,
                                      bayes_plot_opts = list(),
                                      horizon_groups = NULL) {

  # Set random seed for reproducibility
  set.seed(seed)

  bayes_opts <- modifyList(list(
    fig_width        = 10,
    fig_height       = 7,
    plot_range_sigma = 4,
    posterior_lwd    = 2,
    combined_lwd     = 3,
    legend_pos       = "topright"
  ), bayes_plot_opts)

  # Initialize Bayesian plot data collector
  bayesian_plot_data <- list()

  # Get all unique horizon names from synchronicity results
  all_horizons <- unique(unlist(lapply(event_stats$processed, names)))

  excl_parsed         <- parse_excluded_horizons(excluded_horizons, names(event_stats$processed))
  per_record_excluded <- excl_parsed$per_record_excluded
  global_excluded     <- excl_parsed$global_excluded

  all_real_horizons <- unique(unlist(lapply(event_stats$processed, names)))

  if (length(global_excluded) > 0) {
    global_excluded <- unique(unlist(lapply(global_excluded, function(excl) {
      c(excl, all_real_horizons[grepl(paste0("^", excl, "[a-zA-Z0-9]"), all_real_horizons)])
    })))
  }

  # Remove globally excluded horizons from all records' processed data immediately
  for (rec in names(event_stats$processed)) {
    cols_to_keep <- setdiff(names(event_stats$processed[[rec]]), global_excluded)
    event_stats$processed[[rec]] <- event_stats$processed[[rec]][, cols_to_keep, drop = FALSE]
  }

  # --------------------------------------------------------------------------
  # EXCLUDE NON-SYNCHRONOUS HORIZONS
  # --------------------------------------------------------------------------
  if (length(global_excluded) > 0) {
    all_horizons <- setdiff(all_horizons, global_excluded)
  }

  # Validate age_value and age_error for method "age"
  if ("age" %in% method) {

    # Determine horizons using method "age"
    if (is.null(names(method))) {
      age_horizons <- if (method[1] == "age") horizons else character(0)
    } else {
      age_horizons <- names(method)[method == "age"]
    }

    if (length(age_horizons) > 0) {

      if (is.null(age_value) || is.null(age_error)) {
        stop("For method 'age', both 'age_value' and 'age_error' must be provided.")
      }

      # Check format: must be named vectors (or single values)
      if (is.null(names(age_value)) || is.null(names(age_error))) {
        if (length(age_value) == 1 && length(age_error) == 1) {
          # Single values - OK, will apply to all "age" horizons
        } else if (length(age_value) == length(age_horizons) &&
                   length(age_error) == length(age_horizons)) {
          # Correct length but unnamed - add names
          warning("age_value and age_error are unnamed. Assuming order matches age_horizons.")
          names(age_value) <- age_horizons
          names(age_error) <- age_horizons
        } else {
          stop(sprintf(
            "When using method 'age' for multiple horizons, 'age_value' and 'age_error' must be named vectors.\nExpected names: %s",
            paste(age_horizons, collapse = ", ")
          ))
        }
      } else {
        # Named vectors - check all required horizons present
        missing_value <- setdiff(age_horizons, names(age_value))
        missing_error <- setdiff(age_horizons, names(age_error))

        if (length(missing_value) > 0) {
          stop(sprintf(
            "age_value is missing entries for horizons using method 'age': %s",
            paste(missing_value, collapse = ", ")
          ))
        }

        if (length(missing_error) > 0) {
          stop(sprintf(
            "age_error is missing entries for horizons using method 'age': %s",
            paste(missing_error, collapse = ", ")
          ))
        }
      }
    }
  }

  # Similar validation for mean_fixederror
  if ("mean_fixederror" %in% method) {

    if (is.null(names(method))) {
      fixederror_horizons <- if (method[1] == "mean_fixederror") horizons else character(0)
    } else {
      fixederror_horizons <- names(method)[method == "mean_fixederror"]
    }

    if (length(fixederror_horizons) > 0) {

      if (is.null(age_error)) {
        stop("For method 'mean_fixederror', 'age_error' must be provided.")
      }

      if (is.null(names(age_error))) {
        if (length(age_error) == 1) {
          # Single value applies to all
        } else if (length(age_error) == length(fixederror_horizons)) {
          warning("age_error is unnamed for mean_fixederror. Assuming order matches horizons.")
          names(age_error) <- fixederror_horizons
        } else {
          stop(sprintf(
            "For method 'mean_fixederror', 'age_error' must be a named vector.\nExpected names: %s",
            paste(fixederror_horizons, collapse = ", ")
          ))
        }
      } else {
        missing_error <- setdiff(fixederror_horizons, names(age_error))

        if (length(missing_error) > 0) {
          stop(sprintf(
            "age_error is missing entries for horizons using method 'mean_fixederror': %s",
            paste(missing_error, collapse = ", ")
          ))
        }
      }
    }
  }

  # Initialize output: list of horizon → data.frame
  adjusted_ages <- list()

  # Track which generic groups have been processed
  processed_generic_groups <- c()

  # Remove globally excluded horizons from iteration list
  all_horizons <- setdiff(all_horizons, global_excluded)

  # Build reverse lookup: horizon_name → group_name (from user-supplied horizon_groups)
  group_of <- build_group_lookup(horizon_groups)

  if (!is.null(horizons)) {
    horizons_to_keep <- character(0)
    for (h in horizons) {
      horizons_to_keep <- c(horizons_to_keep, h)
      if (!is.null(horizon_groups) && h %in% names(horizon_groups))
        horizons_to_keep <- c(horizons_to_keep, as.character(horizon_groups[[h]]))
    }
    all_horizons <- intersect(all_horizons, unique(horizons_to_keep))
  }

  # Main processing loop for each horizon
  for (horizon in all_horizons) {

    if (!is.null(group_of) && horizon %in% names(group_of)) {
      event_base        <- group_of[[horizon]]
      relevant_horizons <- as.character(horizon_groups[[event_base]])
      has_generic       <- length(relevant_horizons) > 1
    } else {
      event_base        <- horizon
      relevant_horizons <- horizon
      has_generic       <- FALSE
    }

    # Skip if this generic group was already processed
    if (has_generic && event_base %in% processed_generic_groups) {
      next
    }

    #--------------------------------------------------------------------------
    # GENERIC CASE: Process horizon group once
    #--------------------------------------------------------------------------
    if (has_generic) {

      matching_horizons <- intersect(as.character(horizon_groups[[event_base]]), all_horizons)

      relevant_horizons_existing <- intersect(relevant_horizons, all_horizons)

      # Check that each record contains at most one variant
      records_with_multiple <- names(event_stats$processed)[
        vapply(names(event_stats$processed), function(rec) {
          df <- event_stats$processed[[rec]]
          present <- intersect(relevant_horizons_existing, names(df))

          # Apply exclusions
          excl <- per_record_excluded[[rec]]
          if (!is.null(excl)) present <- setdiff(present, excl)
          present <- setdiff(present, global_excluded)

          length(present) > 1
        }, logical(1))
      ]

      if (length(records_with_multiple) > 0) {
        stop(sprintf(
          "Generic horizon '%s': multiple matching horizons found in records: %s",
          event_base, paste(records_with_multiple, collapse = ", ")
        ))
      }

      # Find records containing ANY of the horizons
      records_with_horizon <- names(event_stats$processed)[
        vapply(names(event_stats$processed), function(rec) {
          df <- event_stats$processed[[rec]]
          present <- intersect(relevant_horizons_existing, names(df))

          excl <- per_record_excluded[[rec]]
          if (!is.null(excl)) present <- setdiff(present, excl)
          present <- setdiff(present, global_excluded)

          length(present) >= 1
        }, logical(1))
      ]

      # Mark this group as processed
      processed_generic_groups <- c(processed_generic_groups, event_base)
      group_name <- event_base

    } else {
      #------------------------------------------------------------------------
      # NON-GENERIC CASE: Single horizon
      #------------------------------------------------------------------------

      relevant_horizons_existing <- horizon

      records_with_horizon <- names(event_stats$processed)[
        vapply(event_stats$processed, function(df) horizon %in% names(df), logical(1))
      ]

      group_name <- horizon
    }

    # Determine method for this horizon
    if (is.null(method)) {
      method_h <- "mean"
    } else if (!is.null(names(method)) && group_name %in% names(method)) {
      method_h <- method[[group_name]]
    } else if (length(method) == 1) {
      method_h <- method[1]
    } else {
      method_h <- "mean"
    }

    # Build per-record original age sample list (needed for most methods)
    original_ages_list <- collect_horizon_samples(
      records_with_horizon, relevant_horizons_existing,
      event_stats$processed, per_record_excluded, global_excluded
    )

    # Resolve per-horizon record exclusions for pooled methods
    excluded_recs_h <- character(0)
    if (!is.null(excluded_records) && method_h %in% c("mean", "mean_fixederror", "Bayesian")) {
      if (is.list(excluded_records) && group_name %in% names(excluded_records)) {
        excluded_recs_h <- as.character(excluded_records[[group_name]])
      } else if (is.character(excluded_records)) {
        excluded_recs_h <- excluded_records
      }
      excluded_present <- intersect(excluded_recs_h, names(original_ages_list))
      if (length(excluded_present) > 0) {
        cat(sprintf("  Excluding record(s) from %s calculation: %s\n",
                    method_h, paste(excluded_present, collapse = ", ")))
      }
    }
    calc_ages_list <- original_ages_list[setdiff(names(original_ages_list), excluded_recs_h)]

    #==========================================================================
    # METHOD: MEAN (uses ALR covariance for error)
    #==========================================================================
    if (method_h == "mean") {

      # Calculate mean of all original ages
      mean_age <- mean(unlist(calc_ages_list))

      # Choose reference record (first alphabetically among non-excluded records)
      ref_record <- sort(names(calc_ages_list))[1]
      ref_samples <- calc_ages_list[[ref_record]]

      if (is.null(ref_samples) || length(ref_samples) == 0) {
        stop(sprintf("Reference record '%s' has no posterior samples for horizon '%s'",
                     ref_record, group_name))
      }

      # Calculate ALR for each non-reference record
      alr_matrix <- list()

      for (rec in names(calc_ages_list)) {
        if (rec == ref_record) next

        rec_samples <- calc_ages_list[[rec]]

        if (is.null(rec_samples) || length(rec_samples) == 0) {
          warning(sprintf("Record '%s' has no posterior samples. Skipping from ALR.", rec))
          next
        }

        # Sample to ensure consistent length
        n_samp <- min(length(ref_samples), length(rec_samples), n_samples)
        ref_samp <- sample(ref_samples, n_samp, replace = TRUE)
        rec_samp <- sample(rec_samples, n_samp, replace = TRUE)

        # Calculate log-ratio
        alr_vec <- log(rec_samp / ref_samp)
        alr_matrix[[rec]] <- alr_vec
      }

      if (length(alr_matrix) == 0) {
        stop(sprintf("No valid ALR comparisons could be computed for horizon '%s'", group_name))
      }

      # Combine into matrix
      alr_mat <- do.call(cbind, alr_matrix)

      # Test multivariate normality
      alr_mean <- colMeans(alr_mat)
      alr_cov_mat <- cov(alr_mat)

      mahal_dist <- mahalanobis(alr_mat, center = alr_mean, cov = alr_cov_mat)
      df <- ncol(alr_mat)

      test_data <- if (length(mahal_dist) > 5000) sample(mahal_dist, 5000) else mahal_dist
      ks_test <- suppressWarnings(ks.test(test_data, "pchisq", df = df))

      is_mvn <- ks_test$p.value >= 0.05
      if (!is_mvn) {
        warning(sprintf(
          "\n  WARNING: log-ratios are not normally distributed for horizon '%s' (p=%.4f)!\n  The error on this age might not accurately represent the age differences between the records.\n",
          group_name, ks_test$p.value
        ))
      }

      # Calculate multivariate statistics
      alr_cov <- cov(alr_mat)
      n_records <- length(calc_ages_list)

      # Total variance in log-space from ALR
      sigma_log_squared <- sum(diag(alr_cov)) / (n_records - 1)
      sigma_log <- sqrt(sigma_log_squared)

      # Convert to coefficient of variation
      cv <- sqrt(exp(sigma_log^2) - 1)

      # Calculate adjusted age and error
      adj_age_value <- mean_age - offset
      adj_error_value <- abs(mean_age) * cv

      # Build output dataframe
      adjusted_ages_df <- data.frame(
        record = records_with_horizon,
        adjusted_age  = rep(adj_age_value, length(records_with_horizon)),
        adjusted_error = rep(adj_error_value, length(records_with_horizon)),
        cc = 0,
        stringsAsFactors = FALSE
      )

      # Assign to all relevant horizon names
      for (horizon_name in relevant_horizons_existing) {
        records_with_this_horizon <- vapply(adjusted_ages_df$record, function(rec) {
          df <- event_stats$processed[[rec]]
          present <- intersect(horizon_name, names(df))

          excl <- per_record_excluded[[rec]]
          if (!is.null(excl)) present <- setdiff(present, excl)
          present <- setdiff(present, global_excluded)

          length(present) > 0
        }, logical(1))

        adjusted_ages[[horizon_name]] <- adjusted_ages_df[records_with_this_horizon, , drop=FALSE]
      }
    }

    #==========================================================================
    # METHOD: MEAN_FIXEDERROR (mean age with user-specified error)
    #==========================================================================
    else if (method_h == "mean_fixederror") {

      # Calculate mean
      mean_age <- mean(unlist(calc_ages_list))

      # Get user-specified error
      if (!is.null(names(age_error)) && group_name %in% names(age_error)) {
        fixed_error <- age_error[[group_name]]
      } else if (is.null(names(age_error)) && length(age_error) == 1) {
        fixed_error <- age_error
      } else {
        stop(sprintf("No age_error provided for horizon '%s' when using method 'mean_fixederror'",
                     group_name))
      }

      # Calculate adjusted values
      adj_age_value <- mean_age - offset
      adj_error_value <- fixed_error

      # Build output
      adjusted_ages_df <- data.frame(
        record = records_with_horizon,
        adjusted_age  = rep(adj_age_value, length(records_with_horizon)),
        adjusted_error = rep(adj_error_value, length(records_with_horizon)),
        cc = 0,
        stringsAsFactors = FALSE
      )

      # Assign to horizon names
      for (horizon_name in relevant_horizons_existing) {
        records_with_this_horizon <- vapply(adjusted_ages_df$record, function(rec) {
          df <- event_stats$processed[[rec]]
          present <- intersect(horizon_name, names(df))

          excl <- per_record_excluded[[rec]]
          if (!is.null(excl)) present <- setdiff(present, excl)
          present <- setdiff(present, global_excluded)

          length(present) > 0
        }, logical(1))

        adjusted_ages[[horizon_name]] <- adjusted_ages_df[records_with_this_horizon, , drop=FALSE]
      }
    }

    #==========================================================================
    # METHOD: BAYESIAN (combine PDFs using Bayesian product)
    #==========================================================================
    else if (method_h == "Bayesian") {

      cat("Using Bayesian MC combined age for", group_name, "\n")

      valid_records <- names(calc_ages_list)
      samples_list  <- unname(calc_ages_list)

      if (length(samples_list) == 0) {
        warning(sprintf("Horizon '%s' had no usable posterior samples. Skipping.", group_name))
        next
      }

      # Shift samples by offset BEFORE combination
      samples_list_shifted <- lapply(samples_list, function(s) s - offset)

      # Bayesian combination
      comb <- combine_pdfs_mc(samples_list_shifted, return_full_pdf = TRUE)
      mu_comb <- comb$mean
      sigma_comb <- comb$sd

      cat(sprintf("Adjusted age: %.1f +/- %.1f cal yrs BP\n\n",
                  mu_comb, sigma_comb))

      # Store raw data for plotting by set_to_zero via plot_bayesian_combination()
      bayesian_plot_data[[group_name]] <- list(
        samples_list_shifted = samples_list_shifted,
        combined_pdf_x       = comb$pdf$x,
        combined_pdf_vals    = comb$pdf$pdf,
        valid_records        = valid_records,
        mu_comb              = mu_comb,
        sigma_comb           = sigma_comb,
        group_name           = group_name,
        bayes_opts           = bayes_opts
      )

      # Store results (all records at this horizon receive the combined age,
      # including any that were excluded from the combination calculation)
      adjusted_ages_df <- data.frame(
        record = records_with_horizon,
        adjusted_age   = rep(mu_comb, length(records_with_horizon)),
        adjusted_error = rep(max(sigma_comb, 1e-6), length(records_with_horizon)),
        cc = 0,
        stringsAsFactors = FALSE
      )

      for (horizon_name in relevant_horizons_existing) {
        records_with_this_horizon <- vapply(adjusted_ages_df$record, function(rec) {
          df <- event_stats$processed[[rec]]
          present <- intersect(horizon_name, names(df))

          excl <- per_record_excluded[[rec]]
          if (!is.null(excl)) present <- setdiff(present, excl)
          present <- setdiff(present, global_excluded)

          length(present) > 0
        }, logical(1))

        adjusted_ages[[horizon_name]] <- adjusted_ages_df[records_with_this_horizon, , drop=FALSE]
      }
    }

    #==========================================================================
    # METHOD: AGEOFRECORD (use age from specific record)
    #==========================================================================
    else if (method_h == "ageofrecord") {

      if (is.null(age_record)) {
        stop("age_record must be specified.\n")
      }

      if (!(age_record %in% records_with_horizon)) {
        warning(sprintf("%s is not present in %s. Skipping this horizon...\n",
                        group_name, age_record))
        next
      } else {
        cat("Using age of", age_record, "for", group_name, "\n")
      }

      # Extract retained horizon for reference record
      df_ref_proc <- event_stats$processed[[age_record]]
      df_ref_sum  <- event_stats$summaries[[age_record]]

      present_ref <- intersect(relevant_horizons_existing, names(df_ref_proc))
      present_ref <- setdiff(present_ref, global_excluded)

      if (length(present_ref) != 1) {
        stop(sprintf(
          "Reference record '%s' does not contain exactly one retained horizon for group %s",
          age_record, group_name
        ))
      }

      retained_h <- present_ref

      # Extract ages from processed table
      age0_samples <- df_ref_proc[[retained_h]]
      age0 <- mean(age0_samples)

      # Extract sigma from summaries
      sigma_col <- paste0(retained_h, "_sigma")

      if (sigma_col %in% names(df_ref_sum)) {
        sigma0 <- df_ref_sum[[sigma_col]]
      } else {
        warning(sprintf(
          "Sigma column %s not found in summaries for record %s. Estimating SD from samples.",
          sigma_col, age_record
        ))
        sigma0 <- sd(age0_samples)
        if (is.na(sigma0) || sigma0 == 0) sigma0 <- 1
      }

      cat(sprintf("Adjusted age: %.1f +/- %.1f cal yrs BP\n\n",
                  age0 - offset, sigma0))

      # Build data frame
      adjusted_ages_df <- data.frame(
        record         = records_with_horizon,
        adjusted_age   = rep(age0 - offset, length(records_with_horizon)),
        adjusted_error = rep(sigma0, length(records_with_horizon)),
        cc = 0,
        stringsAsFactors = FALSE
      )

      # Assign per-horizon
      for (horizon_name in relevant_horizons_existing) {
        records_with_this_horizon <- vapply(adjusted_ages_df$record, function(rec) {
          df_proc <- event_stats$processed[[rec]]
          present <- intersect(horizon_name, names(df_proc))

          excl <- per_record_excluded[[rec]]
          if (!is.null(excl)) present <- setdiff(present, excl)
          present <- setdiff(present, global_excluded)

          length(present) > 0
        }, logical(1))

        adjusted_ages[[horizon_name]] <- adjusted_ages_df[records_with_this_horizon, , drop=FALSE]
      }
    }

    #==========================================================================
    # METHOD: AGE (use user-specified age and error)
    #==========================================================================
    else if (method_h == "age") {

      # Determine age and error for this horizon
      if (!is.null(names(age_value)) && group_name %in% names(age_value)) {
        age0 <- age_value[[group_name]]
        sigma0 <- if (!is.null(age_error) && group_name %in% names(age_error)) {
          age_error[[group_name]]
        } else {
          0
        }
      } else if (is.null(names(age_value))) {
        # Single universal age
        age0 <- age_value
        sigma0 <- if (!is.null(age_error)) age_error else 0
      } else {
        stop(sprintf("No age_value provided for horizon '%s'", group_name))
      }

      # Determine cc value
      cc0 <- if (!is.null(names(age_cc)) && group_name %in% names(age_cc)) {
        age_cc[[group_name]]
      } else if (!is.null(age_cc) && is.null(names(age_cc))) {
        age_cc
      } else {
        0
      }

      cat(sprintf("Using age of %.1f +/- %.1f (cc = %d) for %s\n\n", age0, sigma0, cc0, group_name))

      adjusted_ages_df <- data.frame(
        record         = records_with_horizon,
        adjusted_age   = rep(age0,   length(records_with_horizon)),
        adjusted_error = rep(sigma0, length(records_with_horizon)),
        cc             = rep(cc0,    length(records_with_horizon)),
        stringsAsFactors = FALSE
      )

      # Assign per-horizon
      for (horizon_name in relevant_horizons_existing) {
        records_with_this_horizon <- vapply(adjusted_ages_df$record, function(rec) {
          df <- event_stats$processed[[rec]]
          present <- intersect(horizon_name, names(df))

          excl <- per_record_excluded[[rec]]
          if (!is.null(excl)) present <- setdiff(present, excl)
          present <- setdiff(present, global_excluded)

          length(present) > 0
        }, logical(1))

        adjusted_ages[[horizon_name]] <- adjusted_ages_df[records_with_this_horizon, , drop=FALSE]
      }
    }
  }

  return(list(adjusted_ages = adjusted_ages, bayesian_plot_data = bayesian_plot_data))
}

#' Compute Synchronized Ages and Save Outputs (I/O Wrapper)
#'
#' Calls \code{compute_synchronized_ages()} and saves Bayesian combination plots to PDF.
#'
#' @inheritParams compute_synchronized_ages
#' @param output_dir Character string specifying the directory in which the
#'   "synchro_test" subfolder will be created for Bayesian plot PDFs.
#' @param seed Integer random seed for reproducibility (default: 5128).
#' @param bayes_plot_opts Named list of Bayesian plot layout options passed to
#'   \code{compute_synchronized_ages()}. Supported keys: \code{fig_width},
#'   \code{fig_height}, \code{plot_range_sigma}, \code{posterior_lwd},
#'   \code{combined_lwd}, \code{legend_pos}.
#'
#' @return Named list of synchronized ages (one element per horizon), as returned
#'   by \code{compute_synchronized_ages()$adjusted_ages}.
#'
#' @export
set_to_zero <- function(event_stats,
                        synchro_results,
                        output_dir,
                        method = NULL,
                        horizons = NULL,
                        excluded_horizons = NULL,
                        excluded_records = NULL,
                        age_record = NULL,
                        age_value = NULL,
                        age_error = NULL,
                        age_cc = NULL,
                        offset = bp_datum(),
                        seed = 5128,
                        n_samples = 10000,
                        bayes_plot_opts = list(),
                        horizon_groups = NULL) {

  result <- compute_synchronized_ages(
    event_stats       = event_stats,
    synchro_results   = synchro_results,
    method            = method,
    horizons          = horizons,
    excluded_horizons = excluded_horizons,
    excluded_records  = excluded_records,
    age_record        = age_record,
    age_value         = age_value,
    age_error         = age_error,
    age_cc            = age_cc,
    offset            = offset,
    seed              = seed,
    n_samples         = n_samples,
    bayes_plot_opts   = bayes_plot_opts,
    horizon_groups    = horizon_groups
  )

  for (pd in result$bayesian_plot_data) {
    plot_bayesian_combination(pd, output_dir = output_dir)
  }

  invisible(result$adjusted_ages)
}

#' Assign Ages to Excluded Horizons
#'
#' Assigns synchronized ages to horizons that were excluded from the main synchronization process.
#'
#' @param adjusted_ages Output from \code{set_to_zero()} containing synchronized ages.
#' @param excluded_horizons Named or unnamed character vector specifying excluded horizons
#'   (same format as in \code{set_to_zero()}).
#' @param event_stats Output from \code{process_event_ages()} containing record metadata.
#'
#' @return Named list with one element per excluded horizon, each containing a data frame
#'   with columns: \code{record}, \code{adjusted_age}, \code{adjusted_error}.
#'
#' @details For each excluded horizon, assigns the age and error from the first available
#'   non-excluded horizon in the same record. Parses exclusions using the same logic as
#'   \code{set_to_zero()} (named for per-record, unnamed for global).
#'
#' @export
assign_excluded_ages <- function(adjusted_ages,
                                 excluded_horizons,
                                 event_stats) {

  excl_parsed         <- parse_excluded_horizons(excluded_horizons, names(event_stats$processed))
  per_record_excluded <- excl_parsed$per_record_excluded
  global_excluded     <- excl_parsed$global_excluded

  all_records <- names(event_stats$processed)

  # Final list of exclusions per record
  exclusions_by_record <- lapply(all_records, function(rec) {
    unique(c(global_excluded, per_record_excluded[[rec]]))
  })
  names(exclusions_by_record) <- all_records

  # Build export structure
  out <- list()

  # For every record
  for (rec in all_records) {

    rec_exclusions <- exclusions_by_record[[rec]]
    if (length(rec_exclusions) == 0) next  # No exclusions for this record

    present_horizons <- names(event_stats$processed[[rec]])

    # Process each excluded horizon
    for (ex_h in rec_exclusions) {

      if (!(ex_h %in% present_horizons)) next  # Excluded horizon not present in record

      # Find all OTHER horizons in the record that have adjusted ages
      other_hz <- setdiff(present_horizons, ex_h)
      other_hz <- other_hz[other_hz %in% names(adjusted_ages)]

      if (length(other_hz) == 0) {
        warning(sprintf(
          "Record '%s': excluded horizon '%s' has no other horizons with adjusted ages.",
          rec, ex_h
        ))
        next
      }

      # Take first reference horizon (alphabetically)
      ref <- other_hz[1]
      ref_df <- adjusted_ages[[ref]]
      ref_row <- ref_df[ref_df$record == rec, ]

      if (nrow(ref_row) == 0) next

      # Initialize output list for this horizon if needed
      if (is.null(out[[ex_h]])) {
        out[[ex_h]] <- data.frame(
          record = character(),
          adjusted_age = numeric(),
          adjusted_error = numeric(),
          stringsAsFactors = FALSE
        )
      }

      # Add row
      out[[ex_h]] <- rbind(out[[ex_h]], data.frame(
        record = rec,
        adjusted_age = ref_row$adjusted_age,
        adjusted_error = ref_row$adjusted_error,
        stringsAsFactors = FALSE
      ))

      cat(sprintf(
        "Assigned %.1f +/- %.1f to excluded '%s' in record '%s' using '%s'\n",
        ref_row$adjusted_age, ref_row$adjusted_error, ex_h, rec, ref
      ))
    }
  }

  return(out)
}

#' Bayesian Combination of PDFs Using Monte Carlo Samples
#'
#' Combines multiple posterior age distributions using Bayesian product of probability
#' density functions.
#'
#' @param samples_list List of numeric vectors containing posterior samples from different records.
#' @param return_full_pdf Logical indicating whether to return the full combined PDF
#'   (default: FALSE).
#' @param n_grid Integer specifying number of grid points for kernel density estimation
#'   (default: 2000).
#' @param bw Character string or numeric value specifying bandwidth method for KDE
#'   (default: "nrd0").
#'
#' @return List containing:
#'   \itemize{
#'     \item \code{mean}: Numeric mean of the combined distribution
#'     \item \code{sd}: Numeric standard deviation of the combined distribution
#'     \item \code{pdf}: (if \code{return_full_pdf = TRUE}) List with \code{x} (grid points)
#'           and \code{pdf} (density values)
#'   }
#'
#' @details Uses kernel density estimation to construct PDFs from samples, multiplies densities
#'   in log-space to avoid underflow, normalizes using trapezoidal integration, and calculates
#'   moment-matched statistics. This is the recommended method for combining ages from multiple
#'   records with overlapping distributions.
#'
#' @export
combine_pdfs_mc <- function(samples_list,
                            return_full_pdf = FALSE,
                            n_grid = 2000,
                            bw = "nrd0") {

  # Remove null/empty entries from samples list
  samples_list <- samples_list[vapply(samples_list, function(x) {
    !is.null(x) && length(x) > 1
  }, logical(1))]

  if (length(samples_list) == 0) {
    stop("No valid sample distributions provided.")
  }

  # Create kernel density estimates for each distribution
  kde_list <- lapply(samples_list, function(samp) {
    density(samp, bw = bw, n = n_grid)
  })

  # Common evaluation grid (use first KDE's grid)
  x_grid <- kde_list[[1]]$x

  # Extract density values for each KDE at grid points
  y_mat <- sapply(kde_list, function(k) k$y)

  # Bayesian product PDF: multiply densities (in log-space to avoid underflow)
  eps <- 1e-300  # Small constant to prevent log(0)
  log_prod <- rowSums(log(y_mat + eps))

  # Exponentiate back to get product PDF
  prod_pdf <- exp(log_prod)

  # Normalize to integrate to 1 using trapezoidal rule
  prod_pdf <- prod_pdf / pracma::trapz(x_grid, prod_pdf)

  # Compute moment-matched normal approximation
  # Mean: integral x*p(x) dx
  mean_comb <- pracma::trapz(x_grid, x_grid * prod_pdf)

  # Variance: integral (x-mu)^2*p(x) dx
  var_comb <- pracma::trapz(x_grid, (x_grid - mean_comb)^2 * prod_pdf)
  sd_comb <- sqrt(var_comb)

  # Prepare output
  out <- list(
    mean = mean_comb,
    sd   = sd_comb
  )

  # Include full PDF if requested
  if (return_full_pdf) {
    out$pdf <- list(
      x   = x_grid,
      pdf = prod_pdf
    )
  }

  return(out)
}
