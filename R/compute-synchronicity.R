#' Compute Synchronicity of Event Deposits
#'
#' Pure-compute core: performs all pairwise comparisons and overall synchronicity
#' calculations without any file I/O or plotting.
#'
#' @param event_stats Output from \code{process_event_ages()}.
#' @param event_names Character vector of event names to test.
#' @param confidence_level Numeric value or named vector giving the confidence level (a ratio)
#'   at which synchronicity is tested (default: \code{0.95}). Only correlations whose
#'   synchronicity score exceeds this value pass the test -- e.g. with
#'   \code{confidence_level = 0.95}, only pairs scoring higher than 0.95 are counted as passes.
#'   Supply a named vector to use different confidence levels per horizon.
#' @param age_difference Numeric value or named vector giving the maximum age difference between
#'   two records that is still considered synchronous (default: \code{0.05}). Interpreted in one
#'   of two ways depending on its magnitude:
#'   \itemize{
#'     \item Values \code{< 1} are treated as a \emph{relative} (proportional) age difference.
#'           For example, \code{age_difference = 0.05} means the synchronicity score reports the
#'           confidence with which the maximum age difference between the two records is 5\%.
#'     \item Values \code{>= 1} are treated as an \emph{absolute} age difference in years
#'           (e.g. \code{10} for 10 yr), converted internally to a relative difference using the
#'           horizon's mean age.
#'   }
#'   Supply a named vector to set different values for different horizons; include one unnamed
#'   (or last) value to act as the default for any horizon not named.
#' @param horizon_groups Named list mapping group names to character vectors of horizon names
#'   that belong to each group (e.g., \code{list(tephra = c("tephra1", "tephra1a"))}). When
#'   \code{NULL} (default) each horizon in \code{event_names} is treated as standalone. Use
#'   this to explicitly declare which horizon names across records belong to the same
#'   depositional event.
#'
#' @param n_samples Integer number of Monte Carlo samples drawn per pairwise comparison
#'   (default: 10000).
#' @param seed Integer random seed for reproducibility (default: 5128).
#'
#' @return Named list with elements:
#'   \itemize{
#'     \item \code{results}: Pairwise log-ratio results
#'     \item \code{overall_scores}: Overall synchronicity scores per horizon group
#'     \item \code{viz_data}: Data frame of visualization inputs (one row per pairwise comparison)
#'     \item \code{all_scores}: Data frame of all pairwise scores
#'     \item \code{summary}: List with total, positives, negatives, min_score, max_score,
#'           positive_pct, negative_pct
#'   }
#'
#' @export
compute_synchronicity <- function(event_stats,
                                  event_names,
                                  confidence_level = 0.95,
                                  age_difference = 0.05,
                                  horizon_groups = NULL,
                                  n_samples = 10000,
                                  seed = 5128) {

  set.seed(seed)

  results           <- list()
  viz_data          <- data.frame()
  overall_scores    <- list()
  all_horizon_stats <- data.frame(
    Event = character(), record1 = character(), record2 = character(),
    variant1 = character(), variant2 = character(),
    mean_LR = numeric(), SD_LR = numeric(),
    confidence_level = numeric(), age_difference = numeric(),
    synchronicity_score = numeric(), synchronicity_precision = numeric(),
    stringsAsFactors = FALSE
  )

  register_pair <- function(pr) {
    if (is.null(results[[pr$lr_name]]))
      results[[pr$lr_name]] <<- list(ABlr = list(),
                                      synchronicity_percentage = list(),
                                      synchronicity_precision  = list())
    results[[pr$lr_name]]$ABlr[[pr$variant_key]]                    <<- pr$ABlr
    results[[pr$lr_name]]$synchronicity_percentage[[pr$variant_key]] <<- pr$score
    results[[pr$lr_name]]$synchronicity_precision[[pr$variant_key]]  <<- pr$precision
    all_horizon_stats <<- rbind(all_horizon_stats, pr$stats_row)
    viz_data          <<- rbind(viz_data,          pr$viz_row)
  }

  for (event_name in event_names) {

    unique_horizons <- unique(unlist(lapply(event_stats$processed, colnames)))

    if (!is.null(horizon_groups) && event_name %in% names(horizon_groups)) {
      event_base        <- event_name
      matching_horizons <- as.character(horizon_groups[[event_name]])
      has_generic       <- length(matching_horizons) > 1
    } else {
      event_base  <- event_name
      has_generic <- FALSE
    }

    if (has_generic) {
      #------------------------------------------------------------------------
      # GENERIC HORIZON GROUP
      #------------------------------------------------------------------------

      records_variants <- list()
      for (rec_name in names(event_stats$processed)) {
        present_variants <- intersect(matching_horizons, names(event_stats$processed[[rec_name]]))
        if (length(present_variants) > 0) {
          if (event_base %in% present_variants)
            present_variants <- c(event_base, setdiff(present_variants, event_base))
          records_variants[[rec_name]] <- present_variants
        }
      }

      relevant_record_names <- names(records_variants)
      if (length(relevant_record_names) == 0) next

      th <- resolve_horizon_thresholds(
        event_base, confidence_level, age_difference,
        event_stats$summaries[relevant_record_names]
      )

      for (i in 1:(length(relevant_record_names) - 1)) {
        for (j in (i + 1):length(relevant_record_names)) {
          rec_i <- relevant_record_names[i]
          rec_j <- relevant_record_names[j]

          for (var_i in records_variants[[rec_i]]) {
            for (var_j in records_variants[[rec_j]]) {
              col1 <- event_stats$processed[[rec_i]][[var_i]]
              col2 <- event_stats$processed[[rec_j]][[var_j]]
              if (is.null(col1) || is.null(col2)) next

              register_pair(compare_pair(col1, col2, rec_i, rec_j, var_i, var_j,
                                         horizon = event_base, thresholds = th,
                                         summaries = event_stats$summaries,
                                         n_samples = n_samples))
            }
          }
        }
      }

      # Overall synchronicity for this generic group
      if (length(records_variants) > 0) {

        for (i in 1:(length(relevant_record_names) - 1)) {
          for (j in (i + 1):length(relevant_record_names)) {
            rec_i <- relevant_record_names[i]
            rec_j <- relevant_record_names[j]

            for (var_i in records_variants[[rec_i]]) {
              for (var_j in records_variants[[rec_j]]) {
                pair_samples          <- list()
                pair_samples[[rec_i]] <- event_stats$processed[[rec_i]][[var_i]]
                pair_samples[[rec_j]] <- event_stats$processed[[rec_j]][[var_j]]

                if (!is.null(pair_samples[[rec_i]]) && !is.null(pair_samples[[rec_j]])) {
                  pair_result <- calculate_overall_synchronicity(
                    pair_samples, th$conf_level_h, th$age_diff_log_bounds,
                    n_samples = n_samples, seed = seed
                  )
                  pair_key <- paste(rec_i, rec_j, var_i, var_j, sep = "_")
                  if (is.null(overall_scores[[event_base]]))
                    overall_scores[[event_base]] <- list()

                  overall_scores[[event_base]][[pair_key]] <- list(
                    variant_combination = paste(rec_i, var_i, rec_j, var_j, sep = ":"),
                    overall_score       = pair_result$overall_score,
                    overall_precision   = pair_result$overall_precision,
                    n_records           = 2,
                    ref_record          = pair_result$ref_record
                  )
                }
              }
            }
          }
        }

        if (length(records_variants) > 2) {
          variant_lists     <- lapply(records_variants, function(v) v)
          all_combos        <- expand.grid(variant_lists, stringsAsFactors = FALSE)
          names(all_combos) <- names(records_variants)

          for (combo_idx in 1:nrow(all_combos)) {
            combo        <- as.character(all_combos[combo_idx, ])
            names(combo) <- names(records_variants)

            combo_samples_list <- list()
            combo_valid <- TRUE
            for (rec in names(combo)) {
              df <- event_stats$processed[[rec]]
              if (combo[[rec]] %in% names(df)) {
                combo_samples_list[[rec]] <- df[[combo[[rec]]]]
              } else {
                combo_valid <- FALSE; break
              }
            }

            if (combo_valid && length(combo_samples_list) > 2) {
              combo_result <- calculate_overall_synchronicity(
                combo_samples_list, th$conf_level_h, th$age_diff_log_bounds,
                n_samples = n_samples, seed = seed
              )
              combo_key <- paste(combo, collapse = "_")
              overall_scores[[event_base]][[combo_key]] <- list(
                variant_combination = paste(names(combo), combo, sep = ":", collapse = " | "),
                overall_score       = combo_result$overall_score,
                overall_precision   = combo_result$overall_precision,
                n_records           = combo_result$n_records,
                ref_record          = combo_result$ref_record
              )
              all_horizon_stats <- rbind(all_horizon_stats, data.frame(
                Event            = event_base,
                record1          = "MULTI",
                record2          = "MULTI",
                variant1         = paste(names(combo), combo, sep = ":", collapse = " | "),
                variant2         = NA,
                mean_LR          = NA,
                SD_LR            = NA,
                confidence_level = th$conf_level_h,
                age_difference   = th$age_diff_h,
                synchronicity_score     = combo_result$overall_score,
                synchronicity_precision = combo_result$overall_precision,
                stringsAsFactors = FALSE
              ))
            }
          }
        }
      }

    } else {
      #------------------------------------------------------------------------
      # SINGLE HORIZON
      #------------------------------------------------------------------------

      matching_horizons <- unique_horizons[unique_horizons == event_name]
      if (length(matching_horizons) == 0) next

      for (horizon in matching_horizons) {

        relevant_records <- event_stats$processed[
          vapply(event_stats$processed, function(df) horizon %in% names(df), logical(1))
        ]
        relevant_record_names <- names(relevant_records)

        th <- resolve_horizon_thresholds(
          horizon, confidence_level, age_difference,
          event_stats$summaries[relevant_record_names]
        )

        if (length(relevant_records) > 1) {
          for (i in 1:(length(relevant_records) - 1)) {
            for (j in (i + 1):length(relevant_records)) {
              col1 <- relevant_records[[i]][[horizon]]
              col2 <- relevant_records[[j]][[horizon]]
              if (is.null(col1) || is.null(col2)) next

              register_pair(compare_pair(col1, col2,
                                         relevant_record_names[i], relevant_record_names[j],
                                         horizon = horizon, thresholds = th,
                                         summaries = event_stats$summaries,
                                         n_samples = n_samples))
            }
          }
        }

        overall_samples_list <- list()
        for (rec in relevant_record_names) {
          df <- event_stats$processed[[rec]]
          if (horizon %in% names(df)) overall_samples_list[[rec]] <- df[[horizon]]
        }

        if (length(overall_samples_list) >= 2) {
          overall_result <- calculate_overall_synchronicity(
            overall_samples_list, th$conf_level_h, th$age_diff_log_bounds,
            n_samples = n_samples, seed = seed
          )
          overall_scores[[horizon]] <- list(
            all = list(
              variant_combination = "single",
              overall_score       = overall_result$overall_score,
              overall_precision   = overall_result$overall_precision,
              n_records           = overall_result$n_records,
              ref_record          = overall_result$ref_record
            )
          )
          all_horizon_stats <- rbind(all_horizon_stats, data.frame(
            Event            = horizon,
            record1          = "OVERALL",
            record2          = "OVERALL",
            variant1         = NA,
            variant2         = NA,
            mean_LR          = NA,
            SD_LR            = NA,
            confidence_level = th$conf_level_h,
            age_difference   = th$age_diff_h,
            synchronicity_score     = overall_result$overall_score,
            synchronicity_precision = overall_result$overall_precision,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }

  all_scores <- data.frame(record_pair = character(), variant_pair = character(),
                            score = numeric(), difference = numeric(), stringsAsFactors = FALSE)
  for (lr_name in names(results)) {
    for (variant_key in names(results[[lr_name]]$synchronicity_percentage)) {
      all_scores <- rbind(all_scores, data.frame(
        record_pair  = lr_name,
        variant_pair = variant_key,
        score        = results[[lr_name]]$synchronicity_percentage[[variant_key]],
        difference   = results[[lr_name]]$synchronicity_precision[[variant_key]],
        stringsAsFactors = FALSE
      ))
    }
  }

  summary_out <- list(total = 0L, positives = 0L, negatives = 0L,
                       min_score = NA_real_, max_score = NA_real_,
                       positive_pct = NA_real_, negative_pct = NA_real_)
  if (nrow(all_scores) > 0) {
    total     <- nrow(all_scores)
    positives <- sum(all_scores$score >= confidence_level[1], na.rm = TRUE)
    negatives <- total - positives
    summary_out <- list(
      total        = total,    positives    = positives,  negatives    = negatives,
      min_score    = min(all_scores$score, na.rm = TRUE),
      max_score    = max(all_scores$score, na.rm = TRUE),
      positive_pct = (positives / total) * 100,
      negative_pct = (negatives / total) * 100
    )
  }

  return(list(results = results, overall_scores = overall_scores, viz_data = viz_data,
               all_scores = all_scores, all_horizon_stats = all_horizon_stats,
               summary = summary_out))
}

#' Write and Plot Synchronicity Results
#'
#' I/O counterpart to \code{compute_synchronicity()}: saves a CSV of horizon statistics,
#' prints a summary to the console, and produces per-horizon PDF visualizations.
#' Call \code{compute_synchronicity()} first and pass its result here.
#'
#' @param synchro_result The list returned by \code{compute_synchronicity()}.
#' @param event_names Character vector of event names being tested. Falls back to
#'   naming the CSV file when \code{group_name} is not supplied.
#' @param output_dir Character string specifying the directory the statistics CSV will be
#'   saved to (default: \code{syncer_output_dir()}, i.e. the \code{SyncER_outputs}
#'   folder in the working directory).
#' @param group_name Character string used to name the CSV file (as
#'   \code{"<group_name>_stats.csv"}), representing the horizon group being tested
#'   as a whole rather than every individual horizon name in \code{event_names}. When
#'   \code{NULL} (default), defaults to \code{"isochron"} if \code{isochron = TRUE}, or
#'   to all of \code{event_names} pasted together otherwise.
#' @param synced Character string appended to output filenames (default: \code{""}).
#' @param isochron Logical; when \code{TRUE} a warning is shown if >10\% of comparisons
#'   fail, and the CSV file is named \code{"isochron"} by default (default: \code{FALSE}).
#' @param offset Numeric offset correction value passed to \code{plot_synchronicity()}
#'   (default: \code{bp_datum()}, i.e. the current year minus 1950).
#' @param fig_width Numeric width of output PDF figures in inches (default: 10).
#' @param fig_height Numeric height of output PDF figures in inches (default: 8).
#' @inheritParams create_visualization
#'
#' @return Invisibly returns \code{synchro_result} unchanged.
#'
#' @export
verify_synchronicity <- function(synchro_result,
                                 event_names,
                                 output_dir = syncer_output_dir(),
                                 group_name = NULL,
                                 synced = "",
                                 isochron = FALSE,
                                 offset = bp_datum(),
                                 fig_width = 10,
                                 fig_height = 8,
                                 plot_opts = list()) {

  synchro <- synchro_result

  if (is.null(group_name)) {
    group_name <- if (isochron) "isochron" else paste(event_names, collapse = "_")
  }

  # Save CSV statistics
  if (nrow(synchro$all_horizon_stats) > 0) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    stats_file <- file.path(output_dir,
                            paste0(group_name, "_stats", synced, ".csv"))
    write.csv(synchro$all_horizon_stats, stats_file, row.names = FALSE)
    cat(sprintf("\nAll horizon statistics saved to: %s\n", stats_file))
  }

  # Print summary to console
  s <- synchro$summary
  if (s$total > 0) {
    cat("\n=== Overall Synchronicity Test Summary ===\n")
    if (isochron && s$negative_pct > 10) {
      cat(sprintf(
        "\n  You have a high percentage (%.2f%%) of ages that do not pass the test.\n   Please verify that your isochrons are correctly labeled!\n\n",
        s$negative_pct
      ))
    }
    cat("Total comparisons:", s$total, "\n")
    cat("Pass:", s$positives, sprintf("(%.1f%%)", s$positive_pct),
        "| Fail:", s$negatives, sprintf("(%.1f%%)\n", s$negative_pct))
    cat("Lowest score:",  sprintf("%.1f%%", s$min_score * 100),
        "| Highest score:", sprintf("%.1f%%\n", s$max_score * 100))
    cat("==========================================\n\n")
  }

  # Generate PDFs and console plots
  plot_synchronicity(
    synchro_result = synchro,
    output_dir     = output_dir,
    offset         = offset,
    synced         = synced,
    fig_width      = fig_width,
    fig_height     = fig_height,
    plot_opts      = plot_opts
  )

  invisible(synchro_result)
}

#' Calculate Validation Thresholds from Synchronized Data
#'
#' Calculates empirically-based validation thresholds and confidence levels from synchronized
#' age data for use in subsequent synchronicity testing.
#'
#' @param adjusted_ages Output from \code{set_to_zero()} containing synchronized ages and errors.
#' @param sigma_multiplier Numeric value or named vector giving the user-defined confidence
#'   level on the age interval, expressed as a standard-deviation multiplier (default: \code{2},
#'   i.e. ~95.4\% confidence on the age; use \code{1} for ~68.3\% or \code{3} for ~99.7\%). This
#'   sets how wide an age difference each isochron's assigned error implies, and therefore what
#'   \code{age_difference} values are reasonable when testing nearby unknown deposits. Supply a
#'   named vector to use a different multiplier per horizon.
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
