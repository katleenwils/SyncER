#' Build Bacon Input Data Structures
#'
#' Prepares radiocarbon and Pb210 data frames in the format required by the Bacon
#' age-depth modeling software, without performing any file I/O.
#'
#' @param record_data A data frame containing record data with at least columns: \code{event},
#'   \code{depth}, \code{C14_age}, \code{C14_error}, and \code{cc} (calibration curve indicator).
#' @param record_name Character string used to identify the record.
#' @param radiocarbon_sample_names Character vector of the event label(s) given to the
#'   radiocarbon samples you want to include (default: \code{"sample"}).
#' @param lead_sample_names Character vector of the event label(s) given to the Pb-210 (lead)
#'   samples you want to include (default: \code{"210Pb_sample"}); leave blank if you only
#'   consider radiocarbon ages.
#' @param adjusted_ages Optional list output from \code{set_to_zero()} containing adjusted ages,
#'   errors, and depths per horizon (default: NULL).
#' @param original_ages Logical or named logical vector indicating whether to include original
#'   radiocarbon ages when \code{adjusted_ages} are provided (default: TRUE).
#'
#' @return A named list with elements:
#'   \itemize{
#'     \item \code{c14}: Data frame formatted for Bacon C14 input
#'     \item \code{pb210}: Data frame formatted for rplum Pb210 input, or \code{NULL} if no
#'           Pb210 data is present
#'   }
#'
#' @importFrom magrittr %>%
#' @export
build_bacon_input <- function(record_data,
                              record_name,
                              radiocarbon_sample_names = "sample",
                              lead_sample_names = "210Pb_sample",
                              adjusted_ages = NULL,
                              original_ages = TRUE) {

  # ── Check if Pb210 data is present ────────────────────────────────────────
  has_lead <- !is.null(lead_sample_names) &&
    length(lead_sample_names) > 0 &&
    any(record_data$event %in% lead_sample_names)

  # ── 1. RADIOCARBON data frame ──────────────────────────────────────────────
  sample_events <- record_data %>%
    dplyr::filter(event %in% radiocarbon_sample_names) %>%
    dplyr::select(depth, C14_age, C14_error, cc,
                  dplyr::any_of(c("t.a", "t.b")))  # include if present

  has_ta_tb <- all(c("t.a", "t.b") %in% colnames(sample_events))

  if (has_ta_tb) {
    colnames(sample_events)[1:4] <- c("depth", "age", "error", "cc")
  } else {
    colnames(sample_events) <- c("depth", "age", "error", "cc")
  }

  sample_events$age   <- as.numeric(sample_events$age)
  sample_events$error <- as.numeric(sample_events$error)
  sample_events$depth <- as.numeric(sample_events$depth)
  sample_events$cc    <- as.numeric(sample_events$cc)
  sample_events$labID  <- seq_len(nrow(sample_events))
  sample_events$d.R    <- 0
  sample_events$d.STD  <- 0

  if (!is.null(adjusted_ages)) {
    if (!has_ta_tb) {
      sample_events$t.a <- 3
      sample_events$t.b <- 4
    } else {
      # Source data had t.a/t.b columns but may contain NAs.
      # These are radiocarbon dates, so fill with the radiocarbon Student-t defaults (3, 4),
      # not the synchronized horizon values (33, 34) used in new_rows below.
      sample_events$t.a[is.na(sample_events$t.a)] <- 3
      sample_events$t.b[is.na(sample_events$t.b)] <- 4
    }
    sample_events <- sample_events %>%
      dplyr::select(labID, age, error, depth, cc, d.R, d.STD, t.a, t.b)
  } else {
    if (has_ta_tb) {
      sample_events <- sample_events %>%
        dplyr::select(labID, age, error, depth, cc, d.R, d.STD, t.a, t.b)
    } else {
      sample_events <- sample_events %>%
        dplyr::select(labID, age, error, depth, cc, d.R, d.STD)
    }
  }

  # Resolve original_ages for this record: named vector → look up by record_name
  use_original_ages <- if (!is.null(names(original_ages)) && record_name %in% names(original_ages)) {
    isTRUE(original_ages[[record_name]])
  } else {
    isTRUE(original_ages[1])
  }

  if (!is.null(adjusted_ages) && !use_original_ages) {
    sample_events <- sample_events %>% dplyr::filter(depth == 0)
  }

  # Add adjusted/synchronised horizon ages if provided
  if (!is.null(adjusted_ages)) {
    new_rows <- list()
    labID_counter <- if (nrow(sample_events) > 0) max(sample_events$labID) + 1 else 1

    for (horizon_name in names(adjusted_ages)) {
      df <- adjusted_ages[[horizon_name]]
      if (record_name %in% df$record) {
        i <- which(df$record == record_name)
        match_row <- record_data %>% dplyr::filter(event == horizon_name)

        if (nrow(match_row) == 0) {
          warning(paste("No matching depth found for horizon:", horizon_name,
                        "in record", record_name))
          next
        }

        new_rows[[length(new_rows) + 1]] <- data.frame(
          labID     = labID_counter,
          age       = df$adjusted_age[i],
          error     = df$adjusted_error[i],
          depth     = match_row$depth[1],
          cc        = if ("cc" %in% names(df)) df$cc[i] else 0,
          d.R       = 0,
          d.STD     = 0,
          t.a       = 33,
          t.b       = 34,
          stringsAsFactors = FALSE
        )
        labID_counter <- labID_counter + 1
      }
    }

    if (length(new_rows) > 0) {
      new_rows_df <- do.call(rbind, new_rows)
      sample_events <- rbind(as.data.frame(sample_events), new_rows_df)
    }
  }

  sample_events <- sample_events[order(sample_events$depth), ]

  # ── 2. Pb210 data frame (rplum / HP1C format) ─────────────────────────────
  lead_df <- NULL
  if (has_lead) {

    lead_events <- record_data %>%
      dplyr::filter(event %in% lead_sample_names)

    has_ra <- all(c("Ra226", "Ra226_sd") %in% colnames(lead_events)) &&
      any(!is.na(lead_events$Ra226))

    if (has_ra) {
      lead_events <- lead_events %>%
        dplyr::select(depth, density, Pb210, Pb210_sd, thickness, Ra226, Ra226_sd, settings)

      colnames(lead_events) <- c("depth(cm)", "density(g/cm^3)", "210Pb(Bq/kg)",
                                 "sd(210Pb)", "thickness(cm)", "226Ra(Bq/kg)", "sd(226Ra)", "settings")

      lead_events$labID <- seq_len(nrow(lead_events))

      lead_events <- lead_events %>%
        dplyr::select(labID, `depth(cm)`, `density(g/cm^3)`, `210Pb(Bq/kg)`,
                      `sd(210Pb)`, `thickness(cm)`, `226Ra(Bq/kg)`, `sd(226Ra)`, settings)

    } else {
      lead_events <- lead_events %>%
        dplyr::select(depth, density, Pb210, Pb210_sd, thickness, settings)

      colnames(lead_events) <- c("depth(cm)", "density(g/cm^3)", "210Pb(Bq/kg)",
                                 "sd(210Pb)", "thickness(cm)", "settings")

      lead_events$labID <- seq_len(nrow(lead_events))

      lead_events <- lead_events %>%
        dplyr::select(labID, `depth(cm)`, `density(g/cm^3)`, `210Pb(Bq/kg)`,
                      `sd(210Pb)`, `thickness(cm)`, settings)
    }

    lead_df <- lead_events[order(lead_events$`depth(cm)`), ]
  }

  return(list(c14 = sample_events, pb210 = lead_df))
}

# ── Internal helper: bake synchronized ages into one record's frame ──────────
# For each horizon in `adjusted_ages`, overwrite the matching event row's age /
# error / cc in `frame` with the synchronized value. Ages carry the generation
# they were computed from in their `record` field (e.g. "core1_synced"), while
# `record_key` may be a later generation (e.g. "core1_synced_1"); they are matched
# when `record_key` equals, or is a child (\code{<record>_...}) of, an age's
# `record`. Records/horizons not found are left untouched.
apply_adjusted_ages_to_frame <- function(frame, record_key, adjusted_ages) {

  for (horizon_name in names(adjusted_ages)) {
    adj <- adjusted_ages[[horizon_name]]

    match_age <- vapply(adj$record, function(r) {
      identical(record_key, r) || startsWith(record_key, paste0(r, "_"))
    }, logical(1))
    cand <- which(match_age)
    if (length(cand) == 0) next
    i <- cand[which.max(nchar(adj$record[cand]))]  # most specific (longest) match

    rows <- which(frame$event == horizon_name)
    if (length(rows) == 0) {
      rows <- which(grepl(paste0("^", horizon_name, "[a-zA-Z0-9]*$"), frame$event))
    }
    if (length(rows) == 0) next

    if ("C14_age"   %in% names(frame)) frame$C14_age[rows]   <- adj$adjusted_age[i]
    if ("C14_error" %in% names(frame)) frame$C14_error[rows] <- adj$adjusted_error[i]
    if ("cc" %in% names(frame))        frame$cc[rows] <- if ("cc" %in% names(adj)) adj$cc[i] else 0
  }

  frame
}

# ── Internal helper: reconstruct the synchronized-age list from a frame ──────
# Any dated event (non-NA C14_age) that is neither a radiocarbon sample nor a
# Pb210 sample is treated as a synchronized horizon, and returned in the
# \code{adjusted_ages} shape build_bacon_input() consumes.
extract_synced_ages <- function(frame, record_key,
                                radiocarbon_sample_names, lead_sample_names) {

  if (!"C14_age" %in% names(frame)) return(list())

  is_dated <- !is.na(frame$C14_age)
  is_rc    <- frame$event %in% radiocarbon_sample_names
  is_lead  <- frame$event %in% lead_sample_names
  idx      <- which(is_dated & !is_rc & !is_lead)

  out <- list()
  for (k in idx) {
    out[[frame$event[k]]] <- data.frame(
      record         = record_key,
      adjusted_age   = frame$C14_age[k],
      adjusted_error = frame$C14_error[k],
      cc             = if ("cc" %in% names(frame)) frame$cc[k] else 0,
      stringsAsFactors = FALSE
    )
  }
  out
}

# ── Internal helper: name the next folder generation for a record ────────────
# The first synchronization of a base record appends "_synced"; every later
# generation increments a trailing number on that base ("_synced_1", "_synced_2",
# ...), choosing the next unused number among sibling folders in `base_dir`.
next_generation_name <- function(record_key, base_dir) {

  if (!grepl("_synced", record_key)) {
    return(paste0(record_key, "_synced"))
  }

  base     <- sub("_[0-9]+$", "", record_key)
  siblings <- basename(list.dirs(base_dir, recursive = FALSE))
  nums     <- suppressWarnings(as.integer(
    sub(paste0("^", base, "_"), "",
        siblings[grepl(paste0("^", base, "_[0-9]+$"), siblings)])
  ))
  n <- if (length(nums) == 0) 1L else max(nums, na.rm = TRUE) + 1L
  paste0(base, "_", n)
}

# ── Internal helper: write a build_bacon_input() result to a folder ──────────
write_bacon_csvs <- function(bacon_input, out_key, base_dir) {

  folder <- file.path(base_dir, out_key)
  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)

  has_lead <- !is.null(bacon_input$pb210)
  c14_path <- file.path(folder, if (has_lead) paste0(out_key, "_C14.csv") else paste0(out_key, ".csv"))
  write.csv(bacon_input$c14, c14_path, row.names = FALSE, quote = FALSE, na = "")
  message("C14 data for ", out_key, " saved to: ", c14_path)

  if (has_lead) {
    pb_path <- file.path(folder, paste0(out_key, ".csv"))
    write.csv(bacon_input$pb210, pb_path, row.names = FALSE, quote = FALSE, na = "")
    message("Pb210 data for ", out_key, " saved to: ", pb_path)
  }
}

#' Write Bacon Input Files and Return the Updated Record Metadata
#'
#' Single entry point for turning record metadata into Bacon-ready CSV input files. It
#' treats \code{record_data} as the source of truth: any new synchronized ages are baked
#' into it, the Bacon CSV(s) are rebuilt from the radiocarbon samples plus \emph{all}
#' synchronized ages now present, and the updated \code{record_data} (re-keyed to the folder
#' generation that was written) is returned -- ready to pass straight into the next call.
#'
#' @param record_data Named list of record metadata (output from \code{load_excel_data()},
#'   or the return value of a previous \code{age_model_input()} call). Each element is a data
#'   frame with at least \code{event}, \code{depth}, \code{C14_age}, \code{C14_error} and
#'   \code{cc} columns.
#' @param adjusted_ages Optional output from \code{set_to_zero()} (or a \code{c()} combination
#'   of several) with the new synchronized ages to add (default: \code{NULL}). When
#'   \code{NULL}, the current \code{record_data} is written as-is (the initial Bacon input).
#' @param update_records Logical; when \code{TRUE} the existing folders that \code{record_data}
#'   is keyed to are overwritten in place. When \code{FALSE} (default) and \code{adjusted_ages}
#'   is supplied, a new folder generation is created (\code{"_synced"} for the first
#'   synchronization, then \code{"_synced_1"}, \code{"_synced_2"}, ...).
#' @param radiocarbon_sample_names Character vector of the event label(s) used for your
#'   radiocarbon samples in the input file (default: \code{"sample"}).
#' @param lead_sample_names Character vector of the event label(s) used for your Pb-210 (lead)
#'   samples (default: \code{""}); leave blank if you only consider radiocarbon ages.
#' @param original_ages Logical or named logical vector indicating whether the new age-depth
#'   models should also use the original (radiocarbon) ages alongside the synchronized ages
#'   (default: \code{TRUE}). Set to \code{FALSE} to build the models from the synchronized ages
#'   only (e.g. a record with poor accuracy); adjust it per record with a named vector such as
#'   \code{c("core1" = FALSE)} -- records not listed fall back to \code{TRUE}. Passed through to
#'   \code{build_bacon_input()}.
#' @param base_dir Character string specifying the output directory (default: \code{"."}, the
#'   working directory, matching the \code{coredir} used by \emph{rbacon}/\emph{rplum}).
#'
#' @return Invisibly, a named list with the same structure as \code{record_data}: the record
#'   frames with the new synchronized ages baked in, keyed by the folder generation that was
#'   written (unchanged keys when \code{update_records = TRUE}). Assign it and pass it to the
#'   next \code{age_model_input()} call.
#'
#' @export
age_model_input <- function(record_data,
                            adjusted_ages = NULL,
                            update_records = FALSE,
                            radiocarbon_sample_names = "sample",
                            lead_sample_names = "",
                            original_ages = TRUE,
                            base_dir = ".") {

  result <- list()

  for (key in names(record_data)) {
    frame <- record_data[[key]]

    # 1. Bake any new synchronized ages into this record's frame.
    if (!is.null(adjusted_ages)) {
      frame <- apply_adjusted_ages_to_frame(frame, key, adjusted_ages)
    }

    # 2. Reconstruct every synchronized age now present in the frame.
    full_adj <- extract_synced_ages(frame, key, radiocarbon_sample_names, lead_sample_names)

    # 3. Rebuild the Bacon-ready input from radiocarbon samples + all synchronized ages.
    bacon_input <- build_bacon_input(
      frame, key,
      radiocarbon_sample_names = radiocarbon_sample_names,
      lead_sample_names        = lead_sample_names,
      adjusted_ages            = if (length(full_adj) > 0) full_adj else NULL,
      original_ages            = original_ages
    )

    # 4. Resolve the output folder: in situ, base (initial build), or new generation.
    out_key <- if (isTRUE(update_records)) {
      key
    } else if (is.null(adjusted_ages)) {
      key
    } else {
      next_generation_name(key, base_dir)
    }

    # 5. Write the CSV(s).
    write_bacon_csvs(bacon_input, out_key, base_dir)

    # 6. Collect the updated frame, keyed by the folder that was written.
    result[[out_key]] <- frame
  }

  invisible(result)
}

#' Process Event Ages
#'
#' Extracts event-specific age columns from Bacon output, applies offset correction,
#' and calculates summary statistics.
#'
#' @param out_data Named list of data frames or data.tables containing age data
#'   (output from \code{extract_ages()}).
#' @param event_deposits Character vector of event deposit names to extract (typically
#'   \code{isochrons} or \code{test_events}). Columns whose name starts with any of these names
#'   are selected.
#' @param offset Numeric value added to every age so that the age reference point becomes, for
#'   example, the year the (most recent) records were retrieved rather than 1950 AD (default:
#'   \code{bp_datum()}, i.e. the current year minus 1950). Because the synchronicity test uses
#'   log-transformations, the age dataset may contain no zero or negative values; this offset
#'   guarantees that for post-1950 (modern) deposits. Set it manually if you want a different
#'   reference point, but use the \emph{same} offset consistently throughout the workflow.
#'
#' @return A named list with two elements:
#'   \itemize{
#'     \item \code{processed}: List of data frames with offset-corrected ages for each record
#'     \item \code{summaries}: List of summary data frames containing min, max, mean, and
#'           sd (sigma) for each event column
#'   }
#'
#' @details Filters columns matching event deposit patterns, applies offset to avoid negative
#'   ages for post-1950 deposits, and ensures the record top is at age 0. This function takes
#'   the direct output from \code{extract_ages()} and processes it in memory without requiring
#'   intermediate file I/O.
#'
#' @importFrom magrittr %>%
#' @export
process_event_ages <- function(out_data,
                               event_deposits,
                               offset = bp_datum()) {
  
  # Validate input
  if (!is.list(out_data) || length(out_data) == 0) {
    stop("out_data must be a non-empty list (output from extract_ages)")
  }
  
  # Create regex pattern to match event deposit columns
  # Example: if event_deposits = c("tephra", "flood"), matches "tephra1", "flood_120", etc.
  pattern_regex <- paste0("^(", paste(event_deposits, collapse = "|"), ")")
  
  # Process each record: select matching columns and apply offset
  processed <- lapply(out_data, function(.x) {
    # Convert to data frame if not already
    if (!is.data.frame(.x)) {
      .x <- as.data.frame(.x)
    }
    
    # Select columns matching event deposits
    matched_cols <- grep(pattern_regex, names(.x), value = TRUE)
    
    if (length(matched_cols) == 0) {
      warning("No matching event columns found in this record")
      return(data.frame())
    }
    
    # Select and apply offset
    .x[, matched_cols, drop = FALSE] %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ . + offset))  # Offset adjusts for modern ages
  })
  
  # Calculate summary statistics for each column in each record
  summaries <- lapply(processed, function(.x) {
    if (ncol(.x) == 0) {
      return(data.frame())
    }
    
    dplyr::summarise(.x,
                     dplyr::across(dplyr::everything(),
                                   list(
                                     min = ~ min(.x, na.rm = TRUE),
                                     max = ~ max(.x, na.rm = TRUE),
                                     mean = ~ mean(.x, na.rm = TRUE),
                                     sigma = ~ sd(.x, na.rm = TRUE)
                                   )
                     )
    )
  })
  
  return(list(processed = processed, summaries = summaries))
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
