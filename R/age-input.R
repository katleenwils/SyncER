#' Build input data structures compatible with Bacon and Plum age-depth modelling 
#'
#' Prepares radiocarbon and Pb210 data frames in the format required by the \emph{rbacon}/\emph{rplum} 
#' age-depth modelling software.
#'
#' @param record_data A data frame containing record data with at least columns: \code{event},
#'   \code{depth}, \code{C14_age}, \code{C14_error}, and \code{cc} (calibration curve indicator).
#' @param record_name Character string used to identify the record.
#' @param radiocarbon_sample_names Character vector of the event label(s) given to the
#'   radiocarbon samples you want to include (default: \code{"sample"}).
#' @param lead_sample_names Character vector of the event label(s) given to the Pb-210 (lead)
#'   samples you want to include (default: \code{"210Pb_sample"}); leave blank if you only
#'   consider radiocarbon ages.
#' @param adjusted_ages Optional list output from \code{synchronize_ages()} containing adjusted ages,
#'   errors, and depths per horizon (default: NULL).
#' @param original_ages Logical or named logical vector indicating whether to include original
#'   radiocarbon ages when \code{adjusted_ages} are provided (default: TRUE).
#'
#' @return A named list with elements:
#'   \itemize{
#'     \item \code{c14}: Data frame formatted for \emph{rbacon} C14 input
#'     \item \code{pb210}: Data frame formatted for \emph{rplum} containing Pb210 and C14 input, or \code{NULL} if no
#'           Pb210 data is present
#'   }
#'
#' @importFrom magrittr %>%
#' @export
build_age_input <- function(record_data,
                              record_name,
                              radiocarbon_sample_names = "sample",
                              lead_sample_names = "210Pb_sample",
                              adjusted_ages = NULL,
                              original_ages = TRUE) {

  # -- Check if Pb210 data is present ----------------------------------------
  has_lead <- !is.null(lead_sample_names) &&
    length(lead_sample_names) > 0 &&
    any(record_data$event %in% lead_sample_names)

  # -- 1. RADIOCARBON data frame ----------------------------------------------
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

  # Resolve original_ages for this record: named vector -> look up by record_name
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

  # -- 2. Pb210 data frame (rplum / HP1C format) -----------------------------
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

#' Add synchronized ages to age data
#'
#' For each horizon in \code{adjusted_ages}, overwrites the matching event row's
#' age / error / cc in \code{frame} with the synchronized value. Ages carry the
#' generation they were computed from in their \code{record} field (e.g.
#' \code{"core1_synced"}), while \code{record_key} may be a later generation (e.g.
#' \code{"core1_synced_1"}); they are matched when \code{record_key} equals, or is
#' a child (\code{<record>_...}) of, an age's \code{record}. Records/horizons not
#' found are left untouched.
#'
#' @param frame Data frame for a single record (\emph{rbacon}/\emph{rplum}-style, with an \code{event}
#'   column plus \code{C14_age}, \code{C14_error} and optionally \code{cc}).
#' @param record_key Character key identifying the record generation \code{frame}
#'   belongs to.
#' @param adjusted_ages Named list of synchronized-age data frames, as returned by
#'   \code{synchronize_ages()}.
#'
#' @return \code{frame} with matching event rows overwritten by the synchronized
#'   ages.
#'
#' @keywords internal
add_adjusted_ages <- function(frame, record_key, adjusted_ages) {

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

#' Reconstruct the synchronized age list from a frame
#'
#' Any dated event (non-NA \code{C14_age}) that is neither a radiocarbon sample nor
#' a Pb210 sample is treated as a synchronized horizon, and returned in the
#' \code{adjusted_ages} shape \code{build_age_input()} consumes.
#'
#' @param frame Data frame for a single record, with an \code{event} column and
#'   \code{C14_age} / \code{C14_error} (and optionally \code{cc}).
#' @param record_key Character key identifying the record the ages belong to.
#' @param radiocarbon_sample_names Character vector of event label(s) marking
#'   radiocarbon samples.
#' @param lead_sample_names Character vector of event label(s) marking Pb-210
#'   samples.
#'
#' @return Named list of one-row data frames (one per synchronized horizon), or an
#'   empty list when \code{frame} has no \code{C14_age} column.
#'
#' @keywords internal
read_adjusted_ages <- function(frame, record_key,
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

#' Determine the next folder generation (after syncing) for a record
#'
#' The first synchronization of a base record appends \code{"_synced"}; every later
#' generation increments a trailing number on that base (\code{"_synced_1"},
#' \code{"_synced_2"}, ...), choosing the next unused number among sibling folders
#' in \code{base_dir}.
#'
#' @param record_key Character key of the record being written.
#' @param base_dir Character path whose sub-folders are scanned for existing
#'   generations.
#'
#' @return Character folder name for the next generation.
#'
#' @keywords internal
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

#' Write a build_age_input() result to a folder
#'
#' Writes the radiocarbon CSV (and, when present, the Pb210 CSV) of a
#' \code{build_age_input()} result into \code{base_dir/out_key}, creating the
#' folder if needed. Follows the \emph{rbacon}/\emph{rplum} naming convention:
#' \code{<key>_C14.csv} plus \code{<key>.csv} when lead data is present, otherwise
#' \code{<key>.csv} for the radiocarbon data.
#'
#' @param age_input List with a \code{c14} data frame and optional \code{pb210}
#'   data frame, as returned by \code{build_age_input()}.
#' @param out_key Character key naming the output folder and file stem.
#' @param base_dir Character path the \code{out_key} folder is created under.
#'
#' @return Invisibly \code{NULL}; called for its side effect of writing CSV files.
#'
#' @keywords internal
write_age_input_data <- function(age_input, out_key, base_dir) {

  folder <- file.path(base_dir, out_key)
  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)

  has_lead <- !is.null(age_input$pb210)
  c14_path <- file.path(folder, if (has_lead) paste0(out_key, "_C14.csv") else paste0(out_key, ".csv"))
  write.csv(age_input$c14, c14_path, row.names = FALSE, quote = FALSE, na = "")
  message("C14 data for ", out_key, " saved to: ", c14_path)

  if (has_lead) {
    pb_path <- file.path(folder, paste0(out_key, ".csv"))
    write.csv(age_input$pb210, pb_path, row.names = FALSE, quote = FALSE, na = "")
    message("Pb210 data for ", out_key, " saved to: ", pb_path)
  }
}

#' Write input files for age-depth modelling and return the updated record metadata
#'
#' Single entry point for turning record metadata into \emph{rbacon}/\emph{rplum}-ready CSV input files. It
#' treats \code{record_data} as the initial input data: any new synchronized ages are included and 
#' the CSV(s) are rebuilt. The updated \code{record_data} (re-keyed to the folder
#' generation that was written) is returned.
#'
#' @param record_data Named list of record metadata (output from \code{read_record_data()},
#'   or the return value of a previous \code{age_model_input()} call). Each element is a data
#'   frame with at least \code{event}, \code{depth}, \code{C14_age}, \code{C14_error} and
#'   \code{cc} columns.
#' @param adjusted_ages Optional output from \code{synchronize_ages()} (or a \code{c()} combination
#'   of several) with the new synchronized ages to add (default: \code{NULL}). When
#'   \code{NULL}, the current \code{record_data} is written as-is (the initial input).
#' @param update_records Logical; when \code{TRUE} the existing folders that \code{record_data}
#'   is keyed to are overwritten in place. When \code{FALSE} (default) and \code{adjusted_ages}
#'   is supplied, a new folder generation is created (\code{"_synced"} for the first
#'   synchronization, then \code{"_synced_1"}, \code{"_synced_2"}, ...).
#' @param radiocarbon_sample_names Character vector of the event label(s) used for your
#'   radiocarbon samples in the input file (default: \code{"sample"}).
#' @param lead_sample_names Character vector of the event label(s) used for your Pb-210 (lead)
#'   samples (default: \code{"210Pb_sample"}); leave blank if you only consider radiocarbon ages.
#' @param original_ages Logical or named logical vector indicating whether the new age-depth
#'   models should also use the original ages alongside the synchronized ages
#'   (default: \code{TRUE}). Set to \code{FALSE} to build the models from the synchronized ages
#'   only (e.g. a record with poor accuracy); adjust it per record with a named vector such as
#'   \code{c("core1" = FALSE)}. Records not listed fall back to \code{TRUE}. Passed through to
#'   \code{build_age_input()}.
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
                            lead_sample_names = "210Pb_sample",
                            original_ages = TRUE,
                            base_dir = ".") {

  result <- list()

  for (key in names(record_data)) {
    frame <- record_data[[key]]

    # 1. Bake any new synchronized ages into this record's frame.
    if (!is.null(adjusted_ages)) {
      frame <- add_adjusted_ages(frame, key, adjusted_ages)
    }

    # 2. Reconstruct every synchronized age now present in the frame.
    full_adj <- read_adjusted_ages(frame, key, radiocarbon_sample_names, lead_sample_names)

    # 3. Rebuild the Bacon-ready input from radiocarbon samples + all synchronized ages.
    age_input <- build_age_input(
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
    write_age_input_data(age_input, out_key, base_dir)

    # 6. Collect the updated frame, keyed by the folder that was written.
    result[[out_key]] <- frame
  }

  invisible(result)
}

#' Assign ages to non-synchronous horizons to use as older/younger than inputs
#'
#' Assigns a reference age to test horizons that were excluded from secondary synchronization so
#' they can be used as "older than"/"younger than" constraints in the age-depth
#' model. Each excluded horizon receives the synchronized age of the horizon
#' \emph{group} it belongs to: e.g. an excluded \code{"synchro-test-wrong"} layer
#' takes the age computed for its group \code{"synchro-test"}. The member-to-group
#' mapping comes from \code{horizon_groups} (as passed to \code{synchronize_ages()}). 
#' Exclusions are parsed with the same logic as
#' \code{synchronize_ages()} (named for per-record, unnamed for global).
#'
#' If an excluded horizon has no synchronized age for its group (it was excluded
#' altogether, so no group age was ever computed), it is skipped with a warning;
#' compute an age for it separately with \code{synchronize_ages()} using the method of
#' your choice (\code{"mean"}, \code{"ageofrecord"}, \code{"Bayesian"}, ...).
#'
#' @param adjusted_ages Output from \code{synchronize_ages()} containing synchronized ages.
#' @param nonsynchro_horizons Named or unnamed character vector specifying non-synchronized horizons
#'   (same format as in \code{synchronize_ages()}).
#' @param event_stats Output from \code{process_event_ages()} containing record metadata.
#' @param horizon_groups Named list mapping each horizon group to its member labels, as
#'   passed to \code{synchronize_ages()} / \code{compute_synchronized_ages()}. Used to map an
#'   excluded member label to the group whose synchronized age it should take (default:
#'   \code{NULL}).
#'
#' @return Named list with one element per excluded horizon, each containing a data frame
#'   with columns: \code{record}, \code{adjusted_age}, \code{adjusted_error}.
#'
#' @export
assign_nonsynchro_age <- function(adjusted_ages,
                                 nonsynchro_horizons,
                                 event_stats,
                                 horizon_groups = NULL) {
  
  excl_parsed         <- parse_excluded_horizons(nonsynchro_horizons, names(event_stats$processed))
  per_record_excluded <- excl_parsed$per_record_excluded
  global_excluded     <- excl_parsed$global_excluded

  # Reverse lookup: member horizon label -> group name (e.g. "synchro-test-wrong"
  # -> "synchro-test"), so an excluded member takes its group's synchronized age.
  group_of <- build_group_lookup(horizon_groups)

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
      
      # Resolve the horizon group this excluded horizon belongs to, and take that
      # group's synchronized age (e.g. "synchro-test-wrong" -> "synchro-test").
      grp <- if (!is.null(group_of) && ex_h %in% names(group_of)) group_of[[ex_h]] else ex_h

      if (!(grp %in% names(adjusted_ages))) {
        warning(sprintf(
          "Record '%s': excluded horizon '%s' has no synchronized age for its group '%s'. Skipping (assign one with synchronize_ages() if needed).",
          rec, ex_h, grp
        ))
        next
      }

      ref_df  <- adjusted_ages[[grp]]
      ref_row <- ref_df[ref_df$record == rec, , drop = FALSE]
      # Records whose member was fully excluded are absent from the group's data
      # frame; fall back to the group's combined age (uniform across its records).
      if (nrow(ref_row) == 0) ref_row <- ref_df[1, , drop = FALSE]
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
        adjusted_age = ref_row$adjusted_age[1],
        adjusted_error = ref_row$adjusted_error[1],
        stringsAsFactors = FALSE
      ))

      cat(sprintf(
        "Assigned %.1f +/- %.1f to excluded '%s' in record '%s' using group '%s'\n",
        ref_row$adjusted_age[1], ref_row$adjusted_error[1], ex_h, rec, grp
      ))
    }
  }
  
  return(out)
}

#' Load age input data and event depths
#'
#' Loads age input dates (radiocarbon and/or 210Pb) and event depths of all records from a folder of
#' CSV files and calculates maximum depths and sedimentation rates for each record. The function reads
#' every \code{.csv} file in the folder, where each file represents one sediment record. For records
#' missing required columns ('depth' or 'C14_age'), warnings are issued and NA values are returned.
#'
#' @param folder_path Character string specifying the path to the folder containing your input file
#'   (default: \code{"."}, i.e. the working directory set by \code{syncer_setup()}).
#' @param file_name Character string specifying the name of the subfolder (inside \code{folder_path})
#'   that holds your input CSV files (default: \code{"record_data_input"}). This folder holds one
#'   \strong{CSV file per record} (named \code{<record>.csv}), and each file gives, for every dated
#'   sample and for the record top, the \code{depth}, \code{C14_age} and \code{C14_error}, together
#'   with the radiocarbon calibration curve code in the \code{cc} column (\code{0} = no calibration /
#'   calendar ages, \code{1} = IntCal20, \code{2} = Marine20, \code{3} = SHCal20). The file must also
#'   list the depths of the considered event deposits; these depths may be given either as
#'   event-free depth or as total depth.
#'
#' @return A named list with three elements:
#'   \itemize{
#'     \item \code{record_data}: Named list of data frames, one per record (each CSV file becomes one element)
#'     \item \code{max_depths}: Named numeric vector containing the maximum depth value for each record
#'     \item \code{sedrates}: Named numeric vector containing an estimated sedimentation rate for each record
#'           (calculated as oldest C14 age divided by its depth)
#'   }
#'
#' @export
read_record_data <- function(folder_path = ".",
                            file_name="record_data_input") {

  # Construct full path to the folder holding one CSV file per record
  input_dir <- file.path(folder_path, file_name)

  if (!dir.exists(input_dir)) {
    stop(paste("Input folder not found:", input_dir))
  }

  # Get all CSV files (each file = one sediment record)
  csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
  record_names <- tools::file_path_sans_ext(basename(csv_files))

  if (length(csv_files) == 0) {
    stop(paste("No CSV files found in:", input_dir))
  }

  # Initialize output structures
  record_data <- list()
  n_records <- length(record_names)
  max_depths <- numeric(n_records)
  sedrates <- numeric(n_records)
  names(max_depths) <- record_names
  names(sedrates) <- record_names

  # Process each record (CSV file)
  for (i in seq_along(record_names)) {
    record <- record_names[i]
    # Read CSV data
    df <- utils::read.csv(csv_files[i], stringsAsFactors = FALSE)
    record_data[[record]] <- df

    # Extract maximum depth if 'depth' column exists
    if ("depth" %in% colnames(df)) {
      max_depths[[record]] <- max(df$depth, na.rm = TRUE)
    } else {
      max_depths[[record]] <- NA_real_
      warning(paste("Record", record, "does not have a 'depth' column"))
    }

    # Calculate sedimentation rate (oldest age / depth at oldest age)
    if (all(c("depth", "C14_age") %in% colnames(df))) {
      # Find oldest (maximum) 14C age
      oldest_index <- which.max(df$`C14_age`)
      oldest_age <- as.numeric(df$`C14_age`[oldest_index])
      depth_at_oldest <- as.numeric(df$depth[oldest_index])

      # Calculate rate if values are valid
      if (is.finite(oldest_age) && is.finite(depth_at_oldest) && depth_at_oldest != 0) {
        sedrates[[record]] <- oldest_age / depth_at_oldest
      } else {
        sedrates[[record]] <- NA_real_
        warning(paste("Invalid values for sedrate in record", record))
      }
    } else {
      sedrates[[record]] <- NA_real_
      warning(paste("Record", record, "is missing 'depth' or 'C14_age' column"))
    }
  }

  return(list(record_data = record_data,
              max_depths = max_depths,
              sedrates = sedrates))
}

#' Read age results of age-depth modelling for each of the event horizons
#'
#' Reads a folder of CSV files containing all age information for events, one file per record,
#' and returns a named list of data frames. This is the inverse operation of
#' \code{write_age_output_data()}, and allows you to reload previously saved results without recalculating or 
#' load age info into SyncER in case the age-depth models were not constructed using \emph{rbacon}/\emph{rplum}.
#'
#' @param folder_path Character string specifying the location of the folder to be read
#'   (default: \code{syncer_output_dir()}, i.e. the \code{SyncER_outputs} folder in the working directory).
#' @param synced Character string suffix for the input folder name (default: "");
#'   use "_synced" for synchronized data.
#'
#' @return A named list of data frames, where each element corresponds to one record. 
#' List names match the CSV file names (without extension).
#'
#' @examples
#' \dontrun{
#' # Read non-synchronized data
#' out_data <- read_age_data("path/to/folder")
#'
#' # Read synchronized data
#' out_data_synced <- read_age_data("path/to/folder", synced = "_synced")
#'
#' # Use with process_event_ages
#' event_stats <- process_event_ages(out_data, event_deposits = c("tephra", "flood"))
#' }
#'
#' @export
read_age_data <- function(folder_path = syncer_output_dir(),
                            synced = "") {

  # Construct input directory
  input_dir <- file.path(folder_path, paste0("out_data_ages", synced))

  # Check if folder exists
  if (!dir.exists(input_dir)) {
    stop(paste("Input folder not found:", input_dir))
  }

  # Get all CSV files (each file = one record)
  csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
  record_names <- tools::file_path_sans_ext(basename(csv_files))

  if (length(csv_files) == 0) {
    stop(paste("No CSV files found in:", input_dir))
  }

  # Read all CSV files into a named list
  out_data <- setNames(
    lapply(csv_files, function(f) {
      df <- utils::read.csv(f, stringsAsFactors = FALSE)
      if (!is.data.frame(df)) {
        stop(paste("File", f, "did not return a data frame"))
      }
      return(df)
    }),
    record_names
  )

  cat("Data successfully read from:", input_dir, "\n")
  cat("Loaded", length(out_data), "record(s):", paste(names(out_data), collapse = ", "), "\n")

  return(out_data)
}
