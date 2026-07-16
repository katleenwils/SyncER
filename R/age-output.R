#' Read Bacon/Plum age-depth modelling output files
#'
#' Reads raw lines from \emph{rbacon}/\emph{rplum} \code{.out} files for a set of records into a named list.
#' Call this function to obtain raw age data, then pass the result to \code{compute_event_ages()} for computation of event ages, 
#' and finally call \code{write_age_output_data()} to save the processed ages in your folder.
#'
#' @param folder_path Character string giving the parent directory that contains
#'   the per-record \emph{rbacon}/\emph{rplum} output sub-folders (default: \code{"."}, i.e. the working
#'   directory itself, matching the \code{coredir} used by \emph{rbacon}/\emph{rplum}
#'   and \code{age_model_input()}).
#' @param synced Character string suffix used to identify synchronized output folders
#'   (default: \code{""}; non-empty values such as \code{"_synced"} select only folders
#'   ending with that suffix).
#' @param max_depths Named numeric vector of original core depths (cm), used only for
#'   progress reporting. Pass \code{NULL} to skip depth messages.
#' @param instantaneous_event_depths Optional named list of excluded depth intervals per record,
#'   used only for progress reporting.
#'
#' @return Named list of character vectors (one element per folder / record), where each
#'   element contains the raw text lines of the corresponding \code{.out} file.
#'
#' @export
read_age_model_output <- function(folder_path = ".",
                              synced = "",
                              max_depths = NULL,
                              instantaneous_event_depths = NULL) {

  # Discover valid folders based on synced status
  all_folders <- list.dirs(folder_path, recursive = FALSE, full.names = FALSE)

  if (synced != "") {
    pattern <- paste0(synced, "$")
    valid_folders <- all_folders[grepl(pattern, all_folders)]
  } else {
    valid_folders <- all_folders[!grepl("synced", all_folders)]
  }

  # Find all .out files ending with a number
  out_files <- list.files(
    folder_path,
    pattern = "[0-9]\\.out$",
    full.names = TRUE,
    recursive = TRUE
  )

  # Keep only files from valid folders
  out_files <- out_files[vapply(out_files, function(x) {
    folder_name <- basename(dirname(x))
    folder_name %in% valid_folders
  }, logical(1))]

  # Read raw lines from each .out file; report progress
  raw_out_data <- list()
  for (out_file in out_files) {
    folder_name <- basename(dirname(out_file))
    cat("Processing:", folder_name, "\n")

    temp_data <- readr::read_lines(out_file)

    # Report detected separator
    first_line <- temp_data[nchar(trimws(temp_data)) > 0][1]
    sep_label  <- if (grepl(",", first_line)) "comma" else if (grepl("\t", first_line)) "tab" else "whitespace"
    cat("  Detected separator:", sep_label, "\n")

    # Extract record name for depth reporting
    record_name_local <- if (synced != "") {
      sub(paste0(synced, "$"), "", folder_name)
    } else {
      folder_name
    }

    if (!is.null(max_depths)) {
      max_depth_original <- max_depths[record_name_local]
      exclusions_local   <- if (!is.null(instantaneous_event_depths)) instantaneous_event_depths[[record_name_local]] else NULL

      cat("  Original max depth:", max_depth_original, "cm\n")

      if (!is.null(exclusions_local) && length(exclusions_local) > 0) {
          exclusion_matrix <- matrix(exclusions_local, ncol = 2, byrow = TRUE)
          total_excl <- sum(exclusion_matrix[, 2] - exclusion_matrix[, 1])
          if (total_excl > 0) {
            cat("  Total event depth:", total_excl, "cm\n")
            cat("  Event-free depth:", max_depth_original - total_excl, "cm\n")
          }
        }
      }

    raw_out_data[[folder_name]] <- temp_data
  }

  return(raw_out_data)
}

#' Prepare age data for log-ratio transformation and calculate basic statistics
#'
#' Extracts event-specific age columns from \emph{rbacon}/\emph{rplum} output, applies offset correction,
#' and calculates summary statistics. Filters columns matching event deposit patterns, 
#' applies offset to avoid negative ages for post-1950 deposits, and ensures the record top is at age 0. 
#' This function takes the direct output from \code{load_event_ages()} and processes it in memory.
#'
#' @param out_data Named list of data frames or data.tables containing age data
#'   (output from \code{load_event_ages()}).
#' @param event_deposits Character vector of event deposit names to extract (typically
#'   \code{isochrons} or \code{test_events}). Columns whose name starts with any of these names
#'   are selected.
#' @param offset Numeric value added to every age so that the age reference point becomes, for
#'   example, the year the (most recent) records were retrieved rather than 1950 AD (default:
#'   \code{bp_datum()}, i.e. the current year minus 1950). Because the synchronicity test uses
#'   log-transformations, the age dataset may not contain zero or negative values. Set it manually 
#'   if you want a different reference point, but use the \emph{same} offset consistently throughout the workflow.
#'
#' @return A named list with two elements:
#'   \itemize{
#'     \item \code{processed}: List of data frames with offset-corrected ages for each record
#'     \item \code{summaries}: List of summary data frames containing min, max, mean, and
#'           sd (sigma) for each event column
#'   }
#'
#' @importFrom magrittr %>%
#' @export
process_event_ages <- function(out_data,
                               event_deposits,
                               offset = bp_datum()) {
  
  # Validate input
  if (!is.list(out_data) || length(out_data) == 0) {
    stop("out_data must be a non-empty list (output from load_event_ages)")
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

#' Interpolate event ages from MCMC runs of Raw Bacon/Plum output files
#'
#' Performs all age interpolation and computation logic from raw \emph{rbacon}/\emph{rplum} .out file content.
#'
#' @param raw_out_data Named list of character vectors (raw lines from each .out file),
#'   one element per record folder.
#' @param record_data Named list of data frames representing each record's metadata
#'   (output from \code{read_record_data()}).
#' @param event_types Character vector of all event type names present in your dataset.
#' @param max_depths Named vector containing the depths to which the age models are calculated for each record.
#' @param isochrons Character vector of deposit names that are known to be synchronous.
#' @param test_horizons Character vector of event deposits for which synchronicity should be tested.
#' @param instantaneous_event_depths Optional named list of depth intervals classified as instantaneous deposits per record
#' that should not be considered for event-free depths.
#' @param thick The \emph{rbacon}/\emph{rplum} section thickness; single numeric, named list per record, or NULL
#'   (default), in which case it is derived per record from the model as
#'   \code{max_depths[record] / (n_cols - 3)}, i.e. the modelled depth range divided
#'   by the number of \emph{rbacon}/\emph{rplum} sections.
#' @param synced Character string suffix to identify synchronized folders (default: "").
#'
#' @return Named list of data frames, one per processed record folder, containing
#'   depth columns plus one column per event with interpolated ages.
#'
#' @export
compute_event_ages <- function(raw_out_data,
                               record_data,
                               event_types,
                               max_depths,
                               isochrons,
                               test_horizons,
                               instantaneous_event_depths = NULL,
                               thick = NULL,
                               synced = "") {

  out_data <- list()

  # Helper function: Calculate total excluded depth up to a given depth
  calculate_excluded_depth <- function(depth, exclusions) {
    if (is.null(exclusions) || length(exclusions) == 0) {
      return(0)
    }

    if (length(exclusions) %% 2 != 0) {
      stop("instantaneous_event_depths must contain pairs of values (top, bottom, top, bottom, ...)")
    }

    exclusion_matrix <- matrix(exclusions, ncol = 2, byrow = TRUE)
    colnames(exclusion_matrix) <- c("top", "bottom")

    total_excluded <- 0
    for (i in 1:nrow(exclusion_matrix)) {
      top <- exclusion_matrix[i, "top"]
      bottom <- exclusion_matrix[i, "bottom"]

      if (bottom <= depth) {
        total_excluded <- total_excluded + (bottom - top)
      }
    }

    return(total_excluded)
  }

  for (folder_name in names(raw_out_data)) {

    temp_data <- raw_out_data[[folder_name]]

    # Extract record name
    if (synced != "") {
      record_name <- sub(paste0(synced, "$"), "", folder_name)
    } else {
      record_name <- folder_name
    }

    record_df  <- record_data[[record_name]]
    exclusions <- if (!is.null(instantaneous_event_depths)) instantaneous_event_depths[[record_name]] else NULL

    # Resolve the thickness multiplier for this record. When the caller supplies
    # `thick` we honour it; otherwise it is derived from the model itself below,
    # once n_cols and the record's max depth are known.
    thick_value <- NULL
    if (!is.null(thick)) {
      if (is.list(thick)) {
        if (!record_name %in% names(thick)) {
          stop(paste("thick is a list but contains no entry for record:", record_name))
        }
        thick_value <- thick[[record_name]]
      } else {
        thick_value <- thick
      }
    }

    # Detect separator by checking first non-empty line
    first_line <- temp_data[nchar(trimws(temp_data)) > 0][1]
    sep <- if (grepl(",", first_line)) "," else if (grepl("\t", first_line)) "\t" else "\\s+"

    temp_data_dt <- as.data.frame(
      do.call(rbind, strsplit(trimws(temp_data), sep)),
      stringsAsFactors = FALSE
    )
    names(temp_data_dt) <- as.character(0:(ncol(temp_data_dt) - 1))

    n_cols <- ncol(temp_data_dt)

    # Extract maximum depth for this record
    max_depth_original <- max_depths[record_name]

    # Derive the section thickness from the model when the caller did not supply
    # one. A Bacon/Plum .out for this record has (n_cols - 3) depth sections spanning
    # 0..max_depth, so each section is max_depth / (n_cols - 3) deep. This replaces
    # the previous heuristic, which read the last number of the folder name (e.g.
    # "core1" -> 1) and thus used the core index as the thickness -- collapsing
    # max_depth for low-index cores and silently dropping their deeper events.
    if (is.null(thick_value)) {
      if (is.na(max_depth_original) || (n_cols - 3) <= 0) {
        stop(paste("Cannot derive section thickness for record:", record_name,
                   "- supply `thick` explicitly."))
      }
      thick_value <- max_depth_original / (n_cols - 3)
    }

    # Convert column names to actual depths
    new_headers <- (1:(n_cols - 1)) * thick_value
    max_depth <- new_headers[length(new_headers) - 2]

    # Update column names to depths
    names(temp_data_dt)[2:n_cols] <- as.character(new_headers)

    # Remove last two columns
    temp_data_dt <- temp_data_dt[, 1:(n_cols - 2)]

    # Convert all columns to numeric
    temp_data_dt[] <- lapply(temp_data_dt, as.numeric)

    # Initialize data frame for event ages
    event_results_dt <- data.frame(matrix(nrow = nrow(temp_data_dt), ncol = 0))

    # Process each event type
    for (event in event_types) {

      event_data <- record_df[
        stringr::str_detect(record_df$event, paste0("^", event, "($|[_-]|[a-zA-Z]|\\d)")),
      ]

      if (nrow(event_data) > 0) {
        for (i in 1:nrow(event_data)) {
          event_depth_original <- event_data$depth[i]

          excluded_above <- calculate_excluded_depth(event_depth_original, exclusions)
          event_depth <- event_depth_original - excluded_above

          if (event_depth > max_depth) {
            warning(paste(event_data$event[i], "at depth", event_depth,
                          "exceeds max depth", max_depth, "- skipping"))
            next
          }

          if (event %in% isochrons || event %in% test_horizons) {
            column_name <- event_data$event[i]
          } else {
            column_name <- paste0(event, "_", event_depth)
          }

          event_col_values <- vapply(seq_len(nrow(temp_data_dt)), function(j) {
            row_values <- as.numeric(temp_data_dt[j, ])
            row_depths <- as.numeric(names(temp_data_dt))

            valid_columns <- which(row_depths <= event_depth)
            last_column <- min(which(row_depths > event_depth))

            if (length(valid_columns) > 1) {
              age_values <- row_values[valid_columns[-1]] * thick_value
              sum_values <- sum(age_values)
            } else {
              sum_values <- 0
            }

            last_value_sum <- row_values[last_column] *
              (event_depth - row_depths[length(valid_columns)])

            row_values[valid_columns[1]] + sum_values + last_value_sum
          }, numeric(1))

          event_results_dt[[column_name]] <- event_col_values
        }
      }
    }

    # Combine depth columns with event age columns
    temp_data_dt <- cbind(temp_data_dt, event_results_dt)

    out_data[[folder_name]] <- temp_data_dt
  }

  return(out_data)
}

#' Load and store event ages from Bacon/Plum output files
#'
#' Wrapper function that reads \emph{rbacon}/\emph{rplum} \code{.out} files and returns a named list of
#' processed event ages. Combines \code{read_age_model_output()} and
#' \code{compute_event_ages()} into a single call. When \code{reload_existing}
#' is \code{TRUE}, the previously exported Excel file is returned instead.
#'
#' @param folder_path Character string giving the parent directory that contains
#'   per-record \emph{rbacon}/\emph{rplum} output sub-folders (default: \code{"."}, i.e. the working
#'   directory itself, matching the \code{coredir} used by \emph{rbacon}/\emph{rplum}
#'   and \code{age_model_input()}). Ignored when \code{reload_existing = TRUE}.
#' @param record_data Named list of data frames representing each record's
#'   metadata (output from \code{read_record_data()}).
#' @param event_types Character vector of all event type names present in your
#'   dataset.
#' @param max_depths Named numeric vector of maximum depths per record.
#' @param isochrons Character vector of deposit names known to be synchronous.
#' @param test_horizons Character vector of event deposits for which
#'   synchronicity should be tested.
#' @param synced Character string suffix identifying synchronized output folders
#'   (default: \code{""}).
#' @param reload_existing Logical; when \code{TRUE} reads from an already-exported
#'   Excel file instead of re-processing \code{.out} files (default: \code{FALSE}).
#' @param instantaneous_event_depths Optional named list of depth intervals classified as instantaneous deposits per record
#' that should not be considered for event-free depths.
#' @param thick The \emph{rbacon}/\emph{rplum} section thickness used during age-depth modelling;
#'   single numeric, named list per record, or \code{NULL} (default), in which case
#'   it is derived per record from the model as \code{max_depths[record] / (n_cols - 3)}.
#' @param output_dir Character string specifying where the previously-exported
#'   \code{out_data_ages.xlsx} lives; only used when \code{reload_existing = TRUE}
#'   (default: \code{syncer_output_dir()}, i.e. the \code{SyncER_outputs} folder in
#'   the working directory).
#'
#' @return Named list of data frames, one per processed record, containing depth
#'   columns plus one column per event with interpolated ages.
#'
#' @export
load_event_ages <- function(folder_path = ".",
                         record_data,
                         event_types,
                         max_depths,
                         isochrons,
                         test_horizons,
                         synced = "",
                         reload_existing = FALSE,
                         instantaneous_event_depths = NULL,
                         thick = NULL,
                         output_dir = syncer_output_dir()) {

  if (reload_existing) {
    return(read_age_data(output_dir, synced = synced))
  }

  raw_out_data <- read_age_model_output(
    folder_path     = folder_path,
    synced          = synced,
    max_depths      = max_depths,
    instantaneous_event_depths = instantaneous_event_depths
  )

  compute_event_ages(
    raw_out_data    = raw_out_data,
    record_data     = record_data,
    event_types     = event_types,
    max_depths      = max_depths,
    isochrons       = isochrons,
    test_horizons   = test_horizons,
    instantaneous_event_depths = instantaneous_event_depths,
    thick           = thick,
    synced          = synced
  )
}

#' Write age data (including event ages) into an Excel file
#'
#' Each record age information (including event ages) is written to a separate sheet in the output Excel file, with
#'   sheet names matching the list element names, and each row represents a single MCMC simulation.
#'
#' @param out_data Named list of data frames where each element will become a separate worksheet.
#' @param folder_path Character string specifying the location where the Excel file should be saved
#'   (default: \code{syncer_output_dir()}, i.e. the \code{SyncER_outputs} folder in the working directory).
#' @param synced Character string suffix for the output filename (default: "");
#'   use "_synced" for synchronized data.
#'
#' @return No return value. Writes an Excel file and prints a success message with the file path.
#'
#' @importFrom writexl write_xlsx
#' @export
write_age_output_data <- function(out_data,
                            folder_path = syncer_output_dir(),
                            synced="") {

  # Convert all list elements to data frames
  records_list <- lapply(out_data, as.data.frame)

  # Construct output filename
  output_file <- file.path(folder_path, paste0("out_data_ages", synced, ".xlsx"))

  # Write to Excel (one sheet per record)
  writexl::write_xlsx(records_list, path = output_file)
  cat("Data successfully exported to:", output_file, "\n")
}
