#' Read Raw Bacon Output Files
#'
#' Reads raw lines from Bacon \code{.out} files for a set of records into a named list.
#' This is the I/O counterpart to \code{extract_event_ages()}: call this first to obtain
#' raw data, then pass the result to \code{extract_event_ages()} for computation, and
#' finally call \code{export_to_excel()} to save the processed ages.
#'
#' @param folder_path Character string giving the parent directory that contains
#'   the per-record Bacon output sub-folders (default: \code{"."}, i.e. the working
#'   directory itself, matching the \code{coredir} used by \emph{rbacon}/\emph{rplum}
#'   and \code{age_model_input()}).
#' @param synced Character string suffix used to identify synchronized output folders
#'   (default: \code{""}; non-empty values such as \code{"_synced"} select only folders
#'   ending with that suffix).
#' @param max_depths Named numeric vector of original core depths (cm), used only for
#'   progress reporting. Pass \code{NULL} to skip depth messages.
#' @param excluded_depths Optional named list of excluded depth intervals per record,
#'   used only for progress reporting.
#'
#' @return Named list of character vectors (one element per folder / record), where each
#'   element contains the raw text lines of the corresponding \code{.out} file.
#'
#' @export
read_bacon_output <- function(folder_path = ".",
                              synced = "",
                              max_depths = NULL,
                              excluded_depths = NULL) {

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
      exclusions_local   <- if (!is.null(excluded_depths)) excluded_depths[[record_name_local]] else NULL

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

#' Compute Event Ages from Raw Bacon Output Lines
#'
#' Performs all age interpolation and computation logic from raw Bacon .out file content,
#' without any file reading or writing.
#'
#' @param raw_out_data Named list of character vectors (raw lines from each .out file),
#'   one element per record folder.
#' @param record_data Named list of data frames representing each record's metadata
#'   (output from \code{load_excel_data()}).
#' @param event_types Character vector of all event type names present in your dataset.
#' @param max_depths Named vector containing the depths to which the age models are calculated for each record.
#' @param isochrons Character vector of deposit names that are known to be synchronous.
#' @param test_horizons Character vector of event deposits for which synchronicity should be tested.
#' @param excluded_depths Optional named list of excluded depth intervals per record.
#' @param thick The Bacon section thickness; single numeric, named list per record, or NULL
#'   (default), in which case it is derived per record from the model as
#'   \code{max_depths[record] / (n_cols - 3)}, i.e. the modelled depth range divided
#'   by the number of Bacon sections.
#' @param synced Character string suffix to identify synchronized folders (default: "").
#'
#' @return Named list of data frames, one per processed record folder, containing
#'   depth columns plus one column per event with interpolated ages.
#'
#' @export
extract_event_ages <- function(raw_out_data,
                               record_data,
                               event_types,
                               max_depths,
                               isochrons,
                               test_horizons,
                               excluded_depths = NULL,
                               thick = NULL,
                               synced = "") {

  out_data <- list()

  # Helper function: Calculate total excluded depth up to a given depth
  calculate_excluded_depth <- function(depth, exclusions) {
    if (is.null(exclusions) || length(exclusions) == 0) {
      return(0)
    }

    if (length(exclusions) %% 2 != 0) {
      stop("excluded_depths must contain pairs of values (top, bottom, top, bottom, ...)")
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
    exclusions <- if (!is.null(excluded_depths)) excluded_depths[[record_name]] else NULL

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
    # one. A Bacon .out for this record has (n_cols - 3) depth sections spanning
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

    # Calculate total excluded depth
    total_excluded <- calculate_excluded_depth(max_depth_original, exclusions)

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

#' Extract Event Ages from Bacon Output Files
#'
#' I/O wrapper that reads Bacon \code{.out} files and returns a named list of
#' processed event ages. Combines \code{read_bacon_output()} and
#' \code{extract_event_ages()} into a single call. When \code{reload_existing}
#' is \code{TRUE} the previously exported Excel file is returned instead.
#'
#' @param folder_path Character string giving the parent directory that contains
#'   per-record Bacon output sub-folders (default: \code{"."}, i.e. the working
#'   directory itself, matching the \code{coredir} used by \emph{rbacon}/\emph{rplum}
#'   and \code{age_model_input()}). Ignored when \code{reload_existing = TRUE}.
#' @param record_data Named list of data frames representing each record's
#'   metadata (output from \code{load_excel_data()}).
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
#' @param excluded_depths Optional named list of excluded depth intervals per record.
#' @param thick The Bacon section thickness used during age-depth modelling;
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
extract_ages <- function(folder_path = ".",
                         record_data,
                         event_types,
                         max_depths,
                         isochrons,
                         test_horizons,
                         synced = "",
                         reload_existing = FALSE,
                         excluded_depths = NULL,
                         thick = NULL,
                         output_dir = syncer_output_dir()) {

  if (reload_existing) {
    return(read_from_excel(output_dir, synced = synced))
  }

  raw_out_data <- read_bacon_output(
    folder_path     = folder_path,
    synced          = synced,
    max_depths      = max_depths,
    excluded_depths = excluded_depths
  )

  extract_event_ages(
    raw_out_data    = raw_out_data,
    record_data     = record_data,
    event_types     = event_types,
    max_depths      = max_depths,
    isochrons       = isochrons,
    test_horizons   = test_horizons,
    excluded_depths = excluded_depths,
    thick           = thick,
    synced          = synced
  )
}
