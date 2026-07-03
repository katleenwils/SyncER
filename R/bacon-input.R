#' Build Bacon Input Data Structures
#'
#' Prepares radiocarbon and Pb210 data frames in the format required by the Bacon
#' age-depth modeling software, without performing any file I/O.
#'
#' @param record_data A data frame containing record data with at least columns: \code{event},
#'   \code{depth}, \code{C14_age}, \code{C14_error}, and \code{cc} (calibration curve indicator).
#' @param record_name Character string used to identify the record.
#' @param radiocarbon_sample_names All labels given to the radiocarbon samples you want to include
#' @param lead_sample_names All labels given to the Pb210 samples you want to include
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

#' Save Records as CSV Files for Bacon
#'
#' Saves Bacon-formatted radiocarbon (and optionally Pb210) input data from a
#' \code{build_bacon_input()} result to CSV file(s) in the appropriate output folder.
#'
#' @param bacon_input Named list with elements \code{c14} and \code{pb210} as returned by
#'   \code{build_bacon_input()}.
#' @param record_name Character string used to name the output file and folder.
#' @param base_dir Character string specifying the output directory where Bacon input files will be saved.
#' @param adjusted_ages Optional list output from \code{set_to_zero()} containing synchronized ages
#'   per horizon (default: NULL). When provided, output is written to a \code{"_synced"} subfolder.
#'
#' @return No return value. Writes one or two CSV files and prints the saved paths.
#'
#' @details If Pb210 data is present in \code{bacon_input}, two CSV files are saved: one for
#'   C14 data (with \code{_C14} suffix) and one for Pb210 data (plain name). If only C14 data
#'   is present, a single CSV file is saved under the plain record name.
#'
#' @export
age_model_input <- function(bacon_input,
                            record_name,
                            base_dir,
                            adjusted_ages = NULL) {

  # ── Helper: resolve output folder & file path ──────────────────────────────
  make_path <- function(suffix = "") {

    synced_tag <- if (!is.null(adjusted_ages)) "_synced" else ""

    if (!is.null(adjusted_ages)) {
      folder <- file.path(base_dir, paste0(record_name, "_synced"))
    } else {
      folder <- file.path(base_dir, record_name)
    }

    if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)

    file.path(folder, paste0(record_name, synced_tag, suffix, ".csv"))
  }

  sample_events <- bacon_input$c14
  has_lead      <- !is.null(bacon_input$pb210)

  # If Pb210 data is present, save with _C14 suffix; otherwise save as plain record.csv
  c14_path <- if (has_lead) make_path("_C14") else make_path()
  write.csv(sample_events, c14_path, row.names = FALSE, quote = FALSE, na = "")
  message("C14 data for ", record_name, " saved to: ", c14_path)

  if (has_lead) {
    pb_path <- make_path()
    write.csv(bacon_input$pb210, pb_path, row.names = FALSE, quote = FALSE, na = "")
    message("Pb210 data for ", record_name, " saved to: ", pb_path)
  }
}
