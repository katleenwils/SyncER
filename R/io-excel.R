#' Load Record Data
#'
#' Loads radiocarbon dating records from an Excel file and calculates maximum depths
#' and sedimentation rates for each record.
#'
#' @param folder_path Character string specifying the path to the folder containing your input file.
#' @param file_name Character string specifying the name of your input Excel file (default: "record_data_input.xlsx").
#'
#' @return A named list with three elements:
#'   \itemize{
#'     \item \code{record_data}: Named list of data frames, one per record (each sheet becomes one element)
#'     \item \code{max_depths}: Named numeric vector containing the maximum depth value for each record
#'     \item \code{sedrates}: Named numeric vector containing the sedimentation rate for each record
#'           (calculated as oldest C14 age divided by its depth)
#'   }
#'
#' @details The function reads all sheets from the Excel file, where each sheet represents one
#'   sediment record. For records missing required columns ('depth' or 'C14_age'), appropriate
#'   warnings are issued and NA values are returned.
#'
#' @export
load_excel_data <- function(folder_path,
                            file_name="record_data_input.xlsx") {

  # Construct full path to Excel file
  excel_file <- file.path(folder_path, file_name)

  # Get all sheet names (each sheet = one sediment record)
  record_names <- readxl::excel_sheets(excel_file)

  # Initialize output structures
  record_data <- list()
  n_records <- length(record_names)
  max_depths <- numeric(n_records)
  sedrates <- numeric(n_records)
  names(max_depths) <- record_names
  names(sedrates) <- record_names

  # Process each record (sheet)
  for (record in record_names) {
    # Read sheet data
    df <- readxl::read_excel(excel_file, sheet = record)
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

#' Export List of Data Frames to an Excel File
#'
#' Exports a named list of data frames to a single Excel file with multiple worksheets.
#'
#' @param out_data Named list of data frames where each element will become a separate worksheet.
#' @param folder_path Character string specifying the location where the Excel file should be saved.
#' @param synced Character string suffix for the output filename (default: "");
#'   use "_synced" for synchronized data.
#'
#' @return No return value. Writes an Excel file and prints a success message with the file path.
#'
#' @details Each list element is written to a separate sheet in the output Excel file, with
#'   sheet names matching the list element names.
#'
#' @importFrom writexl write_xlsx
#' @export
export_to_excel <- function(out_data,
                            folder_path,
                            synced="") {

  # Convert all list elements to data frames
  records_list <- lapply(out_data, as.data.frame)

  # Construct output filename
  output_file <- file.path(folder_path, paste0("out_data_ages", synced, ".xlsx"))

  # Write to Excel (one sheet per record)
  writexl::write_xlsx(records_list, path = output_file)
  cat("Data successfully exported to:", output_file, "\n")
}

#' Read List of Data Frames from an Excel File
#'
#' Reads an Excel file containing multiple worksheets and returns a named list of data frames.
#' This is the inverse operation of \code{export_to_excel()}.
#'
#' @param folder_path Character string specifying the location of the Excel file to be read.
#' @param synced Character string suffix for the input filename (default: "");
#'   use "_synced" for synchronized data.
#'
#' @return A named list of data frames, where each element corresponds to one worksheet
#'   from the Excel file. List names match the worksheet names.
#'
#' @details Reads all sheets from the Excel file "out_data_ages.xlsx" (or
#'   "out_data_ages_synced.xlsx" if synced="_synced") and returns them as a named list,
#'   identical in structure to the \code{out_data} parameter used in \code{export_to_excel()}.
#'   This allows you to reload previously saved results without recalculating.
#'
#' @examples
#' \dontrun{
#' # Read non-synchronized data
#' out_data <- read_from_excel("path/to/folder")
#'
#' # Read synchronized data
#' out_data_synced <- read_from_excel("path/to/folder", synced = "_synced")
#'
#' # Use with process_event_ages
#' event_stats <- process_event_ages(out_data, event_deposits = c("tephra", "flood"))
#' }
#'
#' @importFrom readxl read_excel excel_sheets
#' @export
read_from_excel <- function(folder_path,
                            synced = "") {

  # Construct input filename
  input_file <- file.path(folder_path, paste0("out_data_ages", synced, ".xlsx"))

  # Check if file exists
  if (!file.exists(input_file)) {
    stop(paste("Excel file not found:", input_file))
  }

  # Get all sheet names
  sheet_names <- readxl::excel_sheets(input_file)

  # Read all sheets into a named list
  out_data <- setNames(
    lapply(sheet_names, function(sheet) {
      df <- readxl::read_excel(input_file, sheet = sheet)
      if (!is.data.frame(df)) {
        stop(paste("Sheet", sheet, "did not return a data frame"))
      }
      return(df)
    }),
    sheet_names
  )

  cat("Data successfully read from:", input_file, "\n")
  cat("Loaded", length(out_data), "record(s):", paste(names(out_data), collapse = ", "), "\n")

  return(out_data)
}
