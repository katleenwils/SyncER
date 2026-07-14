#' Process Event Ages
#'
#' Extracts event-specific age columns from Bacon output, applies offset correction,
#' and calculates summary statistics.
#'
#' @param out_data Named list of data frames or data.tables containing age data
#'   (output from \code{extract_ages()}).
#' @param event_deposits Character vector of event deposit names to extract.
#' @param offset Numeric value to add to ages for offset correction (default:
#'   \code{bp_datum()}, i.e. the current year minus 1950).
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
