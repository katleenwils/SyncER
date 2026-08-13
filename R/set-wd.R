#' Set up SyncER working directory and output folder
#'
#' Sets the working directory \emph{SyncER} should use, and ensures a \code{SyncER_outputs} folder exists there
#' for output files. Call it yourself once, after \code{library(SyncER)}, at the top of your analysis
#' script.
#'
#' @param wd Character string giving the working directory to use. If
#'   omitted, first checks \code{options(SyncER.wd = ...)}; if that is unset
#'   and the session is interactive, asks for a path via \code{readline()}.
#'   Leave blank (or pass \code{""}) to keep the current working directory.
#'
#' @return Invisibly returns the path to the \code{SyncER_outputs} directory.
#'
#' @export
syncer_setup <- function(wd = getOption("SyncER.wd")) {
  if (is.null(wd)) {
    wd <- if (interactive()) {
      readline(prompt = "SyncER: enter the working directory to use (leave blank to keep the current one): ")
    } else {
      ""
    }
  }

  if (nzchar(wd)) {
    setwd(wd)
  }
  message("Working directory: ", getwd())

  out_dir <- syncer_output_dir()
  message("SyncER output will be saved to: ", out_dir)
  invisible(out_dir)
}

#' Set up SyncER output directory
#'
#' Returns the path to the \code{SyncER_outputs} folder inside the current
#' working directory, creating it if it does not already exist. Every \emph{SyncER}
#' function that writes files (CSVs, PDFs) defaults its
#' output-location argument to this folder, so that all package output ends
#' up in one predictable place unless the user explicitly overrides it.
#'
#' @return Character string: path to the \code{SyncER_outputs} directory.
#'
#' @keywords internal
syncer_output_dir <- function() {
  out_dir <- file.path(getwd(), "SyncER_outputs")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_dir
}
