#' Set up SyncER working directory and output folder
#'
#' Optional first step of a \emph{SyncER} script: choose the folder \emph{SyncER} should work in
#' (holding your \code{record_data_input} folder and the age-depth model folders, and receiving all
#' output in a \code{SyncER_outputs} subfolder). The choice is remembered for the rest of the R session
#' (via \code{options(SyncER.wd = ...)}) and used as the default location by the other \emph{SyncER}
#' functions, so you only need to set it once, at the top of your script. If the folder does not
#' exist yet, it is created. \strong{Your R working directory (\code{getwd()}) is never changed.}
#' If you never call this function, input is read from the current working directory and output is
#' written to a per-session temporary directory (see \code{\link{syncer_output_dir}}), so
#' \emph{SyncER} never writes to your file space without being asked to.
#'
#' @param wd Character string giving the directory \emph{SyncER} should use. If
#'   omitted, first checks \code{options(SyncER.wd = ...)}; if that is unset
#'   and the session is interactive, asks for a path via \code{readline()}.
#'   Leave blank (or pass \code{""}) to use the current working directory.
#' @param verbose Logical; if \code{TRUE} (default), progress and summary messages are
#'   printed via \code{message()}. Set to \code{FALSE} to suppress them.
#'
#' @return Invisibly returns the path to the \code{SyncER_outputs} directory.
#'
#' @examples
#' # Use a temporary directory (never write to your own file space in examples)
#' old <- options(SyncER.wd = getOption("SyncER.wd"))
#' syncer_setup(wd = tempdir())
#' options(old) # restore the previous setting
#'
#' @export
syncer_setup <- function(wd = getOption("SyncER.wd"),
                         verbose = TRUE) {
  if (is.null(wd)) {
    wd <- if (interactive()) {
      readline(prompt = "SyncER: enter the directory to use (leave blank to use the current working directory): ")
    } else {
      ""
    }
  }

  if (!nzchar(wd)) wd <- getwd()
  if (!dir.exists(wd)) {
    dir.create(wd, recursive = TRUE)
    if (verbose) message("Created directory: ", wd)
  }

  # Only remember the path; do not call setwd(), which would change the user's session state.
  options(SyncER.wd = normalizePath(wd, winslash = "/"))
  if (verbose) message("SyncER directory: ", getOption("SyncER.wd"))

  out_dir <- syncer_output_dir()
  if (verbose) message("SyncER output will be saved to: ", out_dir)
  invisible(out_dir)
}

#' Base directory for SyncER input/output folders
#'
#' Returns the directory chosen with \code{\link{syncer_setup}}, or a per-session temporary
#' directory if none was chosen, so that \emph{SyncER} never writes to the user's file space
#' unless asked to.
#'
#' @return Character string: path to the base directory.
#'
#' @keywords internal
syncer_base_dir <- function() {
  getOption("SyncER.wd", default = tempdir())
}

#' Directory to read SyncER input folders from
#'
#' Returns the directory chosen with \code{\link{syncer_setup}}, or the current working
#' directory if none was chosen. Only used for \emph{reading}.
#'
#' @return Character string: path to the input directory.
#'
#' @keywords internal
syncer_input_dir <- function() {
  getOption("SyncER.wd", default = ".")
}

#' Set up SyncER output directory
#'
#' Returns the path to the \code{SyncER_outputs} folder, creating it if it does not already exist.
#' Every \emph{SyncER} function that writes files (CSVs, PDFs) defaults its output-location argument
#' to this folder, so that all package output ends up in one predictable place unless the user
#' explicitly overrides it. If \code{\link{syncer_setup}} has been called earlier in the session, that
#' chosen directory is used; otherwise this defaults to a subfolder of \code{tempdir()}, so that
#' \emph{SyncER} never writes to the user's file space without being asked to.
#'
#' @return Character string: path to the \code{SyncER_outputs} directory.
#'
#' @keywords internal
syncer_output_dir <- function() {
  base_dir <- syncer_base_dir()
  out_dir <- file.path(base_dir, "SyncER_outputs")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_dir
}
