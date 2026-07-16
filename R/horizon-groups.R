#' Derive event-type variables from a defined horizon group
#'
#' Expands a single \code{horizon_groups} definition into the \code{event_types},
#' \code{isochrons}, \code{test_events}, \code{isochron_groups}, and
#' \code{test_horizon_groups} values consumed by the rest of the \emph{SyncER}
#' workflow (\code{load_event_ages()}, \code{compute_synchronicity_values()},
#' \code{synchronize_ages()}, etc.). Combine with \code{list2env()} to unpack the
#' result into your script in a single line, e.g.
#' \code{list2env(load_horizon_names(horizon_groups), environment())}.
#'
#' @param horizon_groups Named list, one entry per event type present in your
#'   records. Each entry is itself a list with:
#'   \itemize{
#'     \item \code{role}: one of \code{"isochron"}, \code{"test"}, or \code{"other"}.
#'     \item \code{members}: optional character vector of the actual per-record
#'           labels that belong to this group (defaults to the entry's own name
#'           when omitted, i.e. a standalone horizon). Only needed when several
#'           differently-named labels in your records should be compared
#'           together as one group (e.g. a deliberately incorrect
#'           \code{"-wrong"} variant).
#'   }
#'
#' @return A named list with elements \code{event_types}, \code{isochrons},
#'   \code{test_events}, \code{isochron_groups}, and \code{test_horizon_groups}.
#'
#' @export
load_horizon_names <- function(horizon_groups) {

  event_types <- names(horizon_groups)

  isochrons   <- character(0)
  test_events <- character(0)
  isochron_groups     <- list()
  test_horizon_groups <- list()

  for (name in names(horizon_groups)) {
    role    <- horizon_groups[[name]]$role
    members <- horizon_groups[[name]]$members
    if (is.null(members)) members <- name # standalone horizon: its only member is itself

    if (role == "isochron") {
      isochrons <- c(isochrons, name)
      isochron_groups[[name]] <- members
    } else if (role == "test") {
      test_events <- c(test_events, name)
      test_horizon_groups[[name]] <- members
    }
  }

  list(
    event_types         = event_types,
    isochrons           = isochrons,
    test_events         = test_events,
    isochron_groups     = isochron_groups,
    test_horizon_groups = test_horizon_groups
  )
}
