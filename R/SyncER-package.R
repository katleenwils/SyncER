#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats cov density dnorm ks.test mahalanobis median pnorm qnorm sd setNames
#' @importFrom graphics hist legend lines
#' @importFrom grDevices dev.off pdf rainbow
#' @importFrom utils modifyList write.csv
## usethis namespace: end
NULL

# Column names referenced via non-standard evaluation (dplyr/ggplot2 pipelines
# and rbacon/rplum-format data frames) rather than as actual global variables.
# Declared here purely to silence R CMD check's "no visible binding for global
# variable" NOTE; these are never assigned or read as real globals.
utils::globalVariables(c(
  "event", "depth", "C14_age", "C14_error", "cc", "labID", "age", "error",
  "d.R", "d.STD", "t.a", "t.b", "Pb210", "Pb210_sd", "thickness",
  "Ra226", "Ra226_sd", "settings",
  "depth(cm)", "density(g/cm^3)", "210Pb(Bq/kg)", "sd(210Pb)",
  "thickness(cm)", "226Ra(Bq/kg)", "sd(226Ra)",
  "xmin", "xmax", "ymin", "ymax", "fill", "x", "y"
))