# =============================================================================
# SyncER METHOD COMPARISON
#   Cross-method synchronicity evaluation across event deposits:
#   SyncER's Synchronicity Score (SS) vs. OxCal-style agreement/difference,
#   Parnell age differences, and the overlap coefficient.
# =============================================================================
# Standalone reproducibility script accompanying the SyncER manuscript.
#
# It runs ONLY the cross-method comparison across event deposits. It does NOT
# run the SyncER synchronisation pipeline (age-depth model input, tie-point
# synchronisation, or re-running Bacon) -- that pipeline is documented in the
# package vignette (workflow.Rmd). All input here is taken from the installed
# package distribution (inst/extdata): the synthetic record data and the Bacon
# age-depth output for each core, both raw (core*) and synchronised (core*_synced).
#
# SyncER's SS is computed through the package (compute_overall_synchronicity);
# the alternative methods live in other_tests.R (sourced below).
#
# To regenerate the synthetic dataset from scratch, see generate_dataset.R.
# =============================================================================

# ── Dependencies ──────────────────────────────────────────────────────────────
library(SyncER)
library(overlapping)
library(dplyr)
library(tidyr)

# Alternative-method implementations that ship alongside this script.
source(system.file("analysis", "other_tests.R", package = "SyncER"))

# Small helper used by the confusion-matrix printer.
`%||%` <- function(a, b) if (!is.null(a) && length(a) == 1 && !is.na(a)) a else b

# ── Paths ─────────────────────────────────────────────────────────────────────
# Input comes from the installed package; results go to a writable directory.
# Change results_dir if you want the CSVs somewhere permanent.
extdata     <- system.file("extdata", package = "SyncER")
results_dir <- file.path(tempdir(), "SyncER_method_comparison")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
cat("Reading bundled input from:\n  ", extdata,
    "\nWriting results to:\n  ", results_dir, "\n\n")

# ── Input data (from the package distribution) ────────────────────────────────
input       <- read_record_data(folder_path = extdata, file_name = "record_data_input")
record_data <- input$record_data
max_depths  <- input$max_depths

# NOTE ON THE INPUT DATASET
# --------------------------
# This comparison expects the *generic* benchmark produced by generate_dataset.R:
# every deposit is labelled either "synchronous" or "non-synchronous" (plus the
# "isochron<n>" tie points and radiocarbon "sample"s). It must NOT contain any
# "synchro-test" / "synchro-test-wrong" horizons.
# If your record_data still has synchro-test rows, relabel them back to the
# generic "synchronous" / "non-synchronous" classes before running this script.
#
# One definition of every event-type role; load_horizon_names() derives
# event_types / isochrons / test_events / isochron_groups / test_horizon_groups.
horizon_groups <- list(
  "isochron1" = list(role = "isochron"), "isochron2" = list(role = "isochron"),
  "isochron3" = list(role = "isochron"), "isochron4" = list(role = "isochron"),
  "isochron5" = list(role = "isochron"), "isochron6" = list(role = "isochron"),
  "isochron7" = list(role = "isochron"), "isochron8" = list(role = "isochron"),
  "isochron9" = list(role = "isochron"),
  "synchronous"     = list(role = "other"),
  "non-synchronous" = list(role = "other")
)
invisible(list2env(load_horizon_names(horizon_groups), environment()))

# ── Comparison settings ───────────────────────────────────────────────────────
confidence_level <- 0.95
age_difference   <- 0.07
age_offset       <- bp_datum()             # datum shifting ages positive for the SS log-ratio
datums_to_test   <- seq(0, 1000, by = 100) # datums at which SS is scored (ss_pass_<d> columns)

# ── Event ages from the bundled Bacon output ──────────────────────────────────
# No Bacon re-run is needed: the raw (core*) and synchronised (core*_synced)
# model output is included in inst/extdata, so load_event_ages() reads it directly.
event_ages <- load_event_ages(
  folder_path = extdata, record_data = record_data, event_types = event_types,
  max_depths = max_depths, isochrons = isochrons, test_horizons = test_events,
  reload_existing = FALSE, output_dir = results_dir
)
event_ages_synced <- load_event_ages(
  folder_path = extdata, record_data = record_data, event_types = event_types,
  max_depths = max_depths, isochrons = isochrons, test_horizons = test_events,
  synced = "_synced", reload_existing = FALSE, output_dir = results_dir
)


# =============================================================================
# COMPARISON MACHINERY  (verbatim from method_evaluation.R, steps 7-9)
# =============================================================================

# ── compute_ss_score ──────────────────────────────────────────────────────────
# SyncER Synchronicity Score (SS) for a set of posterior age vectors, computed by
# the PACKAGE itself via compute_overall_synchronicity() -- the same function
# compute_synchronicity_values() uses for its overall score -- so the confusion-matrix
# evaluation reflects the package's actual scoring, not a reimplementation.
#
# The overall score is the proportion of joint Monte Carlo draws in which every
# record simultaneously falls within the relative tolerance
# +/- t, t = log(1 + age_difference) (matching SyncER's get_horizon_thresholds()).
#
# Arguments:
#   age_list         : named list of numeric age vectors (one per core)
#   age_difference   : relative tolerance, e.g. 0.05 = 5 %
#   confidence_level : SS value at or above which the event is called synchronous
#   n_samples        : Monte Carlo samples used by the package function
#   offset           : datum added to the ages so the log-ratios are defined
#                      (same convention as process_event_ages(offset = ...));
#                      defaults to the pipeline's age_offset.
#   seed             : fixed RNG seed for the package's internal resampling. Held
#                      constant across every datum (and every variant) so the ONLY
#                      thing that moves the score is the offset, not Monte-Carlo
#                      resampling. The caller's RNG state is restored on exit.
#
# Returns:
#   overall_ss  : package overall synchronicity score (NA if < 2 cores)
compute_ss_score <- function(age_list,
                             age_difference   = 0.07,
                             confidence_level = 0.95,
                             n_samples        = 10000,
                             offset           = age_offset,
                             seed             = 20240101L) {
  if (length(age_list) < 2)
    return(list(overall_ss = NA_real_))

  # Shift onto the datum (positive ages for the log-ratio), then score with the package.
  shifted <- lapply(age_list, function(a) a + offset)
  t <- log(1 + age_difference)   # symmetric log-space bounds, as in the package

  # Deterministic Monte Carlo: fix the RNG so compute_overall_synchronicity()
  # resamples the SAME posterior indices at every datum and for every variant.
  # This makes the per-datum SS scores (and hence the max-wins ranking) a pure
  # function of the offset — re-running gives identical PASS/FAIL, and any change
  # across datums is genuinely the datum, not resampling noise. The previous RNG
  # state is saved and restored on exit so the rest of the pipeline is unaffected.
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = .GlobalEnv)), add = TRUE)
  }
  set.seed(seed)

  res <- tryCatch(
    compute_overall_synchronicity(shifted,
                                    conf_level          = confidence_level,
                                    age_diff_log_bounds = c(-t, t),
                                    n_samples           = n_samples),
    error = function(e) list(overall_score = NA_real_)
  )

  list(overall_ss = res$overall_score)
}


# ── run_all_tests ─────────────────────────────────────────────────────────────
# Runs the full battery of synchronicity tests on a named list of age vectors.
#
# Tests performed:
#   • OxCal Combine (overall + pairwise)   — chi-squared agreement index A_comb
#   • OxCal pairwise differences           — significance of age gaps
#   • Parnell pairwise differences         — non-parametric age gap test
#   • Overlap coefficient OV (overall + pairwise)
#   • Synchronicity Score SS               — log-ratio Monte Carlo method
#
# Arguments:
#   age_list         : named list [core_name → age vector]
#   label            : string label used in warning messages
#   age_difference   : relative tolerance for SS
#   confidence_level : SS pass threshold
#   ss_datums        : optional numeric vector of datums (offsets). When supplied,
#                      the SS overall score is also computed at each datum and
#                      returned in $ss_by_datum (named "ss_pass_<d>"). Used by the
#                      combined path to build the per-datum SS methods; the
#                      pairwise path leaves it NULL (SS is not a pairwise method).
#   seed             : RNG seed handed to compute_ss_score, held constant across
#                      every datum so within-run SS variability is due to the datum
#                      only. Vary it between runs to probe Monte-Carlo variability.
#
# Returns a named list with all raw test results (see fields below).
run_all_tests <- function(age_list, label,
                          age_difference   = 0.07,
                          confidence_level = 0.95,
                          ss_datums        = NULL,
                          seed             = 20240101L) {
  n <- length(age_list)
  if (n < 2) {
    warning(paste("Skipping", label, "- fewer than 2 records"))
    return(NULL)
  }

  core_pairs <- combn(names(age_list), 2, simplify = FALSE)

  combine_overall <- tryCatch(
    oxcal_combine_multiple_pdfs(age_list),
    error = function(e) { warning(paste("OxCal Combine failed:", label, e$message)); NULL }
  )

  pair_combines <- list()
  for (pair in core_pairs) {
    pn <- paste(pair[1], "vs", pair[2])
    pair_combines[[pn]] <- tryCatch(
      oxcal_combine_multiple_pdfs(age_list[pair]), error = function(e) NULL
    )
  }

  oxcal_diffs <- tryCatch(
    calculate_pairwise_differences(age_list),
    error = function(e) { warning(paste("OxCal Diff failed:", label, e$message)); list() }
  )

  parnell_diffs <- tryCatch(
    parnell_age_differences(age_list),
    error = function(e) { warning(paste("Parnell failed:", label, e$message)); list() }
  )

  ov_overall <- tryCatch(overlap(age_list), error = function(e) NULL)

  pair_overlaps <- list()
  for (pair in core_pairs) {
    pn <- paste(pair[1], "vs", pair[2])
    pair_overlaps[[pn]] <- tryCatch(overlap(age_list[pair]), error = function(e) NULL)
  }

  # SS is scored ONLY across the datum sweep below (ss_by_datum). There is no
  # single "default datum" SS score anymore — every datum is reported explicitly.

  # Score SS across a sweep of datums (offsets). Each is an independent package
  # call with the datum applied as the age shift; NA on failure.
  ss_by_datum <- NULL
  if (!is.null(ss_datums)) {
    ss_by_datum <- vapply(ss_datums, function(d) {
      tryCatch(
        compute_ss_score(age_list,
                         age_difference   = age_difference,
                         confidence_level = confidence_level,
                         offset           = d,
                         seed             = seed)$overall_ss,
        error = function(e) NA_real_
      )
    }, numeric(1))
    names(ss_by_datum) <- paste0("ss_pass_", ss_datums)
  }

  list(
    label            = label,
    n_cores          = n,
    combine_overall  = combine_overall,
    pair_combines    = pair_combines,
    oxcal_diffs      = oxcal_diffs,
    parnell_diffs    = parnell_diffs,
    ov_overall       = ov_overall,
    pair_overlaps    = pair_overlaps,
    ss_by_datum      = ss_by_datum
  )
}


# ── build_combined_rows ───────────────────────────────────────────────────────
# Converts a group of variant results (one correct + N wrong) into a data frame
# with one decision row per variant.  This implements the MAX-WINS rule:
#
#   SS       → only the variant with the highest SS score gets TRUE (max-wins,
#               no absolute threshold — the best guess is what matters).
#
#   A_comb   → two-step:
#               1. Find the variant with the maximum A_comb across the group.
#               2. That variant gets TRUE only if its score also clears the
#                  chi-squared–based threshold for the number of cores being
#                  compared (computed via compute_dynamic_thresholds).
#               All other variants get FALSE.
#
#   OV       → same two-step logic as A_comb, using its own dynamic threshold.
#
#   OxcalDiff / Parnell → decided per-variant independently (N_significant == 0).
#
# Why MAX-WINS?
#   We are evaluating whether the method correctly picks the BEST combination,
#   not whether it exceeds a threshold in isolation. A method that gives every
#   combination a passing score is not useful; we want it to rank the correct
#   combination highest.
#
# Arguments:
#   correct_result   : run_all_tests() output for the all-correct variant
#   wrong_results    : named list of run_all_tests() outputs for wrong variants
#   confidence_level : passed for reference; not used for TRUE/FALSE here
#   alpha            : significance level for dynamic threshold (default 0.05)
#
# Returns: data.frame with columns comparison, real_sync,
#          Passes_ChiSq, OV_Passes, Synchro_OxcalDiff_PASS,
#          Synchro_ParnellDiff_PASS, raw_Acomb, raw_OV, plus per-datum SS
#          columns ss_pass_<d> (max-wins verdict) and raw_ss_<d> (raw score).
build_combined_rows <- function(correct_result, wrong_results,
                                confidence_level = 0.95) {
  if (is.null(correct_result)) return(NULL)

  all_results    <- c(list(correct_result), wrong_results)
  real_sync_vals <- vapply(all_results, function(r) isTRUE(r$real_sync), logical(1))

  # ── Threshold: read from the combine result object ───────────────────────────
  # oxcal_combine_multiple_pdfs() already returns threshold_chi_sq = 100/sqrt(2*n).
  # This is internally consistent with A_comb = 100*F^(1/sqrt(n)): F_overall is a
  # product of n values in [0,1] so it shrinks with n even for synchronous events,
  # and the threshold must decrease proportionally to stay fair.
  # We read it from the correct_result (all variants share the same n, so the
  # threshold is identical for every row in this group).
  acomb_thresh <- if (!is.null(correct_result$combine_overall))
    correct_result$combine_overall$threshold_chi_sq
  else NA_real_
  # OV has no analogous formula in the overlapping package, so we keep 0.6 as
  # the conventional baseline (same as threshold_standard for A_comb).
  ov_thresh <- 0.5

  # Helper: safely extract a scalar score
  get_score <- function(res, field) {
    val <- switch(field,
                  "Acomb" = if (!is.null(res$combine_overall)) res$combine_overall$A_comb else NULL,
                  "OV"    = if (!is.null(res$ov_overall))      res$ov_overall$OV          else NULL
    )
    if (is.null(val) || !is.finite(val)) NA_real_ else val
  }

  acomb_vals <- sapply(all_results, get_score, "Acomb")
  ov_vals    <- sapply(all_results, get_score, "OV")

  # SS is evaluated per-datum only (see below); there is no single overall SS_Pass.

  # Per-datum SS: max-wins rule applied to each datum's scores. Column names
  # come from ss_by_datum ("ss_pass_0", "ss_pass_10", ...); all variants in a group
  # were scored on the same datum vector, so the names are consistent.
  datum_names <- NULL
  for (r in all_results)
    if (!is.null(r$ss_by_datum)) { datum_names <- names(r$ss_by_datum); break }
  ss_datum_pass <- list()
  ss_datum_raw  <- list()   # raw SS score per variant, per datum (for verifying PASS/FAIL)
  for (dn in datum_names) {
    dvals <- vapply(all_results, function(r) {
      v <- if (!is.null(r$ss_by_datum)) r$ss_by_datum[[dn]] else NA_real_
      if (is.null(v) || !is.finite(v)) NA_real_ else v
    }, numeric(1))
    dmax <- if (all(is.na(dvals))) NA_real_ else max(dvals, na.rm = TRUE)
    # Max-wins with two guards:
    #  - if the highest score is 0 (or all NA) there is no real winner -> everyone FALSE;
    #  - if several variants tie at the (positive) maximum, they ALL get TRUE.
    ss_datum_pass[[dn]] <- !is.na(dvals) & !is.na(dmax) & dmax > 0 & dvals == dmax
    # Companion raw-score column, named raw_ss_<d> alongside the pass column ss_pass_<d>.
    ss_datum_raw[[sub("^ss_pass_", "raw_ss_", dn)]] <- dvals
  }

  # A_comb: max-wins AND must clear the dynamic chi-squared threshold
  max_Acomb    <- if (all(is.na(acomb_vals))) NA_real_ else max(acomb_vals, na.rm = TRUE)
  best_acomb_ok <- !is.na(max_Acomb) && !is.na(acomb_thresh) && max_Acomb >= acomb_thresh
  is_max_acomb  <- !is.na(acomb_vals) & acomb_vals == max_Acomb & best_acomb_ok

  # OV: max-wins AND must clear the dynamic overlap threshold
  max_OV    <- if (all(is.na(ov_vals))) NA_real_ else max(ov_vals, na.rm = TRUE)
  best_ov_ok <- !is.na(max_OV) && !is.na(ov_thresh) && max_OV >= ov_thresh
  is_max_ov  <- !is.na(ov_vals) & ov_vals == max_OV & best_ov_ok

  # Build one row per variant
  rows <- lapply(seq_along(all_results), function(k) {
    res <- all_results[[k]]
    n_sig_oxcal   <- if (length(res$oxcal_diffs)  > 0)
      sum(sapply(res$oxcal_diffs,  function(d) isTRUE(d$significant)))
    else NA_integer_
    n_sig_parnell <- if (length(res$parnell_diffs) > 0)
      sum(sapply(res$parnell_diffs, function(d) isTRUE(d$significant)))
    else NA_integer_

    data.frame(
      comparison               = res$label,
      real_sync                = real_sync_vals[k],
      Passes_ChiSq             = is_max_acomb[k],
      OV_Passes                = is_max_ov[k],
      Synchro_OxcalDiff_PASS   = if (!is.na(n_sig_oxcal))   n_sig_oxcal   == 0 else NA,
      Synchro_ParnellDiff_PASS = if (!is.na(n_sig_parnell)) n_sig_parnell == 0 else NA,
      raw_Acomb = acomb_vals[k],
      raw_OV    = ov_vals[k],
      used_acomb_thresh = acomb_thresh,   # 100/sqrt(2*n), decreases with n
      stringsAsFactors  = FALSE
    )
  })

  out <- do.call(rbind, rows)

  # Append the per-datum SS_Pass columns (ss_pass_0, ss_pass_10, ...) together with
  # their companion raw-score columns (raw_ss_0, raw_ss_10, ...), aligned by variant
  # position with `out`. The raw scores let you verify each PASS/FAIL by hand.
  if (length(ss_datum_pass) > 0) {
    datum_df <- as.data.frame(ss_datum_pass, check.names = FALSE)
    out <- cbind(out, datum_df)
  }
  if (length(ss_datum_raw) > 0) {
    raw_df <- as.data.frame(ss_datum_raw, check.names = FALSE)
    out <- cbind(out, raw_df)
  }
  out
}


# ── Confusion matrix helpers ───────────────────────────────────────────────────
make_confusion <- function(predicted, actual) {
  valid  <- !is.na(predicted) & !is.na(actual)
  pred   <- as.logical(predicted[valid])
  act    <- as.logical(actual[valid])
  TP <- sum( pred &  act); TN <- sum(!pred & !act)
  FP <- sum( pred & !act); FN <- sum(!pred &  act)
  N  <- TP + TN + FP + FN
  acc    <- if (N > 0)          (TP + TN) / N         else NA_real_
  prec   <- if ((TP + FP) > 0)  TP / (TP + FP)        else NA_real_
  recall <- if ((TP + FN) > 0)  TP / (TP + FN)        else NA_real_
  f1     <- if (!is.na(prec) && !is.na(recall) && (prec + recall) > 0)
    2 * prec * recall / (prec + recall)      else NA_real_
  list(TP = TP, TN = TN, FP = FP, FN = FN, N = N,
       accuracy = acc, precision = prec, recall = recall, F1 = f1)
}

print_confusion_block <- function(df, methods, title) {
  cat(sprintf("\n=== %s ===\n", title))
  cat(sprintf("  %-30s %5s %5s %5s %5s %5s %7s %7s %7s %7s\n",
              "Method", "TP", "TN", "FP", "FN", "N", "Acc", "Prec", "Rec", "F1"))
  cat("  ", strrep("-", 90), "\n", sep = "")
  rows <- list()
  for (method in methods) {
    if (!method %in% colnames(df)) next
    cm <- make_confusion(df[[method]], df$real_sync)
    cat(sprintf("  %-30s %5d %5d %5d %5d %5d %7.3f %7.3f %7.3f %7.3f\n",
                method, cm$TP, cm$TN, cm$FP, cm$FN, cm$N,
                cm$accuracy  %||% NA_real_,
                cm$precision %||% NA_real_,
                cm$recall    %||% NA_real_,
                cm$F1        %||% NA_real_))
    rows[[method]] <- data.frame(
      Method = method, Sheet = title,
      TP = cm$TP, TN = cm$TN, FP = cm$FP, FN = cm$FN, N = cm$N,
      Accuracy  = round(cm$accuracy  %||% NA_real_, 3),
      Precision = round(cm$precision %||% NA_real_, 3),
      Recall    = round(cm$recall    %||% NA_real_, 3),
      F1        = round(cm$F1        %||% NA_real_, 3),
      stringsAsFactors = FALSE)
  }
  invisible(do.call(rbind, rows))
}

pairwise_methods <- c("Passes_ChiSq", "OV_Passes",
                      "Synchro_OxcalDiff_PASS", "Synchro_ParnellDiff_PASS")
# SS is reported as one method per datum (ss_pass_0, ss_pass_10, ...); there is no
# single collapsed SS_Pass. method_performance thus scores SS at every datum.
combined_methods <- c(paste0("ss_pass_", datums_to_test),
                      "Passes_ChiSq", "OV_Passes",
                      "Synchro_OxcalDiff_PASS", "Synchro_ParnellDiff_PASS")


# =============================================================================
# MASTER COMPARISON FUNCTION  (Steps 7–9)
# =============================================================================
# Builds all test sets, runs Other_tests, extracts decisions, and saves output.
#
# Arguments:
#   ev_ages          : event_ages or event_ages_synced
#   rec_data         : record_data list (named list of tibbles with columns
#                      'depth' and 'event'); used to look up isochron depths
#   suffix           : "" for non-synced, "_synced" for synced
#   age_difference   : relative tolerance for SS
#   confidence_level : SS pass threshold
#   ev_ages_struct   : optional structural template (see note below)
#   seed             : RNG seed for this run. When non-NULL, (a) the whole run is
#                      made reproducible via set.seed(seed) up front, (b) the seed
#                      is threaded into the SS scoring, and (c) "_seed<NNN>" is
#                      appended to every output filename. Run the same comparison
#                      under several seeds to separate the two sources of spread:
#                      datum variability WITHIN one _seed file, and Monte-Carlo
#                      variability BETWEEN _seed files at the same datum.
#
# ev_ages_struct note:
#   For the synced run, event_ages_synced may have core names like "core1_synced"
#   while the column/event structure was defined by the non-synced event_ages.
#   Pass event_ages as ev_ages_struct in that case to ensure horizon discovery
#   uses the correct column names, while actual age values come from ev_ages.
# =============================================================================

run_comparisons <- function(ev_ages, rec_data, suffix = "",
                            age_difference   = 0.07,
                            confidence_level = 0.95,
                            ev_ages_struct   = NULL,
                            seed             = NULL) {

  if (is.null(ev_ages_struct)) ev_ages_struct <- ev_ages

  # Reproducible run + filename tag. compute_ss_score() saves/restores the RNG
  # around its own seeding, so this global seed drives the other (Parnell/overlap)
  # Monte-Carlo consistently without disturbing the datum-isolated SS resampling.
  seed_tag <- if (is.null(seed)) "" else sprintf("_seed%03d", seed)
  run_seed <- if (is.null(seed)) 20240101L else seed
  if (!is.null(seed)) set.seed(seed)

  label_tag <- if (nchar(suffix) == 0) "NON-SYNCED" else "SYNCED"
  cat(sprintf("\n\n%s\n  RUNNING COMPARISONS: %s\n%s\n\n",
              strrep("=", 70), label_tag, strrep("=", 70)))

  all_cores    <- names(ev_ages)
  struct_cores <- names(ev_ages_struct)

  # Build mapping: struct_core → actual core name in ev_ages
  strip_suffix <- function(nm) sub("_synced.*$", "", nm)
  core_map <- setNames(
    vapply(struct_cores, function(sc) {
      match <- all_cores[strip_suffix(all_cores) == sc]
      if (length(match) == 0) match <- all_cores[strip_suffix(all_cores) == strip_suffix(sc)]
      if (length(match) == 0) NA_character_ else match[1]
    }, character(1)),
    struct_cores
  )
  core_map         <- core_map[!is.na(core_map)]
  struct_cores_used <- names(core_map)

  # ── get_ages ──────────────────────────────────────────────────────────────────
  # Look up the age vector for a given struct-core and event name.
  # sc = struct core name (e.g. "core1")
  # ev = event column name (e.g. "synchronous_12.5")
  get_ages <- function(sc, ev) {
    ac <- core_map[[sc]]
    if (is.null(ac) || is.na(ac)) return(NULL)
    if (is.null(ev) || length(ev) != 1 || is.na(ev)) return(NULL)
    ev_ages[[ac]][[ev]]
  }

  get_base         <- function(ev) sub("-wrong$", "", ev)
  is_correct_combo <- function(variants, horizon = NULL) {
    if (!is.null(horizon) && grepl("non-synchronous", horizon, fixed = TRUE)) return(FALSE)
    !any(grepl("-wrong$", variants))
  }

  # ── run_pair ──────────────────────────────────────────────────────────────────
  # Run a two-core pairwise test and return a single decision row.
  run_pair <- function(cA, vA, cB, vB, horizon, real_sync) {
    aA <- get_ages(cA, vA); aB <- get_ages(cB, vB)
    if (is.null(aA) || is.null(aB) || length(aA) == 0 || length(aB) == 0) {
      warning(sprintf("run_pair: NULL/empty ages for %s:%s vs %s:%s", cA, vA, cB, vB))
      return(NULL)
    }
    key      <- paste(cA, "vs", cB)
    age_list <- setNames(list(aA, aB), c(cA, cB))
    res <- tryCatch(
      run_all_tests(age_list, key, age_difference = age_difference,
                    confidence_level = confidence_level, seed = run_seed),
      error = function(e) {
        warning(sprintf("run_pair failed for %s vs %s: %s", cA, cB, e$message))
        NULL
      }
    )
    if (is.null(res)) return(NULL)
    pc    <- res$pair_combines[[key]]
    acomb <- if (!is.null(pc)) pc$A_comb           else NA_real_
    # For a pairwise (n=2) call, threshold_chi_sq = 100/sqrt(4) = 50
    acomb_thresh <- if (!is.null(pc)) pc$threshold_chi_sq else 50
    ov    <- if (!is.null(res$pair_overlaps[[key]])) res$pair_overlaps[[key]]$OV else NA_real_
    od    <- res$oxcal_diffs[[key]]
    pd    <- res$parnell_diffs[[key]]

    data.frame(
      horizon                  = horizon,
      record1                  = cA, variant1 = vA,
      record2                  = cB, variant2 = vB,
      real_sync                = isTRUE(real_sync),
      Passes_ChiSq             = !is.na(acomb) & acomb >= acomb_thresh,
      OV_Passes                = !is.na(ov)    & ov    >= 0.5,
      Synchro_OxcalDiff_PASS   = if (!is.null(od)) !isTRUE(od$significant)  else NA,
      Synchro_ParnellDiff_PASS = if (!is.null(pd)) !isTRUE(pd$significant)  else NA,
      raw_Acomb = acomb, raw_OV = ov,
      used_acomb_thresh = acomb_thresh,
      stringsAsFactors  = FALSE
    )
  }

  # ── run_combined_variant ───────────────────────────────────────────────────
  run_combined_variant <- function(cores, variants, label) {
    age_list <- setNames(
      lapply(seq_along(cores), function(k) get_ages(cores[k], variants[k])),
      cores
    )
    if (any(sapply(age_list, is.null)) || any(sapply(age_list, length) == 0)) {
      warning(sprintf("run_combined_variant: NULL/empty ages for label: %s", label))
      return(NULL)
    }
    tryCatch(
      run_all_tests(age_list, label, age_difference = age_difference,
                    confidence_level = confidence_level,
                    ss_datums = datums_to_test, seed = run_seed),
      error = function(e) {
        warning(sprintf("run_combined_variant failed for '%s': %s", label, e$message))
        NULL
      }
    )
  }

  pw_rows    <- list()
  cmb_groups <- list()

  # ── isochron_depth_map ────────────────────────────────────────────────
  # Built FIRST so it is available for depth-sorting all_event_cols below.
  # Nested lookup: isochron_depth_map[[core]][[event_base]] -> depth (cm).
  # Built from rec_data by finding the row where 'event' matches the base name,
  # then reading the 'depth' value. Only the first matching row is used.
  isochron_depth_map <- list()
  for (core in struct_cores_used) {
    rc_key <- core
    if (!rc_key %in% names(rec_data)) {
      rc_key <- names(rec_data)[sub("_synced.*$", "", names(rec_data)) == core]
      if (length(rc_key) == 0) next
      rc_key <- rc_key[1]
    }
    tbl <- rec_data[[rc_key]]
    if (is.null(tbl) || !all(c("depth", "event") %in% names(tbl))) next
    named_rows <- tbl[!is.na(tbl$event) &
                        is.na(suppressWarnings(as.numeric(tbl$event))), ]
    depth_vec  <- setNames(named_rows$depth, named_rows$event)
    depth_vec  <- depth_vec[!duplicated(names(depth_vec))]
    isochron_depth_map[[core]] <- as.list(depth_vec)
  }

  # Look up depth for a named event in a given core.
  # Returns NA_real_ if the event is not found.
  lookup_named_event_depth <- function(core, ev_name) {
    base_ev <- sub("-wrong$", "", ev_name)
    depth   <- isochron_depth_map[[core]][[base_ev]]
    if (is.null(depth) || length(depth) == 0) NA_real_ else depth[[1]]
  }

  # Diagnostic: confirm the map is populated before any lookups are attempted
  cat("isochron depth map:\n")
  for (core in names(isochron_depth_map)) {
    entries <- isochron_depth_map[[core]]
    named_entries <- entries[grepl("^isochron", names(entries))]
    if (length(named_entries) > 0)
      cat(sprintf("  %s: %s\n", core,
                  paste(sprintf("%s=%.1f", names(named_entries),
                                unlist(named_entries)), collapse = ", ")))
  }
  cat("\n")

  # ── Build event column lists, depth-sorted within each event class ─────────
  # For each core the non-numeric column names are split into three groups,
  # each sorted by depth before being recombined:
  #   synchronous_<d>     : depth parsed from name suffix
  #   non-synchronous_<d> : depth parsed from name suffix
  #   named events        : depth looked up from isochron_depth_map;
  #                         events with no depth entry are appended last
  all_event_cols <- lapply(struct_cores_used, function(core) {
    nms <- names(ev_ages_struct[[core]])
    nms <- nms[is.na(suppressWarnings(as.numeric(nms)))]

    syn_evs   <- grep("^synchronous_",     nms, value = TRUE)
    ns_evs    <- grep("^non-synchronous_", nms, value = TRUE)
    named_evs <- nms[!grepl("^synchronous_|^non-synchronous_", nms)]

    syn_depths <- as.numeric(sub("^synchronous_",     "", syn_evs))
    syn_evs    <- syn_evs[order(syn_depths)]

    ns_depths  <- as.numeric(sub("^non-synchronous_", "", ns_evs))
    ns_evs     <- ns_evs[order(ns_depths)]

    named_depths <- vapply(named_evs,
                           function(ev) lookup_named_event_depth(core, ev),
                           numeric(1))
    has_depth  <- !is.na(named_depths)
    named_evs  <- c(named_evs[has_depth][order(named_depths[has_depth])],
                    named_evs[!has_depth])

    c(syn_evs, ns_evs, named_evs)
  })
  names(all_event_cols) <- struct_cores_used

  # ==========================================================================
  # FIX-2 + FIX-3 : Build synchro_by_core and nonsync_by_core BEFORE any loop
  # so they are available to Part 1 cross-type comparisons and Part 3.
  # Depth order is guaranteed by all_event_cols sorting above.
  # ==========================================================================

  synchro_by_core <- list()
  for (core in struct_cores_used) {
    synchro_by_core[[core]] <-
      grep("^synchronous_", all_event_cols[[core]], value = TRUE)
  }

  nonsync_by_core <- list()
  for (core in struct_cores_used) {
    nonsync_by_core[[core]] <-
      grep("^non-synchronous_", all_event_cols[[core]], value = TRUE)
  }

  # ==========================================================================
  # FIX-1 : Define nearest_synchro and nearest_nonsync BEFORE any loop that
  # calls them (Part 1 cross-type block called them before they were defined).
  # ==========================================================================

  # Return up to 2 synchronous events in a core that are closest in depth
  # to target_depth.
  nearest_synchro <- function(core, target_depth, synchro_by_core) {
    syn_evs <- synchro_by_core[[core]]
    if (length(syn_evs) == 0) return(character(0))
    syn_depths <- suppressWarnings(as.numeric(sub("^synchronous_", "", syn_evs)))
    if (all(is.na(syn_depths))) return(character(0))
    diffs     <- abs(syn_depths - target_depth)
    order_idx <- order(diffs)
    syn_evs[order_idx[seq_len(min(2L, length(order_idx)))]]
  }

  # Return up to 2 non-synchronous events in a core closest in depth.
  nearest_nonsync <- function(core, target_depth) {
    ns_evs <- nonsync_by_core[[core]]
    if (length(ns_evs) == 0) return(character(0))
    ns_depths <- as.numeric(sub("^non-synchronous_", "", ns_evs))
    diffs     <- abs(ns_depths - target_depth)
    order_idx <- order(diffs)
    ns_evs[order_idx[seq_len(min(2L, length(order_idx)))]]
  }

  # Helper: event at positional index idx in a core's sorted synchro list
  idx_to_ev <- function(core, idx) {
    evs <- synchro_by_core[[core]]
    if (idx > length(evs) || idx < 1) NA_character_ else evs[idx]
  }



  # ==========================================================================
  # PART 1: NAMED TEST EVENTS (synchro-test*, isochron*)
  # ==========================================================================

  # Collect named bases and sort by median depth across cores so that Part 1
  # processes isochrons and synchro-tests in stratigraphic (shallow-first)
  # order, consistent with the positional indexing used for synchronous events.
  named_bases <- unique(unlist(lapply(all_event_cols, function(evs) {
    bases <- sapply(evs, get_base)
    bases[!grepl("^synchronous_|^non-synchronous_", bases)]
  })))
  named_base_depths <- vapply(named_bases, function(b) {
    depths <- vapply(struct_cores_used,
                     function(core) lookup_named_event_depth(core, b),
                     numeric(1))
    median(depths, na.rm = TRUE)
  }, numeric(1))
  named_bases <- named_bases[order(named_base_depths)]

  for (base in named_bases) {
    vpc <- list()
    for (core in struct_cores_used) {
      evs     <- all_event_cols[[core]]
      matches <- evs[sapply(evs, get_base) == base]
      if (length(matches) > 0) vpc[[core]] <- matches
    }
    cores_with_base <- names(vpc)
    if (length(cores_with_base) < 2) next

    # ── Pairwise ──
    # Phase 1: all explicit variant combinations from the data (correct + -wrong)
    core_pairs <- combn(cores_with_base, 2, simplify = FALSE)
    for (pair in core_pairs) {
      cA <- pair[1]; cB <- pair[2]
      for (vA in vpc[[cA]]) {
        for (vB in vpc[[cB]]) {
          key <- paste(base, cA, vA, cB, vB, sep = "|")
          row <- run_pair(cA, vA, cB, vB, base, is_correct_combo(c(vA, vB), horizon = base))
          if (!is.null(row)) pw_rows[[key]] <- row
        }
      }
    }

    # Phase 2: cross-type pairwise swaps — for each core pair, hold one core's
    # correct event fixed and substitute the other core's event with its nearest
    # synchronous or non-synchronous event at a similar depth.
    # Mirrors the cross-type pairwise logic used for synchronous events in Part 2.
    correct_variants <- lapply(vpc, function(vs) vs[!grepl("-wrong$", vs)])
    cores_with_correct <- names(Filter(function(v) length(v) > 0, correct_variants))

    if (length(cores_with_correct) >= 2) {
      for (pair in combn(cores_with_correct, 2, simplify = FALSE)) {
        cA <- pair[1]; cB <- pair[2]
        vA_correct <- correct_variants[[cA]][1]
        vB_correct <- correct_variants[[cB]][1]

        dA <- lookup_named_event_depth(cA, vA_correct)
        dB <- lookup_named_event_depth(cB, vB_correct)

        # Swap cB's event for its nearest sync/nonsync, keep cA correct
        if (!is.na(dA)) {
          for (v_syn in nearest_synchro(cB, dA, synchro_by_core)) {
            key <- paste("xtype", base, cA, vA_correct, cB, v_syn, sep = "|")
            if (!key %in% names(pw_rows)) {
              row <- run_pair(cA, vA_correct, cB, v_syn,
                              sprintf("%s_vs_sync", base), FALSE)
              if (!is.null(row)) pw_rows[[key]] <- row
            }
          }
          for (v_ns in nearest_nonsync(cB, dA)) {
            key <- paste("xtype", base, cA, vA_correct, cB, v_ns, sep = "|")
            if (!key %in% names(pw_rows)) {
              row <- run_pair(cA, vA_correct, cB, v_ns,
                              sprintf("%s_vs_nonsync", base), FALSE)
              if (!is.null(row)) pw_rows[[key]] <- row
            }
          }
        }

        # Swap cA's event for its nearest sync/nonsync, keep cB correct
        if (!is.na(dB)) {
          for (v_syn in nearest_synchro(cA, dB, synchro_by_core)) {
            key <- paste("xtype", base, cA, v_syn, cB, vB_correct, sep = "|")
            if (!key %in% names(pw_rows)) {
              row <- run_pair(cA, v_syn, cB, vB_correct,
                              sprintf("%s_vs_sync", base), FALSE)
              if (!is.null(row)) pw_rows[[key]] <- row
            }
          }
          for (v_ns in nearest_nonsync(cA, dB)) {
            key <- paste("xtype", base, cA, v_ns, cB, vB_correct, sep = "|")
            if (!key %in% names(pw_rows)) {
              row <- run_pair(cA, v_ns, cB, vB_correct,
                              sprintf("%s_vs_nonsync", base), FALSE)
              if (!is.null(row)) pw_rows[[key]] <- row
            }
          }
        }
      }
    }

    # ── Combined ──────────────────────────────────────────────────────────────
    # Always run combined testing for named events (isochrons and
    # synchro-tests), mirroring the logic used for synchronous events in Part 2.
    #
    # Phase 1: test every explicit variant combination (correct + any -wrong).
    # Phase 2: extend wrong_subset with cross-type swaps (nearest sync /
    #          non-sync in each core), then re-run build_combined_rows.

    has_wrong <- any(unlist(lapply(vpc, function(vs) any(grepl("-wrong$", vs)))))
    variant_grid <- expand.grid(vpc, stringsAsFactors = FALSE)
    colnames(variant_grid) <- cores_with_base

    group_res <- list(); group_correct <- logical()

    # Phase 1: explicit variant combinations
    for (row_i in seq_len(nrow(variant_grid))) {
      chosen        <- unlist(variant_grid[row_i, ], use.names = FALSE)
      names(chosen) <- cores_with_base
      correct       <- is_correct_combo(chosen, horizon = base)
      label         <- paste(sprintf("%s:%s", cores_with_base, chosen), collapse = " | ")

      res <- run_combined_variant(cores_with_base, chosen, label)
      if (is.null(res)) next
      res$real_sync        <- correct
      group_res[[label]]   <- res
      group_correct[label] <- correct
    }

    # Phase 2: cross-type wrong variants — swap the correct event in each core
    # for its nearest synchronous or non-synchronous event at a similar depth.
    # This matches the cross-type extension used for synchronous events in Part 2.
    #
    # correct_chosen is taken directly from correct_variants (already computed
    # for the pairwise Phase 2 above), avoiding any fragile label-parsing.
    if (length(cores_with_correct) >= 2) {
      correct_chosen <- setNames(
        vapply(cores_with_correct, function(c) correct_variants[[c]][1], character(1)),
        cores_with_correct
      )

      for (core in cores_with_correct) {
        base_ev <- correct_chosen[core]
        depth   <- lookup_named_event_depth(core, base_ev)
        if (is.na(depth)) {
          warning(sprintf("Combined Phase 2: no depth found for event '%s' in core '%s' — skipping cross-type swaps for this core", base_ev, core))
          next
        }

        # Nearest synchronous events as wrong alternatives
        for (v_syn in nearest_synchro(core, depth, synchro_by_core)) {
          swapped       <- correct_chosen; swapped[core] <- v_syn
          label_ns      <- paste(sprintf("%s:%s", cores_with_base, swapped), collapse = " | ")
          if (label_ns %in% names(group_res)) next
          res_ns <- run_combined_variant(cores_with_base, swapped, label_ns)
          if (!is.null(res_ns)) {
            res_ns$real_sync        <- FALSE
            group_res[[label_ns]]   <- res_ns
            group_correct[label_ns] <- FALSE
          }
        }

        # Nearest non-synchronous events as wrong alternatives
        for (v_ns in nearest_nonsync(core, depth)) {
          swapped       <- correct_chosen; swapped[core] <- v_ns
          label_ns      <- paste(sprintf("%s:%s", cores_with_base, swapped), collapse = " | ")
          if (label_ns %in% names(group_res)) next
          res_ns <- run_combined_variant(cores_with_base, swapped, label_ns)
          if (!is.null(res_ns)) {
            res_ns$real_sync        <- FALSE
            group_res[[label_ns]]   <- res_ns
            group_correct[label_ns] <- FALSE
          }
        }
      }
    }

    # Build output rows only when there is at least one wrong variant to compare
    # against (a group with only the correct variant yields no MAX-WINS signal).
    if (length(group_res) > 0 && any(!group_correct)) {
      correct_nm   <- names(group_correct)[group_correct]
      wrong_nm     <- names(group_correct)[!group_correct]
      correct_res  <- if (length(correct_nm) > 0) group_res[[correct_nm[1]]] else group_res[[1]]
      wrong_subset <- group_res[wrong_nm]
      rows <- build_combined_rows(correct_res, wrong_subset,
                                  confidence_level = confidence_level)
      if (!is.null(rows)) { rows$horizon <- base; cmb_groups[[base]] <- rows }
    }

  }   # end for (base in named_bases)


  # ==========================================================================
  # PART 2: GENERIC SYNCHRONOUS EVENTS (synchronous_<depth>)
  # ==========================================================================

  counts_per_core     <- sapply(synchro_by_core, length)
  cores_with_synchro  <- struct_cores_used[counts_per_core > 0]
  max_synchro <- if (length(cores_with_synchro) >= 2) max(counts_per_core[cores_with_synchro]) else 0

  cat(sprintf("Synchronous events per core:\n"))
  for (core in struct_cores_used)
    cat(sprintf("  %s: %d events\n", core, counts_per_core[core]))
  cat(sprintf("Maximum synchronous position across all cores: %d\n\n", max_synchro))

  if (max_synchro > 0) {

    for (idx in seq_len(max_synchro)) {
      horizon     <- sprintf("synchronous_pos%d", idx)
      cores_at_idx <- struct_cores_used[counts_per_core >= idx]
      if (length(cores_at_idx) < 2) next

      correct_evs <- setNames(
        sapply(cores_at_idx, idx_to_ev, idx = idx),
        cores_at_idx
      )

      # ── Pairwise ──
      core_pairs <- combn(cores_at_idx, 2, simplify = FALSE)
      for (pair in core_pairs) {
        cA <- pair[1]; cB <- pair[2]

        key <- paste(horizon, cA, correct_evs[cA], cB, correct_evs[cB], sep = "|")
        row <- run_pair(cA, correct_evs[cA], cB, correct_evs[cB], horizon, TRUE)
        if (!is.null(row)) pw_rows[[key]] <- row

        vA_up <- idx_to_ev(cA, idx + 1)
        if (!is.na(vA_up)) {
          key_wu <- paste(horizon, cA, vA_up, cB, correct_evs[cB], sep = "|")
          row_wu <- run_pair(cA, vA_up, cB, correct_evs[cB], horizon, FALSE)
          if (!is.null(row_wu)) pw_rows[[key_wu]] <- row_wu
        }

        vB_up <- idx_to_ev(cB, idx + 1)
        if (!is.na(vB_up)) {
          key_wu2 <- paste(horizon, cA, correct_evs[cA], cB, vB_up, sep = "|")
          row_wu2 <- run_pair(cA, correct_evs[cA], cB, vB_up, horizon, FALSE)
          if (!is.null(row_wu2)) pw_rows[[key_wu2]] <- row_wu2
        }
      }

      # ── Combined ──
      group_res <- list(); group_correct <- logical()

      label_c <- paste(sprintf("%s:%s", cores_at_idx, correct_evs[cores_at_idx]),
                       collapse = " | ")
      res_c   <- run_combined_variant(cores_at_idx, correct_evs[cores_at_idx], label_c)
      if (!is.null(res_c)) {
        res_c$real_sync        <- TRUE
        group_res[[label_c]]   <- res_c
        group_correct[label_c] <- TRUE
      }

      for (wrong_core in cores_at_idx) {
        for (shift in c(-1L, +1L)) {
          v_shift <- idx_to_ev(wrong_core, idx + shift)
          if (is.na(v_shift)) next
          wrong_evs             <- correct_evs
          wrong_evs[wrong_core] <- v_shift
          label_w <- paste(sprintf("%s:%s", cores_at_idx, wrong_evs[cores_at_idx]),
                           collapse = " | ")
          if (label_w %in% names(group_res)) next
          res_w <- run_combined_variant(cores_at_idx, wrong_evs[cores_at_idx], label_w)
          if (!is.null(res_w)) {
            res_w$real_sync        <- FALSE
            group_res[[label_w]]   <- res_w
            group_correct[label_w] <- FALSE
          }
        }
      }

      if (length(group_res) > 0) {
        correct_nm   <- names(group_correct)[group_correct]
        wrong_nm     <- names(group_correct)[!group_correct]
        correct_res  <- if (length(correct_nm) > 0) group_res[[correct_nm[1]]] else group_res[[1]]
        wrong_subset <- group_res[wrong_nm]
        rows <- build_combined_rows(correct_res, wrong_subset,
                                    confidence_level = confidence_level)
        if (!is.null(rows)) { rows$horizon <- horizon; cmb_groups[[horizon]] <- rows }
      }

      # ── Cross-type: synchronous vs nearest non-synchronous ──
      synchro_depths <- setNames(
        sapply(cores_at_idx, function(core) {
          as.numeric(sub("^synchronous_", "", correct_evs[core]))
        }),
        cores_at_idx
      )

      core_pairs_idx <- combn(cores_at_idx, 2, simplify = FALSE)
      for (pair in core_pairs_idx) {
        cA <- pair[1]; cB <- pair[2]
        vA_sync <- correct_evs[cA]; vB_sync <- correct_evs[cB]

        # Up to 2 nearest non-sync in cB paired against the sync in cA
        for (vB_ns in nearest_nonsync(cB, synchro_depths[cA])) {
          key_ns <- paste("synvsnonsync", horizon, cA, vA_sync, cB, vB_ns, sep = "|")
          if (!key_ns %in% names(pw_rows)) {
            row_ns <- run_pair(cA, vA_sync, cB, vB_ns,
                               sprintf("%s_vs_nonsync", horizon), FALSE)
            if (!is.null(row_ns)) pw_rows[[key_ns]] <- row_ns
          }
        }

        # Up to 2 nearest non-sync in cA paired against the sync in cB
        for (vA_ns in nearest_nonsync(cA, synchro_depths[cB])) {
          key_ns2 <- paste("synvsnonsync", horizon, cA, vA_ns, cB, vB_sync, sep = "|")
          if (!key_ns2 %in% names(pw_rows)) {
            row_ns2 <- run_pair(cA, vA_ns, cB, vB_sync,
                                sprintf("%s_vs_nonsync", horizon), FALSE)
            if (!is.null(row_ns2)) pw_rows[[key_ns2]] <- row_ns2
          }
        }
      }

      # Combined: swap one core's sync for its nearest non-sync (up to 2), then re-evaluate
      for (swap_core in cores_at_idx) {
        ns_candidates <- nearest_nonsync(swap_core, synchro_depths[swap_core])
        for (v_ns in ns_candidates) {   # up to 2 nearest, consistent with pairwise
          if (is.na(v_ns) || v_ns == "") next
          swapped_evs            <- correct_evs
          swapped_evs[swap_core] <- v_ns
          label_ns <- paste(sprintf("%s:%s", cores_at_idx, swapped_evs[cores_at_idx]),
                            collapse = " | ")
          if (label_ns %in% names(group_res)) next
          res_ns <- run_combined_variant(cores_at_idx, swapped_evs[cores_at_idx], label_ns)
          if (!is.null(res_ns)) {
            res_ns$real_sync        <- FALSE
            group_res[[label_ns]]   <- res_ns
            group_correct[label_ns] <- FALSE
          }
        }
      }

      # Re-run build_combined_rows with the extended wrong_subset
      if (length(group_res) > 0) {
        correct_nm   <- names(group_correct)[group_correct]
        wrong_nm     <- names(group_correct)[!group_correct]
        correct_res  <- if (length(correct_nm) > 0) group_res[[correct_nm[1]]] else group_res[[1]]
        wrong_subset <- group_res[wrong_nm]
        rows <- build_combined_rows(correct_res, wrong_subset,
                                    confidence_level = confidence_level)
        if (!is.null(rows)) { rows$horizon <- horizon; cmb_groups[[horizon]] <- rows }
      }

    }   # end for (idx in seq_len(max_synchro))
  }   # end if (max_synchro > 0)


  # ==========================================================================
  # PART 3: GENERIC NON-SYNCHRONOUS EVENTS (non-synchronous_<depth>)
  # ==========================================================================

  ns_counts <- sapply(nonsync_by_core, length)
  cat(sprintf("Non-synchronous events per core:\n"))
  for (core in struct_cores_used)
    cat(sprintf("  %s: %d events\n", core, ns_counts[core]))

  seen_pairs <- character()
  core_pairs <- combn(struct_cores_used, 2, simplify = FALSE)

  for (pair in core_pairs) {
    cA <- pair[1]; cB <- pair[2]
    nA <- ns_counts[cA]; nB <- ns_counts[cB]
    if (nA == 0 || nB == 0) next

    all_iA <- seq_len(nA); all_iB <- seq_len(nB)

    # seq_bounded() returns integer(0) when lo > hi, avoiding R's descending
    # `lo:hi` pitfall that would otherwise yield out-of-range (NA) event indices
    # when the two cores have different numbers of non-synchronous events.
    seq_bounded <- function(lo, hi) if (lo <= hi) lo:hi else integer(0)

    for (iA in all_iA) {
      vA <- nonsync_by_core[[cA]][iA]
      for (iB in seq_bounded(max(1L, iA - 1L), min(nB, iA + 1L))) {
        vB        <- nonsync_by_core[[cB]][iB]
        dedup_key <- paste(sort(c(paste(cA, iA), paste(cB, iB))), collapse = "||")
        if (dedup_key %in% seen_pairs) next
        seen_pairs <- c(seen_pairs, dedup_key)
        horizon    <- sprintf("non-synchronous_%s_pos%d_vs_%s_pos%d", cA, iA, cB, iB)
        row <- run_pair(cA, vA, cB, vB, horizon, FALSE)
        if (!is.null(row)) pw_rows[[paste("nonsync", cA, vA, cB, vB, sep = "|")]] <- row
      }
    }

    for (iB in all_iB) {
      vB <- nonsync_by_core[[cB]][iB]
      for (iA in seq_bounded(max(1L, iB - 1L), min(nA, iB + 1L))) {
        vA        <- nonsync_by_core[[cA]][iA]
        dedup_key <- paste(sort(c(paste(cA, iA), paste(cB, iB))), collapse = "||")
        if (dedup_key %in% seen_pairs) next
        seen_pairs <- c(seen_pairs, dedup_key)
        horizon    <- sprintf("non-synchronous_%s_pos%d_vs_%s_pos%d", cA, iA, cB, iB)
        row <- run_pair(cA, vA, cB, vB, horizon, FALSE)
        if (!is.null(row)) pw_rows[[paste("nonsync", cA, vA, cB, vB, sep = "|")]] <- row
      }
    }
  }


  # ==========================================================================
  # ASSEMBLE, PRINT, AND SAVE
  # ==========================================================================

  pairwise_df <- do.call(rbind, Filter(Negate(is.null), pw_rows))
  rownames(pairwise_df) <- NULL
  pairwise_df$real_sync <- as.logical(pairwise_df$real_sync)
  cat(sprintf("\nPairwise table: %d rows | TRUE: %d | FALSE: %d\n",
              nrow(pairwise_df),
              sum( pairwise_df$real_sync, na.rm = TRUE),
              sum(!pairwise_df$real_sync, na.rm = TRUE)))

  combined_df <- do.call(rbind, Filter(Negate(is.null), cmb_groups))
  rownames(combined_df) <- NULL
  combined_df$real_sync <- as.logical(combined_df$real_sync)
  cat(sprintf("Combined table: %d rows | TRUE: %d | FALSE: %d\n\n",
              nrow(combined_df),
              sum( combined_df$real_sync, na.rm = TRUE),
              sum(!combined_df$real_sync, na.rm = TRUE)))

  pw_summary  <- print_confusion_block(pairwise_df, pairwise_methods,
                                       paste("Pairwise", label_tag))
  cmb_summary <- print_confusion_block(combined_df,  combined_methods,
                                       paste("Combined", label_tag))

  # Stamp the seed on every table so files can be pooled and grouped by seed later.
  perf_df <- rbind(pw_summary, cmb_summary)
  if (!is.null(seed)) {
    pairwise_df$seed <- seed
    combined_df$seed <- seed
    perf_df$seed     <- seed
  }

  pw_file   <- file.path(results_dir, paste0("pairwise_decisions", suffix, seed_tag, ".csv"))
  cmb_file  <- file.path(results_dir, paste0("combined_decisions", suffix, seed_tag, ".csv"))
  perf_file <- file.path(results_dir, paste0("method_performance",  suffix, seed_tag, ".csv"))

  write.csv(pairwise_df, pw_file,   row.names = FALSE)
  write.csv(combined_df, cmb_file,  row.names = FALSE)
  write.csv(perf_df,     perf_file, row.names = FALSE)

  cat(sprintf("\nSaved [%s%s]:\n  %s\n  %s\n  %s\n",
              label_tag, seed_tag, pw_file, cmb_file, perf_file))

  invisible(list(pairwise = pairwise_df, combined = combined_df,
                 performance = perf_df))
}

# =============================================================================
# RUN THE COMPARISONS  (raw and synchronised event ages), ACROSS A SET OF SEEDS
# =============================================================================
# The entire comparison is repeated once per seed. Each seed produces its own
# self-consistent set of CSVs tagged "_seed<NNN>". This separates the two sources
# of spread:
#   • WITHIN one seed: the resampling is fixed, so any change across the per-datum
#     SS columns (ss_pass_0, ss_pass_100, ...) is genuine datum variability.
#   • BETWEEN seeds (same datum column across _seed files): the datum is fixed, so
#     any change is Monte-Carlo / random variability of the SS score.
# Add or remove seeds here; each is an independent full run.
seeds_to_test <- c(1, 2, 3)

all_files     <- character()
all_summaries <- list()   # per-seed performance summaries, pooled for the across-seeds table

for (s in seeds_to_test) {
  cat(sprintf("\n\n##############  SEED %03d  ##############\n", s))

  results_nonsynced <- run_comparisons(event_ages, rec_data = record_data,
                                       suffix = "", seed = s)
  results_synced    <- run_comparisons(event_ages_synced, rec_data = record_data,
                                       suffix = "_synced", ev_ages_struct = event_ages,
                                       seed = s)

  # Per-seed performance summary across all four sheets (pairwise/combined x
  # non-synced/synced), keeping every method row -- including the per-datum SS
  # methods ss_pass_0, ss_pass_100, ...
  seed_tag <- sprintf("_seed%03d", s)
  method_performance_summary <- rbind(results_nonsynced$performance,
                                      results_synced$performance)
  summary_file <- file.path(results_dir, paste0("method_performance_summary", seed_tag, ".csv"))
  write.csv(method_performance_summary, summary_file, row.names = FALSE)

  all_summaries[[seed_tag]] <- method_performance_summary
  all_files <- c(all_files,
                 paste0("pairwise_decisions",        seed_tag, ".csv"),
                 paste0("combined_decisions",        seed_tag, ".csv"),
                 paste0("method_performance",        seed_tag, ".csv"),
                 paste0("pairwise_decisions_synced", seed_tag, ".csv"),
                 paste0("combined_decisions_synced", seed_tag, ".csv"),
                 paste0("method_performance_synced", seed_tag, ".csv"),
                 paste0("method_performance_summary", seed_tag, ".csv"))
}


# =============================================================================
# ACROSS-SEED AGGREGATE
# =============================================================================
# Pool the per-seed performance summaries and, for each (Method, Sheet), report
# the mean and sd of every metric over the seeds. This is the between-seed
# (Monte-Carlo) view: for a given method — e.g. a per-datum SS method ss_pass_500
# — <metric>_mean is its typical performance and <metric>_sd is how much random
# resampling moves it. A large sd flags a method/datum whose PASS/FAIL is unstable.
pooled_summary <- do.call(rbind, all_summaries)
rownames(pooled_summary) <- NULL

metric_cols <- c("TP", "TN", "FP", "FN", "N",
                 "Accuracy", "Precision", "Recall", "F1")

across_seeds <- pooled_summary %>%
  group_by(Method, Sheet) %>%
  summarise(
    n_seeds = n(),
    across(all_of(metric_cols),
           list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))

across_file <- file.path(results_dir, "method_performance_across_seeds.csv")
write.csv(across_seeds, across_file, row.names = FALSE)
all_files <- c(all_files, "method_performance_across_seeds.csv")


cat("\n\n=== DONE ===\nResults written to:\n  ", results_dir, "\n")
cat(sprintf("Seeds: %s\n\nFiles:\n", paste(sprintf("%03d", seeds_to_test), collapse = ", ")))
for (f in all_files) cat("  ", f, "\n")
