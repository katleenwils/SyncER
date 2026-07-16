# Method resolution in compute_synchronized_ages(): which method each horizon gets.
#
# The methods are told apart by their error: "mean_fixederror" returns the supplied
# age_error verbatim, while "mean" derives the error from the data. Only $processed
# is needed here ($summaries is used by "ageofrecord" only).

make_event_stats <- function(seed = 42, n = 400) {
  set.seed(seed)
  mk <- function(a1, a2) data.frame(isochron1 = rnorm(n, a1, 20),
                                    isochron2 = rnorm(n, a2, 30))
  list(processed = list(core1 = mk(1000, 2000),
                        core2 = mk(1010, 2020),
                        core3 = mk(1005, 1990)))
}

# Runs compute_synchronized_ages() quietly and returns its adjusted_ages.
sync_quietly <- function(...) {
  res <- NULL
  invisible(capture.output(res <- suppressWarnings(compute_synchronized_ages(...))))
  res$adjusted_ages
}

test_that("a named method vector only applies to the horizons it names; the rest fall back to mean", {
  # Regression: a length-1 named vector used to leak its method to every horizon,
  # because the "single method applies to all" branch did not check for names.
  adj <- sync_quietly(
    make_event_stats(),
    method    = c(isochron1 = "mean_fixederror"),
    age_error = c(isochron1 = 5),
    horizons  = c("isochron1", "isochron2"),
    offset    = 0
  )

  expect_true(all(c("isochron1", "isochron2") %in% names(adj)))

  # isochron1 is named -> mean_fixederror -> the fixed error, verbatim
  expect_equal(unique(adj$isochron1$adjusted_error), 5)

  # isochron2 is not named -> "mean" -> error derived from the data, not the fixed 5
  expect_false(isTRUE(all.equal(unique(adj$isochron2$adjusted_error), 5)))
  expect_true(all(adj$isochron2$adjusted_error > 0))
})

test_that("a named method vector with several entries still falls back to mean for unnamed horizons", {
  adj <- sync_quietly(
    make_event_stats(),
    method    = c(isochron1 = "mean_fixederror", isochron3 = "mean_fixederror"),
    age_error = c(isochron1 = 5, isochron3 = 7),
    horizons  = c("isochron1", "isochron2"),
    offset    = 0
  )

  expect_equal(unique(adj$isochron1$adjusted_error), 5)
  expect_false(isTRUE(all.equal(unique(adj$isochron2$adjusted_error), 5)))
})

test_that("a single unnamed method applies to every horizon", {
  adj <- sync_quietly(
    make_event_stats(),
    method    = "mean_fixederror",
    age_error = 5,
    horizons  = c("isochron1", "isochron2"),
    offset    = 0
  )

  expect_equal(unique(adj$isochron1$adjusted_error), 5)
  expect_equal(unique(adj$isochron2$adjusted_error), 5)
})

test_that("a single unnamed method supplied via c() applies to every horizon", {
  # The form used in the vignette: matching_method <- c("Bayesian")
  adj <- sync_quietly(
    make_event_stats(),
    method    = c("mean_fixederror"),
    age_error = 5,
    horizons  = c("isochron1", "isochron2"),
    offset    = 0
  )

  expect_equal(unique(adj$isochron1$adjusted_error), 5)
  expect_equal(unique(adj$isochron2$adjusted_error), 5)
})

test_that("NULL method synchronizes every horizon with mean", {
  es  <- make_event_stats()
  adj <- sync_quietly(es, method = NULL, horizons = c("isochron1", "isochron2"), offset = 0)

  expect_true(all(c("isochron1", "isochron2") %in% names(adj)))
  expect_true(all(adj$isochron1$adjusted_error > 0))
  expect_true(all(adj$isochron2$adjusted_error > 0))

  # "mean" centres on the pooled mean of all records' samples (offset = 0 here)
  expect_equal(unique(adj$isochron1$adjusted_age),
               mean(unlist(lapply(es$processed, `[[`, "isochron1"))))
})


# nonsynchro_horizons matches horizon names exactly: name variants belonging to one
# event are declared via horizon_groups, never inferred from a shared name prefix.

make_numbered_event_stats <- function(seed = 7, n = 300) {
  set.seed(seed)
  mk <- function() data.frame(isochron1  = rnorm(n, 1000, 20),
                              isochron2  = rnorm(n, 2000, 20),
                              isochron10 = rnorm(n, 3000, 20))
  list(processed = list(core1 = mk(), core2 = mk(), core3 = mk()))
}

test_that("a global exclusion matches exactly and does not drop same-prefix horizons", {
  # Regression: excluding "isochron1" used to prefix-match and silently also drop
  # "isochron10" (and isochron11, isochron12, ...).
  adj <- sync_quietly(
    make_numbered_event_stats(),
    method              = NULL,
    nonsynchro_horizons = list(global = "isochron1"),
    offset              = 0
  )

  expect_false("isochron1" %in% names(adj))    # excluded, as asked
  expect_true("isochron10" %in% names(adj))    # must survive: a different horizon
  expect_true("isochron2"  %in% names(adj))
})

test_that("a per-record exclusion matches exactly and does not drop same-prefix horizons", {
  adj <- sync_quietly(
    make_numbered_event_stats(),
    method              = NULL,
    nonsynchro_horizons = list(core2 = "isochron1"),
    offset              = 0
  )

  # isochron1 still synchronized from the other records, but core2 is dropped from it
  expect_true("isochron1" %in% names(adj))
  expect_false("core2" %in% adj$isochron1$record)
  expect_setequal(adj$isochron1$record, c("core1", "core3"))

  # core2's other horizons are untouched, isochron10 included
  expect_true("core2" %in% adj$isochron10$record)
})
