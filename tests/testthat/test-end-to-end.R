# End-to-end tests for the synchronicity compute pipeline.
#
# These drive the full statistical core of the method --
#   process_event_ages() -> compute_synchronicity()
# -- on small synthetic posterior age-sample tables where the correct verdict is
# known *by construction*: records drawn from the same age distribution must score
# as synchronous, and a record shifted well outside the age-difference tolerance
# must not.
#
# Bacon age-modelling (the raw-Excel -> MCMC front end) is deliberately NOT run
# here: it is slow, calls an external sampler, and is not reproducible inside
# R CMD check. The synthetic samples below stand in for its output, which is the
# input the method actually consumes.

# Build an `out_data`-style list: one data frame per record, each with posterior
# age samples for the given horizon column(s). Named `record_means` -> record names.
make_out_data <- function(record_means, n = 3000, sd = 30, seed = 1,
                          horizon = "tephra1") {
  set.seed(seed)
  lapply(record_means, function(mu) {
    df <- data.frame(rnorm(n, mean = mu, sd = sd))
    names(df) <- horizon
    df
  })
}

test_that("process_event_ages yields the structure compute_synchronicity expects", {
  out_data    <- make_out_data(c(recA = 5000, recB = 5000, recC = 5000))
  event_stats <- process_event_ages(out_data, event_deposits = "tephra", offset = 0)

  expect_named(event_stats, c("processed", "summaries"))
  expect_equal(names(event_stats$processed), c("recA", "recB", "recC"))
  # Summaries must carry the <horizon>_{min,max,mean,sigma} columns the
  # pairwise comparison and threshold resolution read back out.
  expect_true(all(c("tephra1_min", "tephra1_max", "tephra1_mean", "tephra1_sigma")
                  %in% names(event_stats$summaries$recA)))
})

test_that("end-to-end: three co-located records score as highly synchronous", {
  # All three records share the same age distribution => every pairwise and the
  # joint overall score should be ~1.
  out_data    <- make_out_data(c(recA = 5000, recB = 5000, recC = 5000))
  event_stats <- process_event_ages(out_data, event_deposits = "tephra", offset = 0)

  res <- compute_synchronicity(event_stats, event_names = "tephra1",
                               n_samples = 2000, seed = 5128)

  overall <- res$overall_scores$tephra1$all
  expect_equal(overall$n_records, 3L)
  expect_gt(overall$overall_score, 0.8)

  # Three records => three unordered pairwise comparisons, each scoring high.
  expect_equal(nrow(res$all_scores), 3L)
  expect_true(all(res$all_scores$score > 0.8))

  # Precision, when estimable, is a small positive relative tolerance.
  prec <- overall$overall_precision
  expect_true(is.na(prec) || (prec > 0 && prec < 0.1))
})

test_that("end-to-end: one strongly offset record breaks synchronicity", {
  # recC sits ~1000 yr (20%) away from the others -- far outside the default 5%
  # tolerance -- so the joint (all-records-simultaneously) score collapses to ~0.
  out_data    <- make_out_data(c(recA = 5000, recB = 5000, recC = 6000))
  event_stats <- process_event_ages(out_data, event_deposits = "tephra", offset = 0)

  res <- compute_synchronicity(event_stats, event_names = "tephra1",
                               n_samples = 2000, seed = 5128)

  expect_lt(res$overall_scores$tephra1$all$overall_score, 0.2)

  # The recA-vs-recC and recB-vs-recC pairs fail; recA-vs-recB stays synchronous.
  scores <- setNames(res$all_scores$score, res$all_scores$record_pair)
  expect_gt(scores[["recA_recB_lr"]], 0.8)
  expect_lt(scores[["recA_recC_lr"]], 0.2)
  expect_lt(scores[["recB_recC_lr"]], 0.2)
})

test_that("end-to-end: results are reproducible for a fixed seed", {
  # The reproducibility guarantee the method rests on: identical input + seed
  # must give byte-identical scores across independent runs.
  out_data    <- make_out_data(c(recA = 5000, recB = 5010, recC = 4990))
  event_stats <- process_event_ages(out_data, event_deposits = "tephra", offset = 0)

  res1 <- compute_synchronicity(event_stats, "tephra1", n_samples = 2000, seed = 5128)
  res2 <- compute_synchronicity(event_stats, "tephra1", n_samples = 2000, seed = 5128)

  expect_equal(res1$overall_scores$tephra1$all$overall_score,
               res2$overall_scores$tephra1$all$overall_score)
  expect_equal(res1$all_horizon_stats, res2$all_horizon_stats)

  # A different seed should generally move the Monte Carlo estimate a little.
  res3 <- compute_synchronicity(event_stats, "tephra1", n_samples = 2000, seed = 99)
  expect_false(isTRUE(all.equal(res1$all_horizon_stats$mean_LR,
                                res3$all_horizon_stats$mean_LR)))
})

test_that("end-to-end: horizon grouping compares variant columns across records", {
  # recA carries the event as "tephra1", recB as the variant "tephra1a".
  # Declaring them one group must produce a scored comparison rather than
  # treating them as two unrelated standalone horizons.
  set.seed(7)
  out_data <- list(
    recA = data.frame(tephra1  = rnorm(3000, 5000, 30)),
    recB = data.frame(tephra1a = rnorm(3000, 5000, 30))
  )
  event_stats <- process_event_ages(out_data, event_deposits = "tephra", offset = 0)

  res <- compute_synchronicity(
    event_stats,
    event_names    = "tephra1",
    horizon_groups = list(tephra1 = c("tephra1", "tephra1a")),
    n_samples      = 2000, seed = 5128
  )

  expect_true(length(res$overall_scores$tephra1) > 0)
  expect_gt(nrow(res$all_horizon_stats), 0)
})
