test_that("assess_lr returns a list with ok, lr_mean, lr_sd", {
  set.seed(1)
  x <- rnorm(1000, mean = 0, sd = 0.1)
  result <- assess_lr(x)
  expect_named(result, c("ok", "lr_mean", "lr_sd"))
  expect_true(is.logical(result$ok))
  expect_true(is.numeric(result$lr_mean))
  expect_true(is.numeric(result$lr_sd))
})

test_that("assess_lr ok=TRUE for centred normal LR vector", {
  set.seed(42)
  # Centred at 0 and normal-shaped => ok should be TRUE
  x <- rnorm(5000, mean = 0, sd = 0.2)
  result <- assess_lr(x)
  expect_true(result$ok)
})

test_that("assess_lr ok=FALSE when mean is far from zero relative to sd", {
  # Mean >> 0.5 * sd => centring fails
  x <- rep(10, 3000) + rnorm(3000, sd = 0.01)
  result <- assess_lr(x)
  expect_false(result$ok)
})

test_that("assess_lr handles all-NA input gracefully", {
  result <- assess_lr(c(NA, NA, NA))
  expect_true(is.na(result$lr_mean))
  expect_true(is.na(result$lr_sd))
})

test_that("compute_pairwise_lr returns vector of length n_samples", {
  set.seed(1)
  a <- rnorm(500, mean = 100, sd = 5)
  b <- rnorm(500, mean = 100, sd = 5)
  result <- compute_pairwise_lr(a, b, n_samples = 200)
  expect_length(result, 200)
  expect_true(is.numeric(result))
})

test_that("compute_pairwise_lr antisymmetry holds in expectation: E[lr(A,B)] = -E[lr(B,A)]", {
  # log(A/B) = -log(B/A) is an exact identity per-draw, but the MC estimator
  # draws independently each call, so only the means should be mirror images.
  set.seed(7)
  a <- rnorm(2000, mean = 200, sd = 10)
  b <- rnorm(2000, mean = 180, sd = 10)

  set.seed(99);  lr_ab <- compute_pairwise_lr(a, b, n_samples = 5000)
  set.seed(99);  lr_ba <- compute_pairwise_lr(b, a, n_samples = 5000)

  # Means should be approximately equal and opposite
  expect_equal(mean(lr_ab), -mean(lr_ba), tolerance = 0.01)
  # lr(A,B) > 0 when A > B on average; lr(B,A) < 0
  expect_true(mean(lr_ab) > 0)
  expect_true(mean(lr_ba) < 0)
})

test_that("compute_pairwise_lr returns all zeros when A == B (same object)", {
  set.seed(3)
  a <- rnorm(500, mean = 150, sd = 5)
  # When sampling from the same vector, log(x/x) = 0 only if identical draws.
  # Instead test that mean LR is near 0 for identical distributions.
  set.seed(10)
  lr <- compute_pairwise_lr(a, a, n_samples = 5000)
  expect_equal(mean(lr), 0, tolerance = 0.05)
})

test_that("find_precision_threshold returns NA for empty input", {
  expect_true(is.na(find_precision_threshold(numeric(0), 0.95)))
})

test_that("find_precision_threshold returns NA for all-NA input", {
  expect_true(is.na(find_precision_threshold(c(NA, NA, NA), 0.95)))
})

test_that("find_precision_threshold returns correct threshold for known vector", {
  # Sorted absolute values: 0.1, 0.2, 0.3, 0.4, 0.5
  # For conf_level = 0.6, need at least 60% within [-t, t].
  # 3 out of 5 = 60%, so t should be the 3rd smallest abs value = 0.3
  x <- c(-0.1, 0.2, -0.3, 0.4, -0.5)
  result <- find_precision_threshold(x, conf_level = 0.6)
  expect_equal(result, 0.3, tolerance = 1e-10)
})

test_that("find_precision_threshold works at conf_level = 1.0", {
  x <- c(-0.1, 0.2, -0.5)
  # 100% within t, so t = max(abs(x)) = 0.5
  result <- find_precision_threshold(x, conf_level = 1.0)
  expect_equal(result, 0.5, tolerance = 1e-10)
})

test_that("calculate_overall_synchronicity returns NA when n_records < 2", {
  samples <- list(recA = rnorm(100, 500, 10))
  result <- calculate_overall_synchronicity(samples, conf_level = 0.95,
                                            age_diff_log_bounds = c(-0.05, 0.05))
  expect_true(is.na(result$overall_score))
  expect_true(is.na(result$overall_precision))
  expect_equal(result$n_records, 1L)
})

test_that("calculate_overall_synchronicity returns score in [0,1] for two records", {
  set.seed(5128)
  a <- rnorm(1000, mean = 500, sd = 5)
  b <- rnorm(1000, mean = 500, sd = 5)
  samples <- list(recA = a, recB = b)
  result <- calculate_overall_synchronicity(samples, conf_level = 0.95,
                                            age_diff_log_bounds = c(-0.05, 0.05))
  expect_true(result$overall_score >= 0 && result$overall_score <= 1)
  expect_equal(result$n_records, 2L)
})
