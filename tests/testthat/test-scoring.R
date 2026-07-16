test_that("get_horizon_thresholds: scalar confidence_level used for any horizon", {
  summaries <- list(
    recA = data.frame(h1_mean = 500),
    recB = data.frame(h1_mean = 510)
  )
  result <- get_horizon_thresholds("h1", confidence_level = 0.95,
                                       age_difference = 0.05, summaries = summaries)
  expect_equal(result$conf_level_h, 0.95)
})

test_that("get_horizon_thresholds: named confidence_level lookup by horizon", {
  summaries <- list(
    recA = data.frame(h1_mean = 500),
    recB = data.frame(h1_mean = 510)
  )
  cl <- c(h1 = 0.90, h2 = 0.99)
  result <- get_horizon_thresholds("h1", confidence_level = cl,
                                       age_difference = 0.05, summaries = summaries)
  expect_equal(result$conf_level_h, 0.90)
})

test_that("get_horizon_thresholds: named confidence_level falls back to scalar for unknown horizon", {
  summaries <- list(recA = data.frame(h1_mean = 500))
  cl <- c(h1 = 0.90)
  # "unknown_h" is not in names(cl), so it should fall back to cl[1]
  result <- get_horizon_thresholds("unknown_h", confidence_level = cl,
                                       age_difference = 0.05, summaries = summaries)
  expect_equal(result$conf_level_h, 0.90)
})

test_that("get_horizon_thresholds: relative age_diff (< 1) passed through as-is", {
  summaries <- list(recA = data.frame(h1_mean = 500))
  result <- get_horizon_thresholds("h1", confidence_level = 0.95,
                                       age_difference = 0.05, summaries = summaries)
  expect_equal(result$age_diff_h, 0.05)
})

test_that("get_horizon_thresholds: absolute age_diff (>= 1) converted to relative", {
  # Mean age ~500, absolute diff = 50 -> relative = 50/500 = 0.1
  summaries <- list(
    recA = data.frame(h1_mean = 500),
    recB = data.frame(h1_mean = 500)
  )
  result <- get_horizon_thresholds("h1", confidence_level = 0.95,
                                       age_difference = 50, summaries = summaries)
  expect_equal(result$age_diff_h, 50 / 500, tolerance = 1e-6)
})

test_that("get_horizon_thresholds: log_bounds are symmetric around 0", {
  summaries <- list(recA = data.frame(h1_mean = 500))
  result <- get_horizon_thresholds("h1", confidence_level = 0.95,
                                       age_difference = 0.05, summaries = summaries)
  bounds <- result$age_diff_log_bounds
  expect_equal(bounds[1], -bounds[2], tolerance = 1e-10)
  # t = log(1 + 0.05)
  expect_equal(bounds[2], log(1 + 0.05), tolerance = 1e-10)
})

test_that("compute_minimal_precision_threshold: sorted known vector at different confidence levels", {
  x <- c(-0.4, -0.2, 0.1, 0.3, 0.5)
  # abs values sorted: 0.1, 0.2, 0.3, 0.4, 0.5
  # conf = 0.2 -> need >=1/5 within t: first abs_val 0.1 -> proportion = 1/5 = 0.2 >= 0.2 => t=0.1
  expect_equal(compute_minimal_precision_threshold(x, 0.2), 0.1, tolerance = 1e-10)
  # conf = 0.4 -> need >=2/5 within t: abs_val 0.2 -> proportion = 2/5 = 0.4 >= 0.4 => t=0.2
  expect_equal(compute_minimal_precision_threshold(x, 0.4), 0.2, tolerance = 1e-10)
  # conf = 1.0 -> need >=1.0 within t: abs_val 0.5 -> proportion = 5/5 = 1.0 >= 1.0 => t=0.5
  expect_equal(compute_minimal_precision_threshold(x, 1.0), 0.5, tolerance = 1e-10)
})
