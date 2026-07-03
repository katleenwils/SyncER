test_that("shape_ok returns NA for length < 3", {
  expect_identical(shape_ok(c(1, 2)), NA)
  expect_identical(shape_ok(numeric(0)), NA)
  expect_identical(shape_ok(c(5)), NA)
})

test_that("shape_ok returns NA for constant vector (sd = 0)", {
  expect_identical(shape_ok(c(3, 3, 3, 3, 3)), NA)
})

test_that("shape_ok returns TRUE for a roughly normal distribution", {
  set.seed(42)
  x <- rnorm(5000)
  expect_true(isTRUE(shape_ok(x)))
})

test_that("shape_ok returns FALSE for highly skewed vector", {
  # Exponential distribution is strongly right-skewed
  set.seed(42)
  x <- rexp(2000, rate = 0.1)
  result <- shape_ok(x)
  # Should be FALSE or NA, but not TRUE
  expect_false(isTRUE(result))
})

test_that("shape_ok ignores NAs", {
  set.seed(42)
  x <- c(rnorm(500), NA, NA)
  # Should not error and should process without the NAs
  result <- shape_ok(x)
  expect_true(is.logical(result) || is.na(result))
})

test_that("bp_datum returns an integer", {
  result <- bp_datum()
  expect_true(is.integer(result))
})

test_that("bp_datum equals current year minus 1950", {
  expected <- as.integer(format(Sys.Date(), "%Y")) - 1950L
  expect_equal(bp_datum(), expected)
})
