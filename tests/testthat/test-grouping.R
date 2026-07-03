test_that("build_group_lookup returns NULL for NULL input", {
  expect_null(build_group_lookup(NULL))
})

test_that("build_group_lookup maps horizons to group names correctly", {
  groups <- list(
    groupA = c("h1", "h2"),
    groupB = c("h3", "h4")
  )
  result <- build_group_lookup(groups)
  expect_equal(result[["h1"]], "groupA")
  expect_equal(result[["h2"]], "groupA")
  expect_equal(result[["h3"]], "groupB")
  expect_equal(result[["h4"]], "groupB")
})

test_that("build_group_lookup: horizon in multiple groups gets last group assignment", {
  # When a horizon appears in multiple groups, the last one wins (loop overwrites)
  groups <- list(
    groupA = c("h1", "shared"),
    groupB = c("shared", "h2")
  )
  result <- build_group_lookup(groups)
  # "shared" was set by groupA then overwritten by groupB
  expect_equal(result[["shared"]], "groupB")
})

test_that("parse_excluded_horizons returns empty lists for NULL input", {
  result <- parse_excluded_horizons(NULL, c("rec1", "rec2"))
  expect_equal(result$per_record_excluded, list())
  expect_equal(result$global_excluded, character(0))
})

test_that("parse_excluded_horizons: named list maps per-record correctly", {
  result <- parse_excluded_horizons(
    list(rec1 = "h1", rec2 = c("h1", "h2")),
    all_records = c("rec1", "rec2")
  )
  expect_equal(result$per_record_excluded$rec1, "h1")
  expect_equal(sort(result$per_record_excluded$rec2), c("h1", "h2"))
  expect_equal(result$global_excluded, character(0))
})

test_that("parse_excluded_horizons: unnamed character vector goes to global_excluded", {
  result <- parse_excluded_horizons(c("h1", "h2"), all_records = c("rec1", "rec2"))
  expect_equal(result$global_excluded, c("h1", "h2"))
  expect_equal(result$per_record_excluded, list())
})

test_that("parse_excluded_horizons: named character vector routes named to per-record, unnamed to global", {
  excluded <- c(rec1 = "h1", "h2")
  result <- parse_excluded_horizons(excluded, all_records = c("rec1", "rec2"))
  expect_equal(result$per_record_excluded$rec1, "h1")
  expect_true("h2" %in% result$global_excluded)
})

test_that("parse_excluded_horizons: ambiguous key matching multiple records errors", {
  # Both "rec1" and "rec1b" match key "rec1" via grepl
  expect_error(
    parse_excluded_horizons(list(rec1 = "h1"), all_records = c("rec1", "rec1b")),
    "Ambiguous"
  )
})
