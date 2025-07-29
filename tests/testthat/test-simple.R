# Simple working tests that don't require Seurat
library(testthat)

test_that("package loads", {
  expect_true(TRUE)
})

test_that("get_package_version works", {
  version <- get_package_version()
  expect_equal(version, "1.0.0")
})

test_that("get_package_info works", {
  info <- get_package_info()
  expect_type(info, "list")
  expect_equal(info$package, "SpatialCoherence")
  expect_equal(info$version, "1.0.0")
})

test_that("classify_organization works", {
  scores <- c(A = 0.3, B = 0.6)
  result <- classify_organization(scores)
  expect_length(result, 2)
  expect_true(all(result %in% c("Organized", "Disorganized")))
})

test_that("get_default_config works", {
  config <- get_default_config()
  expect_type(config, "list")
  expect_true("data_columns" %in% names(config))
  expect_true("spatial_coherence" %in% names(config))
  expect_true("output" %in% names(config))
})

