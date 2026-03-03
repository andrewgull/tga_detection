library(testthat)
library(here)
library(dplyr)
library(readr)

# Source the script to test its functions
source(here("workflow", "scripts", "filter_by_distance_between_fr.R"))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

make_fr_ru <- function(subjects, distances) {
  tibble(subject = subjects, distance.btw.FR = distances)
}

make_bla_counts <- function(subjects, counts) {
  tibble(subject = subjects, n.blaSHV.merged = counts)
}

# ---------------------------------------------------------------------------
# filter_by_distance() — core logic
# ---------------------------------------------------------------------------

test_that("filter_by_distance keeps exact matches", {
  # base_len = 4000, ru_len = 1000
  # dist = 4000 -> exp = (4000-4000)/1000 + 1 = 1
  # observed = 1 -> diff = 0
  fr_ru <- make_fr_ru("read1", 4000)
  bla_counts <- make_bla_counts("read1", 1)
  
  result <- filter_by_distance(fr_ru, bla_counts, base_len = 4000, ru_len = 1000)
  
  expect_equal(nrow(result), 1)
  expect_equal(result$n.blaSHV.merged, 1)
})

test_that("filter_by_distance keeps difference of 1", {
  # dist = 5000 -> exp = (5000-4000)/1000 + 1 = 2
  # observed = 1 -> diff = 1 -> Keep
  fr_ru <- make_fr_ru("read1", 5000)
  bla_counts <- make_bla_counts("read1", 1)
  
  result <- filter_by_distance(fr_ru, bla_counts, base_len = 4000, ru_len = 1000)
  
  expect_equal(nrow(result), 1)
})

test_that("filter_by_distance filters out difference > 1", {
  # dist = 6000 -> exp = (6000-4000)/1000 + 1 = 3
  # observed = 1 -> diff = 2 -> Filter
  fr_ru <- make_fr_ru("read1", 6000)
  bla_counts <- make_bla_counts("read1", 1)
  
  result <- filter_by_distance(fr_ru, bla_counts, base_len = 4000, ru_len = 1000)
  
  expect_equal(nrow(result), 0)
})

test_that("filter_by_distance handles rounding correctly", {
  # base_len = 0, ru_len = 10
  # dist = 5 -> exp = (5-0)/10 + 1 = 1.5
  # observed = 1 -> diff = 0.5 -> round(0.5) = 0 -> Keep
  fr_ru <- make_fr_ru("read1", 5)
  bla_counts <- make_bla_counts("read1", 1)
  
  result <- filter_by_distance(fr_ru, bla_counts, base_len = 0, ru_len = 10)
  expect_equal(nrow(result), 1)
  
  # dist = 15 -> exp = (15-0)/10 + 1 = 2.5
  # observed = 1 -> diff = 1.5 -> round(1.5) = 2 -> Filter
  fr_ru <- make_fr_ru("read2", 15)
  bla_counts <- make_bla_counts("read2", 1)
  
  result <- filter_by_distance(fr_ru, bla_counts, base_len = 0, ru_len = 10)
  expect_equal(nrow(result), 0)
})

test_that("filter_by_distance filters out NAs from missing join", {
  fr_ru <- make_fr_ru(c("read1", "read2"), c(4000, 4000))
  bla_counts <- make_bla_counts("read1", 1) # read2 is missing
  
  result <- filter_by_distance(fr_ru, bla_counts, base_len = 4000, ru_len = 1000)
  
  expect_equal(nrow(result), 1)
  expect_equal(result$subject, "read1")
})

# ---------------------------------------------------------------------------
# filter_by_distance() — parameter validation
# ---------------------------------------------------------------------------

test_that("filter_by_distance errors on invalid parameters", {
  fr_ru <- make_fr_ru("read1", 4000)
  bla_counts <- make_bla_counts("read1", 1)
  
  expect_error(filter_by_distance(fr_ru, bla_counts, base_len = NA))
  expect_error(filter_by_distance(fr_ru, bla_counts, ru_len = NA))
  expect_error(suppressWarnings(filter_by_distance(fr_ru, bla_counts, base_len = "abc")))
})

test_that("filter_by_distance accepts numeric strings for parameters", {
  fr_ru <- make_fr_ru("read1", 4000)
  bla_counts <- make_bla_counts("read1", 1)

  result <- filter_by_distance(fr_ru, bla_counts, base_len = "4000", ru_len = "1000")
  expect_equal(nrow(result), 1)
})

# ---------------------------------------------------------------------------
# filter_by_distance() — negative differences (observed > expected)
# ---------------------------------------------------------------------------

test_that("filter_by_distance keeps reads where observed exceeds expected by 1", {
  # dist = 4000 -> exp = (4000-4000)/1000 + 1 = 1; observed = 2 -> diff = -1 -> keep
  fr_ru <- make_fr_ru("read1", 4000)
  bla_counts <- make_bla_counts("read1", 2)

  result <- filter_by_distance(fr_ru, bla_counts, base_len = 4000, ru_len = 1000)

  expect_equal(nrow(result), 1)
  expect_equal(result$n.blaSHV.merged, 2)
})

test_that("filter_by_distance removes reads where observed exceeds expected by 2", {
  # dist = 4000 -> exp = 1; observed = 3 -> diff = -2 -> filter
  fr_ru <- make_fr_ru("read1", 4000)
  bla_counts <- make_bla_counts("read1", 3)

  result <- filter_by_distance(fr_ru, bla_counts, base_len = 4000, ru_len = 1000)

  expect_equal(nrow(result), 0)
})

# ---------------------------------------------------------------------------
# filter_by_distance() — NA in n.blaSHV.merged
# ---------------------------------------------------------------------------

test_that("filter_by_distance removes rows where n.blaSHV.merged is NA", {
  # NA n.blaSHV.merged (e.g. from count_merged_gene_hits) -> difference = NA
  # abs(NA) <= 1 is NA -> dplyr::filter() drops NA rows
  fr_ru <- make_fr_ru("read1", 4000)
  bla_counts <- tibble(subject = "read1", n.blaSHV.merged = NA_integer_)

  result <- filter_by_distance(fr_ru, bla_counts, base_len = 4000, ru_len = 1000)

  expect_equal(nrow(result), 0)
})

# ---------------------------------------------------------------------------
# filter_by_distance() — output schema
# ---------------------------------------------------------------------------

test_that("filter_by_distance returns exactly 2 named columns", {
  fr_ru <- make_fr_ru("read1", 4000)
  bla_counts <- make_bla_counts("read1", 1)

  result <- filter_by_distance(fr_ru, bla_counts, base_len = 4000, ru_len = 1000)

  expect_equal(ncol(result), 2)
  expect_named(result, c("subject", "n.blaSHV.merged"))
})

# ---------------------------------------------------------------------------
# filter_by_distance() — edge: empty input
# ---------------------------------------------------------------------------

test_that("filter_by_distance returns zero rows when fr_ru_filt is empty", {
  fr_ru <- tibble(subject = character(0), distance.btw.FR = numeric(0))
  bla_counts <- make_bla_counts("read1", 1)

  result <- filter_by_distance(fr_ru, bla_counts, base_len = 4000, ru_len = 1000)

  expect_equal(nrow(result), 0)
})

# ---------------------------------------------------------------------------
# filter_by_distance() — edge: ru_len = 0
# ---------------------------------------------------------------------------

test_that("filter_by_distance errors when ru_len is zero", {
  # ru_len = 0 is rejected explicitly to prevent silent division by zero
  # (which would otherwise produce NaN and silently drop all rows)
  fr_ru <- make_fr_ru("read1", 4000)
  bla_counts <- make_bla_counts("read1", 1)

  expect_error(
    filter_by_distance(fr_ru, bla_counts, base_len = 4000, ru_len = 0)
  )
})
