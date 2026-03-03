library(testthat)
library(here)
library(readr)
library(dplyr)

# Source the script to test its functions
source(here("workflow", "scripts", "make_bed.R"))

test_that("is_filtered detects 'filtered' in filename", {
  expect_true(is_filtered("blast_results_filtered.tsv"))
  expect_false(is_filtered("blast_results.tsv"))
  expect_true(is_filtered("path/to/filtered_hits.txt"))
})

test_that("read_blast reads filtered and non-filtered tables correctly", {
  # Create a mock non-filtered blast table (12 columns, no headers)
  non_filtered_file <- tempfile(fileext = ".tsv")
  withr::defer(unlink(non_filtered_file))
  mock_data_non_filtered <- tibble(
    X1 = "query1", X2 = "subject1", X3 = 99.0, X4 = 100,
    X5 = 0, X6 = 0, X7 = 1, X8 = 100, X9 = 500, X10 = 600,
    X11 = 1e-10, X12 = 200.0
  )
  write_tsv(mock_data_non_filtered, non_filtered_file, col_names = FALSE)

  # Create a mock filtered blast table (with headers)
  filtered_file <- tempfile(pattern = "filtered", fileext = ".tsv")
  withr::defer(unlink(filtered_file))
  mock_data_filtered <- tibble(
    query = "query1", subject = "subject1", identity = 99.0, length = 100,
    mismatch = 0, gaps = 0, start.query = 1, end.query = 100,
    start.subject = 500, end.subject = 600, e.value = 1e-10, bit.score = 200.0
  )
  write_tsv(mock_data_filtered, filtered_file, col_names = TRUE)

  # Test reading non-filtered
  res_non_filtered <- read_blast(non_filtered_file)
  expect_equal(colnames(res_non_filtered)[1], "query")
  expect_equal(nrow(res_non_filtered), 1)
  expect_equal(res_non_filtered$subject[1], "subject1")

  # Test reading filtered
  res_filtered <- read_blast(filtered_file)
  expect_equal(colnames(res_filtered)[1], "query")
  expect_equal(nrow(res_filtered), 1)
  expect_equal(res_filtered$subject[1], "subject1")
})

test_that("blast2bed converts coordinates and handles strands correctly", {
  mock_blast <- tibble(
    query = c("q1", "q2"),
    subject = c("s1", "s2"),
    identity = c(100, 100),
    length = c(100, 100),
    mismatch = c(0, 0),
    gaps = c(0, 0),
    start.query = c(1, 1),
    end.query = c(100, 100),
    start.subject = c(100, 200), # q1: positive strand, q2: negative strand
    end.subject = c(200, 100),
    e.value = c(0, 0),
    bit.score = c(100, 100)
  )

  bed <- blast2bed(mock_blast)

  expect_equal(nrow(bed), 2)
  expect_equal(colnames(bed), c("chrom", "start", "end", "name", "score", "strand"))

  q1_bed <- bed %>% filter(name == "q1")
  expect_equal(q1_bed$strand, "+")
  # BED is 0-indexed half-open: start = pmin - 1, end = pmax (not pmax - 1)
  # q1: start.subject=100, end.subject=200 -> start=99, end=200
  expect_equal(q1_bed$start, 99)
  expect_equal(q1_bed$end, 200)

  q2_bed <- bed %>% filter(name == "q2")
  expect_equal(q2_bed$strand, "-")
  # q2: start.subject=200, end.subject=100 -> start=pmin(200,100)-1=99, end=pmax(200,100)=200
  expect_equal(q2_bed$start, 99)
  expect_equal(q2_bed$end, 200)
})

test_that("is_filtered is case-sensitive", {
  expect_false(is_filtered("blast_FILTERED.tsv"))
  expect_false(is_filtered("blast_Filtered.tsv"))
})

test_that("is_filtered handles edge cases", {
  expect_false(is_filtered(""))
  expect_false(is_filtered("filter.tsv"))
  expect_true(is_filtered("filteredfiltered.tsv"))
})

test_that("read_blast errors on empty non-filtered table", {
  empty_file <- tempfile(fileext = ".tsv")
  withr::defer(unlink(empty_file))
  write_tsv(tibble(), empty_file, col_names = FALSE)
  # tibble warns about the missing X1 column before stopifnot fires
  expect_error(suppressWarnings(read_blast(empty_file)))
})

test_that("read_blast errors when non-filtered table has multiple unique X1 values", {
  bad_file <- tempfile(fileext = ".tsv")
  withr::defer(unlink(bad_file))
  mock_bad <- tibble(
    X1 = c("query1", "query2"), X2 = "s1", X3 = 99.0, X4 = 100,
    X5 = 0, X6 = 0, X7 = 1, X8 = 100, X9 = 500, X10 = 600,
    X11 = 1e-10, X12 = 200.0
  )
  write_tsv(mock_bad, bad_file, col_names = FALSE)
  expect_error(read_blast(bad_file))
})

test_that("read_blast errors when filtered table first column is not 'query'", {
  bad_filtered <- tempfile(pattern = "filtered", fileext = ".tsv")
  withr::defer(unlink(bad_filtered))
  mock_bad <- tibble(wrong_col = "oops", subject = "s1")
  write_tsv(mock_bad, bad_filtered, col_names = TRUE)
  expect_error(read_blast(bad_filtered))
})

test_that("read_blast errors when non-filtered table has wrong number of columns", {
  bad_file <- tempfile(fileext = ".tsv")
  withr::defer(unlink(bad_file))
  mock_bad <- tibble(X1 = "query1", X2 = "s1", X3 = 99.0)
  write_tsv(mock_bad, bad_file, col_names = FALSE)
  expect_error(read_blast(bad_file))
})

test_that("blast2bed assigns '-' strand when start.subject equals end.subject", {
  mock <- tibble(
    query = "q1", subject = "s1", identity = 100, length = 1,
    mismatch = 0, gaps = 0, start.query = 1, end.query = 1,
    start.subject = 100, end.subject = 100,
    e.value = 0, bit.score = 50
  )
  bed <- blast2bed(mock)
  # start.subject < end.subject is FALSE when equal, so strand = "-"
  expect_equal(bed$strand, "-")
  expect_equal(bed$start, 99)  # pmin(100, 100) - 1
  expect_equal(bed$end, 100)   # pmax(100, 100)
})

test_that("blast2bed always returns exactly 6 named columns", {
  mock <- tibble(
    query = "q1", subject = "s1", identity = 100, length = 100,
    mismatch = 0, gaps = 0, start.query = 1, end.query = 100,
    start.subject = 1, end.subject = 100, e.value = 0, bit.score = 100
  )
  bed <- blast2bed(mock)
  expect_equal(ncol(bed), 6)
  expect_named(bed, c("chrom", "start", "end", "name", "score", "strand"))
})
