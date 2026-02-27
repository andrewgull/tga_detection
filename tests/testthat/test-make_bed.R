library(testthat)
library(readr)
library(dplyr)

# Source the script to test its functions
source("../../workflow/scripts/make_bed.R")

test_that("is_filtered detects 'filtered' in filename", {
  expect_true(is_filtered("blast_results_filtered.tsv"))
  expect_false(is_filtered("blast_results.tsv"))
  expect_true(is_filtered("path/to/filtered_hits.txt"))
})

test_that("read_blast reads filtered and non-filtered tables correctly", {
  # Create a mock non-filtered blast table (12 columns, no headers)
  non_filtered_file <- tempfile(fileext = ".tsv")
  mock_data_non_filtered <- tibble(
    X1 = "query1", X2 = "subject1", X3 = 99.0, X4 = 100,
    X5 = 0, X6 = 0, X7 = 1, X8 = 100, X9 = 500, X10 = 600,
    X11 = 1e-10, X12 = 200.0
  )
  write_tsv(mock_data_non_filtered, non_filtered_file, col_names = FALSE)
  
  # Create a mock filtered blast table (with headers)
  filtered_file <- tempfile(pattern = "filtered", fileext = ".tsv")
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
  
  unlink(non_filtered_file)
  unlink(filtered_file)
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
  
  # Positive strand check (s1: 100-200)
  # BED is 0-indexed, [start, end)
  # Script logic: start = start.subject - 1, end = end.subject - 1
  # For q1: start.subject=100, end.subject=200 -> start=99, end=199
  # Wait, standard BED transformation usually means:
  # if start < end: start = start-1, end = end
  # Let's check script's logic:
  # mutate(strand = if_else(start.subject < end.subject, "+", "-"),
  #        start = start.subject - 1,
  #        end = end.subject - 1)
  # If q1 (pos): start=99, end=199.
  # If q2 (neg): start.subject=200, end.subject=100 -> strand="-", start=199, end=99.
  # Then it swaps for neg: rename("start" = end, "end" = start) -> start=99, end=199.
  
  q1_bed <- bed %>% filter(name == "q1")
  expect_equal(q1_bed$strand, "+")
  expect_equal(q1_bed$start, 99)
  expect_equal(q1_bed$end, 199)
  
  q2_bed <- bed %>% filter(name == "q2")
  expect_equal(q2_bed$strand, "-")
  expect_equal(q2_bed$start, 99)
  expect_equal(q2_bed$end, 199)
})
