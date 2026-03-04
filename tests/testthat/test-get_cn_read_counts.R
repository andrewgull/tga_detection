library(testthat)
library(here)
library(dplyr)
library(purrr)

# Source the script to test its functions
# Note: Sourcing might fail if Biostrings is missing, 
# but we assume the environment is set up.
source(here("workflow", "scripts", "get_cn_read_counts.R"))

# ---------------------------------------------------------------------------
# count_reads_cn()
# ---------------------------------------------------------------------------

test_that("count_reads_cn logic for direct orientation", {
  df <- tibble(end.red = c(1500, 2000, 5000))
  # b=1000, i=1000. 
  # CN 0: 1000 -> 3 reads (all > 1000)
  # CN 1: 2000 -> 1 read (5000 > 2000)
  expect_equal(count_reads_cn(df, cn = 0, b = 1000, i = 1000, direct = TRUE), 3)
  expect_equal(count_reads_cn(df, cn = 1, b = 1000, i = 1000, direct = TRUE), 1)
})

test_that("count_reads_cn logic for reverse orientation", {
  df <- tibble(end.red = c(1000, 1000), read.len = c(5000, 2500))
  # b=1000, i=1000.
  # r1: 5000 - 1000 = 4000
  # r2: 2500 - 1000 = 1500
  # CN 0: 1000 -> 2 reads (4000 > 1000, 1500 > 1000)
  # CN 1: 2000 -> 1 read (4000 > 2000)
  expect_equal(count_reads_cn(df, cn = 0, b = 1000, i = 1000, direct = FALSE), 2)
  expect_equal(count_reads_cn(df, cn = 1, b = 1000, i = 1000, direct = FALSE), 1)
})

# ---------------------------------------------------------------------------
# separate()
# ---------------------------------------------------------------------------

test_that("separate() correctly splits and joins for reverse", {
  df <- tibble(
    subject = c("s1", "s2"),
    orient = c("direct", "reverse"),
    end.red = 1000
  )
  reads_len <- tibble(subject = "s2", read.len = 5000)
  
  res_dir <- separate(df, "direct")
  expect_equal(nrow(res_dir), 1)
  expect_equal(res_dir$subject, "s1")
  
  res_rev <- separate(df, "reverse", reads_len = reads_len)
  expect_equal(nrow(res_rev), 1)
  expect_equal(res_rev$subject, "s2")
  expect_equal(res_rev$read.len, 5000)
})

# ---------------------------------------------------------------------------
# count_reads()
# ---------------------------------------------------------------------------

test_that("count_reads calls count_reads_cn for each variant", {
  df <- tibble(end.red = c(1500, 2500))
  cn_vars <- c(0, 1)
  # b=1000, i=1000. CN 0 (>1000): 2. CN 1 (>2000): 1.
  res <- count_reads(cn_vars, df, min_len = 1000, increment = 1000, direct = TRUE)
  expect_equal(res, c(2, 1))
})

# ---------------------------------------------------------------------------
# main() logic test (without Biostrings call if possible)
# ---------------------------------------------------------------------------

test_that("main() flow works", {
  # We mock read_length to avoid Biostrings dependency
  mock_read_length <- function(fasta) {
    tibble(subject = "s2", read.len = 5000)
  }
  
  # Inject mock into environment
  # We need to be careful if we source the script, we might want to override it
  # in the environment where main runs.
  
  # Create dummy input table
  in_df <- tibble(
    subject = c("s1", "s2"),
    start.red = 100, end.red = 1500,
    start.rep.unit = 200, end.rep.unit = 300,
    dist = 50, orient = c("direct", "reverse")
  )
  in_file <- tempfile(fileext = ".tsv")
  write_tsv(in_df, in_file)
  
  # Since main() calls read_length, we temporarily override it
  orig_read_length <- read_length
  assign("read_length", mock_read_length, envir = .GlobalEnv)
  
  # max_cn=1, min_len=1000, increment=1000
  # Direct hit (s1): end.red=1500. CN0 (>1000): 1. CN1 (>2000): 0.
  # Reverse hit (s2): read.len-end.red = 5000-1500 = 3500. CN0 (>1000): 1. CN1 (>2000): 1.
  # Totals: CN0 = 2, CN1 = 1.
  
  res <- main(in_file, "dummy.fasta", max_cn = 1, min_len = 1000, increment = 1000)
  
  expect_equal(res$CN, c(0, 1))
  expect_equal(res$n_reads_theoretical, c(2, 1))
  
  # Restore
  assign("read_length", orig_read_length, envir = .GlobalEnv)
  unlink(in_file)
})
