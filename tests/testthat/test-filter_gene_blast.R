library(testthat)
library(here)
library(dplyr)
library(readr)
library(tidyr)

# Source the script to test its functions
source(here("workflow", "scripts", "filter_gene_blast.R"))

# ---------------------------------------------------------------------------
# is_within()
# ---------------------------------------------------------------------------

test_that("is_within correctly identifies hits inside and outside range", {
  # Direct: 1000 -> 2000
  expect_true(is_within(1200, 1500, 1000, 2000))
  expect_false(is_within(800, 1200, 1000, 2000))
  expect_false(is_within(1200, 2200, 1000, 2000))
  
  # Reverse: 2000 -> 1000
  expect_true(is_within(1500, 1200, 2000, 1000))
  expect_false(is_within(2200, 1800, 2000, 1000))
  expect_false(is_within(1200, 800, 2000, 1000))
})

# ---------------------------------------------------------------------------
# genes_within_flanks()
# ---------------------------------------------------------------------------

test_that("genes_within_flanks detects aberrant hits", {
  # fr_df mirrors the real pipeline output: has start.subject/end.subject from
  # the green FR join, plus end.red from the red FR join.
  # The join inside genes_within_flanks with bla_df (also start.subject/end.subject)
  # produces the .x/.y suffixes that the function selects on.
  fr_df <- tibble(
    subject = c("s1", "s2"),
    start.subject = c(100, 100),   # green start
    end.subject = c(500, 500),     # green end
    end.red = c(2000, 2000),       # red end (boundary of interesting region)
    orientation = c("direct", "direct")
  )
  # is_within checks: is hit inside [end.red, end.subject] = [2000, 500]?
  # Since fr_start > fr_end, the range is [500, 2000].
  bla_df <- tibble(
    subject = c("s1", "s1", "s2", "s2"),
    start.subject = c(1000, 1500, 1000, 2200),
    end.subject = c(1100, 1600, 1100, 2300)
  )

  expect_false(genes_within_flanks("s1", fr_df, bla_df)) # All hits inside [500, 2000]
  expect_true(genes_within_flanks("s2", fr_df, bla_df))  # Hit at 2200-2300 is outside
})

# ---------------------------------------------------------------------------
# main()
# ---------------------------------------------------------------------------

test_that("main() filters subjects with aberrant hits and e-value", {
  fr <- tibble(
    subject = c("s1", "s2"),
    start.subject = 100,
    end.subject = 500,
    end.red = 2000,
    orientation = "direct"
  )
  
  bla <- tibble(
    query = "q",
    subject = c("s1", "s2", "s1"),
    identity = 100,
    length = 100,
    mismatch = 0,
    gaps = 0,
    start.query = 1,
    end.query = 100,
    start.subject = c(1000, 2200, 1200),
    end.subject = c(1100, 2300, 1300),
    e.value = c(0, 0, 1e-10),
    bit.score = 200
  )
  
  # s2 should be filtered because hit at 2200 is outside [500, 2000]
  # Then bla hits for s1 should be filtered by e-value (if we set evalue = 1e-20)
  
  res <- main(bla, fr, evalue = 1e-5)
  expect_equal(nrow(res), 2)
  expect_true(all(res$subject == "s1"))
  
  res_strict <- main(bla, fr, evalue = 1e-20)
  expect_equal(nrow(res_strict), 1)
  expect_equal(res_strict$start.subject, 1000)
})
