library(testthat)
library(here)
library(dplyr)
library(readr)

# Source the script to test its functions
source(here("workflow", "scripts", "filter_fr_hits.R"))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

make_blast_df <- function(subjects = c("s1"), lengths = c(100), evals = c(0), ids = c(100), starts = c(1), ends = c(101)) {
  tibble(
    query = "q",
    subject = subjects,
    identity = ids,
    length = lengths,
    mismatch = 0,
    gaps = 0,
    start.query = 1,
    end.query = lengths,
    start.subject = starts,
    end.subject = ends,
    e.value = evals,
    bit.score = 100
  ) |>
    mutate(orientation = ifelse(start.subject < end.subject, "direct", "reverse"))
}

# ---------------------------------------------------------------------------
# parse_blast()
# ---------------------------------------------------------------------------

test_that("parse_blast correctly parses outfmt 6 and adds columns", {
  df <- make_blast_df(starts = c(1, 100), ends = c(100, 1)) |> select(-orientation)
  tmp_file <- tempfile(fileext = ".tsv")
  write_tsv(df, tmp_file, col_names = FALSE)
  
  res <- parse_blast(tmp_file, "FR_green")
  
  expect_equal(nrow(res), 2)
  expect_equal(ncol(res), 13)
  expect_equal(res$query, c("FR_green", "FR_green"))
  expect_equal(res$orientation, c("direct", "reverse"))
  
  unlink(tmp_file)
})

# ---------------------------------------------------------------------------
# filter_blast_p1()
# ---------------------------------------------------------------------------

test_that("filter_blast_p1 filters correctly", {
  df <- make_blast_df(
    subjects = c("s1", "s2", "s3"),
    lengths = c(100, 50, 100),
    evals = c(0, 0, 1),
    ids = c(100, 100, 70)
  )
  
  # min_len=80, max_e=0.1, min_id=90 -> only s1 should pass
  res <- filter_blast_p1(df, min_len = 80, max_e_value = 0.1, min_identity = 90)
  
  expect_equal(nrow(res), 1)
  expect_equal(res$subject, "s1")
})

# ---------------------------------------------------------------------------
# get_multihits_ids()
# ---------------------------------------------------------------------------

test_that("get_multihits_ids identifies multihits", {
  df <- tibble(subject = c("s1", "s2", "s1", "s3", "s3", "s3"))
  res <- get_multihits_ids(df)
  expect_setequal(res, c("s1", "s3"))
})

# ---------------------------------------------------------------------------
# filter_red_green()
# ---------------------------------------------------------------------------

test_that("filter_red_green joins, filters multihits, and checks orientation", {
  # Red hits (filtered from previous steps, has 'orient' and 'end.red')
  red <- tibble(
    subject = c("s1", "s2", "s3", "s4"),
    orient = c("direct", "direct", "reverse", "direct"),
    end.red = c(1500, 1500, 1500, 1500)
  )
  
  # Green hits (from blast, has 'orientation' and 'start.subject')
  green <- make_blast_df(
    subjects = c("s1", "s2", "s3", "s4", "s4"), # s4 is multihit in green
    starts = c(1000, 2000, 2000, 1000, 1100),
    ends = c(1100, 2100, 1900, 1100, 1200)
  )
  # s1: red.orient=direct, green.orient=direct (MATCH). dist = abs(1500 - 1000) = 500
  # s2: red.orient=direct, green.orient=direct (MATCH). dist = abs(1500 - 2000) = 500
  # s3: red.orient=reverse, green.orient=reverse (MATCH). dist = abs(1500 - 2000) = 500
  # s4: multihit in green -> EXCLUDED
  
  # Also test multihit in red
  red <- bind_rows(red, tibble(subject = "s5", orient = "direct", end.red = 1500))
  red <- bind_rows(red, tibble(subject = "s5", orient = "direct", end.red = 1600))
  green <- bind_rows(green, make_blast_df(subjects = "s5", starts = 1000, ends = 1100))
  # s5: multihit in red -> EXCLUDED
  
  res <- filter_red_green(red, green)
  
  expect_equal(nrow(res), 3)
  expect_setequal(res$subject, c("s1", "s2", "s3"))
  expect_equal(res$distance.btw.FR, c(500, 500, 500))
})

test_that("filter_red_green filters out mismatching orientations", {
  red <- tibble(subject = "s1", orient = "direct", end.red = 1500)
  green <- make_blast_df(subjects = "s1", starts = 2000, ends = 1900) # reverse
  
  res <- filter_red_green(red, green)
  expect_equal(nrow(res), 0)
})

# ---------------------------------------------------------------------------
# main()
# ---------------------------------------------------------------------------

test_that("main() processes valid data", {
  red <- tibble(subject = "s1", orient = "direct", end.red = 1500)
  green <- make_blast_df(subjects = "s1", starts = 1000, ends = 1100)
  
  res <- main(red, green, len = 50, evalue = 1e-5, identity = 90)
  
  expect_equal(nrow(res), 1)
  expect_equal(res$subject, "s1")
})

test_that("main() stops on invalid parameters", {
  red <- tibble(subject = "s1", orient = "direct", end.red = 1500)
  green <- make_blast_df(subjects = "s1", starts = 1000, ends = 1100)
  
  expect_error(suppressWarnings(main(red, green, len = "abc", evalue = 0, identity = 0)))
})

test_that("main() stops if inputs are empty after filtering", {
  red <- tibble(subject = "s1", orient = "direct", end.red = 1500)
  green <- make_blast_df(subjects = "s2", starts = 1000, ends = 1100, lengths = 10) # will be filtered
  
  expect_error(main(red, green, len = 50, evalue = 1e-5, identity = 90))
})
