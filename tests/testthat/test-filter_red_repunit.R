library(testthat)
library(here)
library(dplyr)
library(readr)

# Source the script to test its functions
source(here("workflow", "scripts", "filter_red_repunit.R"))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

make_blast_df <- function(subjects = c("s1"), lengths = c(100), starts = c(1), ends = c(101)) {
  tibble(
    query = "q",
    subject = subjects,
    identity = 100,
    length = lengths,
    mismatch = 0,
    gaps = 0,
    start.query = 1,
    end.query = lengths,
    start.subject = starts,
    end.subject = ends,
    e.value = 0,
    bit.score = 100
  ) |>
    mutate(orientation = ifelse(start.subject < end.subject, "direct", "reverse"))
}

# ---------------------------------------------------------------------------
# filter_ru_fr()
# ---------------------------------------------------------------------------

test_that("filter_ru_fr joins and calculates distance correctly (direct)", {
  fr <- make_blast_df(subjects = "s1", starts = 1000, ends = 1100) |> rename(orientation.x = orientation)
  ru <- make_blast_df(subjects = "s1", starts = 800, ends = 900) |> rename(orientation.y = orientation)
  
  # direct: start.subject.x - end.subject.y + 1 = 1000 - 900 + 1 = 101
  res <- filter_ru_fr(fr, ru, max_distance = 200)
  
  expect_equal(nrow(res), 1)
  expect_equal(res$distance, 101)
})

test_that("filter_ru_fr joins and calculates distance correctly (reverse)", {
  fr <- make_blast_df(subjects = "s1", starts = 1000, ends = 900) |> rename(orientation.x = orientation)
  ru <- make_blast_df(subjects = "s1", starts = 1200, ends = 1100) |> rename(orientation.y = orientation)
  
  # reverse: end.subject.y - start.subject.x + 1 = 1100 - 1000 + 1 = 101
  res <- filter_ru_fr(fr, ru, max_distance = 200)
  
  expect_equal(nrow(res), 1)
  expect_equal(res$distance, 101)
})

test_that("filter_ru_fr filters by max_distance", {
  fr <- make_blast_df(subjects = "s1", starts = 1000, ends = 1100) |> rename(orientation.x = orientation)
  ru <- make_blast_df(subjects = "s1", starts = 500, ends = 600) |> rename(orientation.y = orientation)
  
  # distance = 1000 - 600 + 1 = 401
  res <- filter_ru_fr(fr, ru, max_distance = 200)
  expect_equal(nrow(res), 0)
})

test_that("filter_ru_fr filters by matching orientation", {
  fr <- make_blast_df(subjects = "s1", starts = 1000, ends = 1100) |> rename(orientation.x = orientation) # direct
  ru <- make_blast_df(subjects = "s1", starts = 900, ends = 800) |> rename(orientation.y = orientation)   # reverse
  
  res <- filter_ru_fr(fr, ru, max_distance = 1000)
  expect_equal(nrow(res), 0)
})

# ---------------------------------------------------------------------------
# main()
# ---------------------------------------------------------------------------

test_that("main() processes valid files and parameters", {
  # Mocking files is hard, so we'll test the core logic via mock paths if possible, 
  # but here we'll just test that main handles the flow if we provide real temp files.
  
  fr_df <- make_blast_df(subjects = "s1", starts = 1000, ends = 1100) |> select(-orientation)
  ru_df <- make_blast_df(subjects = "s1", starts = 900, ends = 950) |> select(-orientation)
  
  fr_file <- tempfile(fileext = ".tsv")
  ru_file <- tempfile(fileext = ".tsv")
  write_tsv(fr_df, fr_file, col_names = FALSE)
  write_tsv(ru_df, ru_file, col_names = FALSE)
  
  # dist = 1000 - 950 + 1 = 51
  res <- main(fr_file, ru_file, ru_len = 40, fr_len = 80, e = 1e-5, identity = 90, maxd = 100)
  
  expect_equal(nrow(res), 1)
  expect_equal(res$subject, "s1")
  expect_equal(res$dist, 51)
  expect_named(res, c("subject", "start.red", "end.red", "start.rep.unit", "end.rep.unit", "dist", "orient"))
  
  unlink(fr_file)
  unlink(ru_file)
})

test_that("main() stops on invalid parameters", {
  expect_error(suppressWarnings(main("dummy", "dummy", ru_len = "abc", fr_len = 10, e = 0, identity = 0, maxd = 0)))
})
