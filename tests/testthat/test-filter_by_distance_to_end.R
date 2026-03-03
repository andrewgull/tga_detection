library(testthat)
library(here)
library(dplyr)
library(readr)

# Source the script to test its functions
source(here("workflow", "scripts", "filter_by_distance_to_end.R"))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

make_test_df <- function() {
  tibble(
    read = c("r1", "r2", "r3", "r4"),
    orient = c("direct", "direct", "reverse", "reverse"),
    end.red = c(1500, 1000, 500, 500),
    start.rep.unit = c(500, 500, 1500, 1000),
    col5 = "dummy",
    col6 = "dummy",
    col7 = "dummy"
  )
}

# ---------------------------------------------------------------------------
# filter_by_distance()
# ---------------------------------------------------------------------------

test_that("filter_by_distance logic for direct orientation", {
  df <- make_test_df() |> filter(orient == "direct")
  # dist = 1200. r1: end.red=1500 > 1200 (KEEP). r2: end.red=1000 < 1200 (FILTER)
  res <- filter_by_distance(df, dist = 1200)
  expect_equal(nrow(res), 1)
  expect_equal(res$read, "r1")
})

test_that("filter_by_distance logic for reverse orientation", {
  df <- make_test_df() |> filter(orient == "reverse")
  # dist = 1200. r3: start.rep.unit=1500 > 1200 (KEEP). r4: start.rep.unit=1000 < 1200 (FILTER)
  res <- filter_by_distance(df, dist = 1200)
  expect_equal(nrow(res), 1)
  expect_equal(res$read, "r3")
})

test_that("filter_by_distance borderline case", {
  df <- tibble(
    read = c("border"),
    orient = c("direct"),
    end.red = c(1200),
    start.rep.unit = c(0)
  )
  # logic is end.red > dist. So 1200 > 1200 is FALSE.
  res <- filter_by_distance(df, dist = 1200)
  expect_equal(nrow(res), 0)
})

# ---------------------------------------------------------------------------
# main()
# ---------------------------------------------------------------------------

test_that("main() processes valid file and parameters", {
  test_df <- make_test_df()
  tmp_file <- tempfile(fileext = ".tsv")
  write_tsv(test_df, tmp_file)
  
  # dist = 1200. Should keep r1 (direct, 1500) and r3 (reverse, 1500)
  res <- main(tmp_file, dist = 1200)
  
  expect_equal(nrow(res), 2)
  expect_setequal(res$read, c("r1", "r3"))
  
  unlink(tmp_file)
})

test_that("main() stops if ncol != 7", {
  bad_df <- tibble(a = 1, b = 2)
  tmp_file <- tempfile(fileext = ".tsv")
  write_tsv(bad_df, tmp_file)
  
  expect_error(main(tmp_file, 1000))
  
  unlink(tmp_file)
})

test_that("main() stops if empty table", {
  empty_df <- make_test_df()[0, ]
  tmp_file <- tempfile(fileext = ".tsv")
  write_tsv(empty_df, tmp_file)
  
  expect_error(main(tmp_file, 1000))
  
  unlink(tmp_file)
})

test_that("main() handles dist as string", {
  test_df <- make_test_df()
  tmp_file <- tempfile(fileext = ".tsv")
  write_tsv(test_df, tmp_file)
  
  res <- main(tmp_file, dist = "1200")
  expect_equal(nrow(res), 2)
  
  unlink(tmp_file)
})
