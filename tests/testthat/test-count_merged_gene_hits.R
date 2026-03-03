library(testthat)
library(here)
library(dplyr)

# Source the script to test its functions
source(here("workflow", "scripts", "count_merged_gene_hits.R"))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

make_merge <- function(X1, X2, X3) {
  tibble(X1 = X1, X2 = X2, X3 = X3)
}

make_filt <- function(subjects) {
  tibble(subject = subjects)
}

# ---------------------------------------------------------------------------
# main() — output schema
# ---------------------------------------------------------------------------

test_that("main returns a tibble with exactly 2 named columns", {
  bed <- make_merge("read1", 0, 99)
  filt <- make_filt("read1")

  result <- main(bed, filt, 100)

  expect_equal(ncol(result), 2)
  expect_named(result, c("subject", "n.blaSHV.merged"))
})

# ---------------------------------------------------------------------------
# main() — basic counting
# ---------------------------------------------------------------------------

test_that("main counts a single interval correctly", {
  # interval length: X3 - X2 + 1 = 99 - 0 + 1 = 100; 100 / 100 = 1
  bed <- make_merge("read1", 0, 99)
  filt <- make_filt("read1")

  result <- main(bed, filt, 100)

  expect_equal(nrow(result), 1)
  expect_equal(result$n.blaSHV.merged[1], 1)
})

test_that("main sums multiple intervals for the same subject", {
  # two intervals of 100 bp each -> sum = 200; 200 / 100 = 2
  bed <- make_merge(c("read1", "read1"), c(0, 200), c(99, 299))
  filt <- make_filt("read1")

  result <- main(bed, filt, 100)

  expect_equal(nrow(result), 1)
  expect_equal(result$n.blaSHV.merged[1], 2)
})

test_that("main handles multiple subjects independently", {
  # read1: 100 bp / 100 = 1; read2: 200 bp / 100 = 2
  bed <- make_merge(c("read1", "read2"), c(0, 0), c(99, 199))
  filt <- make_filt(c("read1", "read2"))

  result <- main(bed, filt, 100)

  expect_equal(nrow(result), 2)
  expect_equal(result$n.blaSHV.merged[result$subject == "read1"], 1)
  expect_equal(result$n.blaSHV.merged[result$subject == "read2"], 2)
})

test_that("main counts single-base intervals (X2 == X3) as length 1", {
  # X3 - X2 + 1 = 100 - 100 + 1 = 1; 1 / 1 = 1
  bed <- make_merge("read1", 100, 100)
  filt <- make_filt("read1")

  result <- main(bed, filt, 1)

  expect_equal(result$n.blaSHV.merged[1], 1)
})

# ---------------------------------------------------------------------------
# main() — right-join filtering behaviour
# ---------------------------------------------------------------------------

test_that("main drops subjects in df_merge that are absent from df_filt", {
  bed <- make_merge(c("read1", "read2"), c(0, 0), c(99, 99))
  filt <- make_filt("read1")  # read2 is excluded

  result <- main(bed, filt, 100)

  expect_equal(nrow(result), 1)
  expect_equal(result$subject, "read1")
})

test_that("main returns NA for subjects in df_filt absent from df_merge", {
  # read2 is in filt but has no merged intervals -> NA after right_join
  bed <- make_merge("read1", 0, 99)
  filt <- make_filt(c("read1", "read2"))

  result <- main(bed, filt, 100)

  expect_equal(nrow(result), 2)
  expect_true(is.na(result$n.blaSHV.merged[result$subject == "read2"]))
})

test_that("main returns zero rows when df_filt is empty", {
  bed <- make_merge("read1", 0, 99)
  filt <- tibble(subject = character(0))

  result <- main(bed, filt, 100)

  expect_equal(nrow(result), 0)
})

# ---------------------------------------------------------------------------
# main() — rounding (R uses banker's rounding for round())
# ---------------------------------------------------------------------------

test_that("main rounds 0.5 down to 0 (banker's rounding)", {
  # interval length 50, len 100 -> 50 / 100 = 0.5 -> rounds to 0
  bed <- make_merge("read1", 0, 49)
  filt <- make_filt("read1")

  result <- main(bed, filt, 100)

  expect_equal(result$n.blaSHV.merged[1], 0)
})

test_that("main rounds 1.5 up to 2 (banker's rounding to even)", {
  # interval length 150, len 100 -> 150 / 100 = 1.5 -> rounds to 2
  bed <- make_merge("read1", 0, 149)
  filt <- make_filt("read1")

  result <- main(bed, filt, 100)

  expect_equal(result$n.blaSHV.merged[1], 2)
})

# ---------------------------------------------------------------------------
# main() — len parameter validation
# ---------------------------------------------------------------------------

test_that("main errors when len is NA", {
  bed <- make_merge("read1", 0, 99)
  filt <- make_filt("read1")
  expect_error(main(bed, filt, NA))
})

test_that("main errors when len is zero", {
  bed <- make_merge("read1", 0, 99)
  filt <- make_filt("read1")
  expect_error(main(bed, filt, 0))
})

test_that("main errors when len is negative", {
  bed <- make_merge("read1", 0, 99)
  filt <- make_filt("read1")
  expect_error(main(bed, filt, -1))
})

test_that("main errors when len is a non-numeric string", {
  bed <- make_merge("read1", 0, 99)
  filt <- make_filt("read1")
  # as.integer("abc") produces NA with a warning; stopifnot then fires
  expect_error(suppressWarnings(main(bed, filt, "abc")))
})

test_that("main accepts a numeric string for len", {
  bed <- make_merge("read1", 0, 99)
  filt <- make_filt("read1")

  result <- main(bed, filt, "100")

  expect_equal(result$n.blaSHV.merged[1], 1)
})

test_that("main truncates a float len via as.integer", {
  # as.integer(1.9) = 1; 100 / 1 = 100
  bed <- make_merge("read1", 0, 99)
  filt <- make_filt("read1")

  result <- main(bed, filt, 1.9)

  expect_equal(result$n.blaSHV.merged[1], 100)
})
