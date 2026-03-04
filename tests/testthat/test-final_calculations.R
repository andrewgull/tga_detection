library(testthat)
library(here)
library(dplyr)
library(tidyr)

# Source the script to test its functions
# Note: final_calculations.R uses %>% which might require dplyr or magrittr explicitly
# but it also loads dplyr.
source(here("workflow", "scripts", "final_calculations.R"))

# ---------------------------------------------------------------------------
# counts_freq_obs()
# ---------------------------------------------------------------------------

test_that("counts_freq_obs calculates frequencies correctly", {
  bla_cn <- tibble(n.blaSHV.merged = c(1, 1, 2, 3, NA))
  
  res <- counts_freq_obs(bla_cn)
  
  expect_equal(nrow(res), 3)
  expect_equal(res$CN, c(1, 2, 3))
  expect_equal(res$counts_obs, c(2, 1, 1))
  expect_equal(res$freq_obs, c(0.5, 0.25, 0.25))
})

# ---------------------------------------------------------------------------
# freq_theor()
# ---------------------------------------------------------------------------

test_that("freq_theor calculates theoretical frequencies based on CN=0", {
  cn_bins <- tibble(
    CN = c(0, 1, 2),
    n_reads_theoretical = c(100, 50, 10)
  )
  
  res <- freq_theor(cn_bins)
  
  expect_equal(nrow(res), 3)
  expect_equal(res$freq_theoretical, c(1.0, 0.5, 0.1))
})

test_that("freq_theor errors if CN=0 is missing or duplicated", {
  expect_error(freq_theor(tibble(CN = 1, n_reads_theoretical = 100)))
  expect_error(freq_theor(tibble(CN = c(0, 0), n_reads_theoretical = c(100, 100))))
})

# ---------------------------------------------------------------------------
# main()
# ---------------------------------------------------------------------------

test_that("main() combines observed and theoretical data correctly", {
  cn_bins <- tibble(
    CN = c(0, 1, 2),
    n_reads_theoretical = c(100, 50, 10)
  )
  # Theoretical freqs: 1, 0.5, 0.1
  
  bla_cn <- tibble(n.blaSHV.merged = c(1, 1, 2))
  # Observed: CN=1: 2 (freq_obs=2/3), CN=2: 1 (freq_obs=1/3)
  
  res <- main(cn_bins, bla_cn)
  
  # Check CN=1
  # counts_obs = 2, freq_theoretical = 0.5 -> counts_corrected = 2 / 0.5 = 4
  # Check CN=2
  # counts_obs = 1, freq_theoretical = 0.1 -> counts_corrected = 1 / 0.1 = 10
  # Total counts_corrected = 14
  # freq_corrected for CN=1 = 4/14 = 0.2857...
  # freq_corrected for CN=2 = 10/14 = 0.7142...
  
  expect_equal(res$counts_corrected[res$CN == 1], 4)
  expect_equal(res$counts_corrected[res$CN == 2], 10)
  expect_equal(res$freq_corrected[res$CN == 1], 4/14)
  expect_equal(res$freq_corrected[res$CN == 2], 10/14)
  
  # Detection limits
  expect_equal(res$detection_limit[res$CN == 1], 1/50)
  expect_equal(res$detection_limit[res$CN == 2], 1/10)
})

test_that("main() handles CN variants not observed", {
  cn_bins <- tibble(
    CN = c(0, 1, 2),
    n_reads_theoretical = c(100, 50, 10)
  )
  bla_cn <- tibble(n.blaSHV.merged = c(1, 1)) # CN=2 not observed
  
  res <- main(cn_bins, bla_cn)
  expect_equal(res$counts_obs[res$CN == 2], 0)
  expect_equal(res$counts_corrected[res$CN == 2], 0)
  expect_equal(res$freq_corrected[res$CN == 2], 0)
})

test_that("main() handles no observed counts", {
  cn_bins <- tibble(
    CN = c(0, 1, 2),
    n_reads_theoretical = c(100, 50, 10)
  )
  bla_cn <- tibble(n.blaSHV.merged = numeric(0))
  
  res <- main(cn_bins, bla_cn)
  expect_true(all(is.na(res$freq_corrected)))
})
