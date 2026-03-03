#' @title Final calculations of gene frequencies
#' @description This script performs final calculations of gene counts,
#' observed/corrected frequencies, and detection limits.
#' @section Input:
#' 1. cn_bins: table with number of reads capable of carrying each CN variant.
#' 2. bla_cn: table with gene counts for reads that passed filters.
#' @section Output:
#' Table with frequencies and detection limits.
#' @md

# --- 1. Load Libraries ---
suppressPackageStartupMessages(library(dplyr))
library(readr)
library(tidyr)

# --- 2. Define Functions ---
#' Calculate observed frequency of gene counts
#'
#' @param bla_cn Data frame/tibble. Observed gene counts.
#'
#' @return A tibble with observed frequencies.
#'
counts_freq_obs <- function(bla_cn) {
  bla_cn |>
    filter(!is.na(n.blaSHV.merged)) |>
    group_by(n.blaSHV.merged) |>
    count(name = "counts") |>
    ungroup() |>
    rename(
      "CN" = n.blaSHV.merged,
      "counts_obs" = counts
    ) |>
    mutate(freq_obs = counts_obs / sum(counts_obs))
}

#' Calculate theoretical frequency from bins
#'
#' @param cn_bins Data frame/tibble. Theoretical bins.
#'
#' @return A tibble with theoretical frequencies.
#'
freq_theor <- function(cn_bins) {
  total_n <- cn_bins |>
    filter(CN == 0) |>
    pull(n_reads_theoretical)
  if (length(total_n) == 0) {
    stop("No row with CN == 0 found in cn_bins. Cannot calculate theoretical frequencies.")
  }
  if (length(total_n) > 1) {
    stop("Multiple rows with CN == 0 found in cn_bins. Data integrity issue.")
  }
  if (is.na(total_n) || total_n <= 0) {
    stop(paste0(
      "n_reads_theoretical for CN == 0 must be a positive number, got: ", total_n
    ))
  }
  cn_bins_theor <- cn_bins |>
    filter(n_reads_theoretical != 0) |>
    mutate(freq_theoretical = n_reads_theoretical / total_n)
  return(cn_bins_theor)
}

#' Main execution function for frequency calculations
#'
#' @param cn_bins Data frame. Theoretical bins table.
#' @param bla_cn Data frame. Observed gene counts table.
#'
#' @return A combined tibble with all calculations.
#'
main <- function(cn_bins, bla_cn) {
  bla_cn_freq <- counts_freq_obs(bla_cn)
  bins_theor <- freq_theor(cn_bins)

  bla_cn_full <-
    bla_cn_freq %>%
    full_join(bins_theor, by = "CN") %>%
    # replace NA with 0
    mutate(
      counts_obs = replace_na(counts_obs, 0),
      freq_obs = replace_na(freq_obs, 0)
    ) %>%
    arrange(CN) %>%
    # do the rest of the calculations
    mutate(
      counts_corrected = counts_obs / freq_theoretical,
      freq_corrected = if (sum(counts_corrected, na.rm = TRUE) > 0) {
        counts_corrected / sum(counts_corrected, na.rm = TRUE)
      } else {
        NA_real_
      },
      detection_limit = 1 / n_reads_theoretical
    )

  return(bla_cn_full)
}

# --- 3. Snakemake Execution Block ---
if (exists("snakemake")) {
  log_file <- file(snakemake@log[[1]], open = "wt")
  sink(log_file, type = "output")
  sink(log_file, type = "message")
  on.exit({
    sink(type = "message")
    sink(type = "output")
    close(log_file)
  }, add = TRUE)

  cn_bins <- read_delim(snakemake@input[[1]], delim = "\t", show_col_types = FALSE)
  bla_cn <- read_delim(snakemake@input[[2]], delim = "\t", show_col_types = FALSE)

  output_table <- main(cn_bins, bla_cn)

  write_delim(output_table, file = snakemake@output[[1]], delim = "\t")
  print("Finished. No errors.")
}
