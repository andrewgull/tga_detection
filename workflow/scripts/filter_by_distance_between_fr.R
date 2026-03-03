#' @title Filter reads by distance between flanking regions (FR)
#' @description This script calculates the expected number of gene copies (e.g., blaSHV)
#' based on the distance between flanking regions using the formula:
#' N = ((D_fr - L_fr) / L_ru) + 1, where:
#' - N: expected number of copies
#' - D_fr: distance between FRs
#' - L_fr: base length (distance with 1 copy)
#' - L_ru: length of the repeat unit (RU)
#' Reads are filtered if the difference between expected and observed (blast/merged)
#' counts is greater than 1.
#' @section Input:
#' 1. Table with distance between FRs
#' 2. Table with bla counts via merging BED coordinates
#' @section Output:
#' Filtered table with subject and n.blaSHV.merged
#' @md

# --- 1. Load Libraries ---
suppressPackageStartupMessages(library(dplyr))
library(readr)

# --- 2. Define Functions ---
#' Filter reads by distance between FRs
#'
#' @param fr_ru_filt Data frame/tibble. Table with distance between flanking regions.
#' @param bla_counts Data frame/tibble. Table with observed gene counts (merged BED).
#' @param base_len Numeric/Integer. Base length (expected distance for 1 copy). Default 4299.
#' @param ru_len Numeric/Integer. Length of the repeat unit. Default 3450.
#'
#' @return A tibble with `subject` and `n.blaSHV.merged` for reads passing the filter.
#' @export
#'
filter_by_distance <- function(fr_ru_filt, bla_counts,
                              base_len = 4299, ru_len = 3450) {
  # convert/check the base_len and ru_len
  base_len <- as.integer(base_len)
  stopifnot(!is.na(base_len))
  ru_len <- as.integer(ru_len)
  stopifnot(!is.na(ru_len), ru_len > 0)

  diff_per_read <- fr_ru_filt |>
    select(subject, distance.btw.FR) |>
    left_join(bla_counts, by = "subject") |>
    mutate(
      "n.bla.exp" = ((distance.btw.FR - base_len) / ru_len) + 1,
      "n.bla.blast" = n.blaSHV.merged,
      "difference" = round(n.bla.exp - n.blaSHV.merged)
    ) |>
    # filter out abs(diff) > 1
    filter(abs(difference) <= 1) |>
    select(subject, n.blaSHV.merged)

  return(diff_per_read)
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

  # read the input tables
  # this FR table contains distance between FRs
  fr_ru_filt <- read_delim(snakemake@input[[1]],
    show_col_types = FALSE, progress = FALSE
  )
  # this one contains bla counts via merging BED coords
  bla_counts <- read_delim(snakemake@input[[2]],
    show_col_types = FALSE, progress = FALSE
  )

  # filter
  filtered_df <- filter_by_distance(fr_ru_filt, bla_counts,
    base_len = snakemake@params[[1]],
    ru_len = snakemake@params[[2]]
  )

  # write the output
  write_delim(filtered_df, snakemake@output[[1]], delim = "\t")
  print("Finished. No errors.")
}
