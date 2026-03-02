#' @title Filter reads by distance to end
#' @description This script filters reads that have both flanking region (FR) and
#' repeat unit (RU) in the correct orientation and proximity, and are at least a
#' minimum distance from the beginning of the FR.
#' @section Input:
#' 1. Table with BLAST results filtered by FR+RU orientation and distance.
#' 2. Minimum distance parameter (e.g., 1320 nt).
#' @section Output:
#' Table with reads meeting the minimum distance requirement.
#' @md

# --- 1. Load Libraries ---
suppressPackageStartupMessages(library(dplyr))
library(readr)

# --- 2. Define Functions ---
#' Filter reads by distance from start of FR towards RU
#'
#' @param df Data frame/tibble. Table with BLAST results and orientation info.
#' @param dist Numeric/Integer. Minimum distance from the start of FR towards RU.
#'
#' @return A filtered tibble.
#'
filter_by_distance <- function(df, dist = 1320) {
  df |>
    mutate(keep = if_else(
      orient == "direct",
      end.red > dist,
      start.rep.unit > dist
    )) |>
    filter(keep) |>
    select(-keep)
}

#' Main execution function for filtering
#'
#' @param df Character. Path to the input TSV table.
#' @param dist Numeric/Integer. Minimum distance value.
#'
#' @return A tibble with filtered results.
#'
main <- function(df, dist) {
  # read the input table
  input_table <- read_tsv(df, show_col_types = FALSE)

  # make size checks
  stopifnot(ncol(input_table) == 7)
  stopifnot(nrow(input_table) != 0)

  # convert/check dist
  dist <- as.integer(dist)
  stopifnot(!is.na(dist))

  # filter
  filtered_df <- filter_by_distance(input_table, dist = dist)
  return(filtered_df)
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

  output_table <- main(
    df = snakemake@input[[1]],
    dist = snakemake@params[[1]]
  )

  # save to file
  write_delim(output_table, file = snakemake@output[[1]], delim = "\t")
  print("Finished. No errors.")
}
