#' @title Filter reads by RED flanking region and repeat unit (RU)
#' @description This script filters reads containing the correct combination of
#' a RED flanking region followed by a repeat unit. It checks orientation,
#' distance, and basic BLAST metrics.
#' @section Input:
#' 1. BLAST table for RED FR.
#' 2. BLAST table for Rep_unit (RU/IS).
#' 3. Snakemake parameters: identity, e-value, FR min length, RU min length, max distance.
#' @section Output:
#' Filtered table with joined RED and RU coordinates.
#' @md

# --- 1. Load Libraries ---
suppressPackageStartupMessages(library(dplyr))
library(readr)

# --- 2. Define Functions ---
#' Parse BLAST outfmt 6 table
#'
#' @param file_path Character. Path to the BLAST output file.
#' @param region_name Character. Name to assign to the query column.
#'
#' @return A tibble with standard BLAST columns and orientation.
#'
parse_blast <- function(file_path, region_name) {
  df <- read_delim(file_path, col_names = FALSE, show_col_types = FALSE)
  stopifnot(nrow(df) > 0)
  stopifnot(ncol(df) == 12)

  names(df) <- c(
    "query", "subject", "identity", "length", "mismatch",
    "gaps", "start.query", "end.query", "start.subject",
    "end.subject", "e.value", "bit.score"
  )
  df$orientation <- ifelse(df$start.subject < df$end.subject, "direct", "reverse")
  df$query <- region_name
  return(df)
}

#' Filter BLAST results by basic metrics
#'
#' @param blast_df Data frame/tibble. BLAST results.
#' @param min_len Integer. Minimal hit length.
#' @param max_e_value Numeric. Maximum e-value.
#' @param min_identity Integer. Minimal identity percentage.
#'
#' @return A filtered tibble.
#'
filter_blast <- function(blast_df, min_len, max_e_value, min_identity) {
  filt_df <- filter(
    blast_df,
    length >= min_len,
    e.value <= max_e_value,
    identity >= min_identity
  )
  return(filt_df)
}

#' Filter joined FR and RU by orientation and distance
#'
#' @param fr_df Data frame/tibble. RED FR table.
#' @param ru_df Data frame/tibble. RU/IS table.
#' @param max_distance Integer. Maximum distance between FR and RU.
#'
#' @return A joined and filtered tibble.
#'
filter_ru_fr <- function(fr_df, ru_df, max_distance) {
  left_join(fr_df, ru_df, by = "subject") |>
    filter(orientation.x == orientation.y) |>
    mutate(
      distance = if_else(
        orientation.x == "reverse",
        end.subject.y - start.subject.x + 1,
        start.subject.x - end.subject.y + 1
      )
    ) |>
    filter(abs(distance) <= max_distance)
}

#' Main execution function for filtering RED+RU
#'
#' @param red Character. Path to RED BLAST file.
#' @param rep Character. Path to RU BLAST file.
#' @param ru_len Integer. RU min length.
#' @param fr_len Integer. FR min length.
#' @param e Numeric. Max e-value.
#' @param identity Integer. Min identity.
#' @param maxd Integer. Max distance.
#'
#' @return A tibble with filtered and renamed columns.
#'
main <- function(red, rep, ru_len, fr_len, e, identity, maxd) {
  ru_len <- as.integer(ru_len)
  stopifnot(!is.na(ru_len))

  fr_len <- as.integer(fr_len)
  stopifnot(!is.na(fr_len))

  e <- as.numeric(e)
  stopifnot(!is.na(e))

  identity <- as.integer(identity)
  stopifnot(!is.na(identity))

  maxd <- as.integer(maxd)
  stopifnot(!is.na(maxd))

  blast_red <- parse_blast(red, "FR_red") |>
    filter_blast(
      min_len = fr_len,
      max_e_value = e,
      min_identity = identity
    )
  blast_rep <- parse_blast(rep, "Rep_unit") |>
    filter_blast(
      min_len = ru_len,
      max_e_value = e,
      min_identity = identity
    )

  blast_joined <- filter_ru_fr(blast_red, blast_rep, max_distance = maxd) |>
    select(
      subject, start.subject.x, end.subject.x,
      start.subject.y, end.subject.y, distance, orientation.x
    )

  names(blast_joined) <- c(
    "subject", "start.red", "end.red",
    "start.rep.unit", "end.rep.unit", "dist", "orient"
  )
  return(blast_joined)
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
    red = snakemake@input[[1]],
    rep = snakemake@input[[2]],
    identity = snakemake@params[[1]],
    e = snakemake@params[[2]],
    fr_len = snakemake@params[[3]],
    ru_len = snakemake@params[[4]],
    maxd = snakemake@params[[5]]
  )

  write_delim(output_table, file = snakemake@output[[1]], delim = "\t")
  print("Finished. No errors.")
}
