#' @title Filter flanking region (FR) BLAST hits
#' @description This script filters FR BLAST hits by checking for the presence of
#' both RED and GREEN flanking regions in the correct orientation. It joins the
#' filtered results and calculates the distance between FRs.
#' @section Input:
#' 1. Table of RED hits (after previous filtering).
#' 2. Table of GREEN hits (raw BLAST).
#' 3. Snakemake parameters: identity, e-value, and minimal length.
#' @section Output:
#' Filtered and joined table of RED+GREEN hits.
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
filter_blast_p1 <- function(blast_df, min_len, max_e_value, min_identity) {
  filt_df <- filter(
    blast_df,
    length >= min_len,
    e.value <= max_e_value,
    identity >= min_identity
  )
  return(filt_df)
}

#' Identify subjects with multiple hits
#'
#' @param df Data frame/tibble. BLAST results.
#'
#' @return Character vector of subject IDs.
#'
get_multihits_ids <- function(df) {
  df |>
    group_by(subject) |>
    count() |>
    filter(n > 1) |>
    pull(subject)
}

#' Join and filter RED and GREEN flanking regions
#'
#' @param df_red_ru Data frame/tibble. RED hits table.
#' @param df_green Data frame/tibble. GREEN hits table.
#'
#' @return A joined and filtered tibble.
#'
filter_red_green <- function(df_red_ru, df_green) {
  reads_multiple_hits <- union(
    get_multihits_ids(df_red_ru),
    get_multihits_ids(df_green)
  )

  df_red_ru_filt <- df_red_ru |>
    filter(!subject %in% reads_multiple_hits)
  df_green_filt <- df_green |>
    filter(!subject %in% reads_multiple_hits)

  df_joined <- full_join(df_red_ru_filt, df_green_filt, by = "subject") |>
    filter(!is.na(query.x), !is.na(query.y), orient == orientation) |>
    mutate(
      green.red.distance = end.red - start.subject,
      distance.btw.FR = abs(green.red.distance)
    )
  return(df_joined)
}

#' Main execution function for filtering FR hits
#'
#' @param red Data frame. RED hits.
#' @param green Data frame. GREEN hits.
#' @param len Integer. Minimal length.
#' @param evalue Numeric. Maximum e-value.
#' @param identity Integer. Minimal identity.
#'
#' @return A tibble with filtered and joined hits.
#'
main <- function(red, green, len, evalue, identity) {
  len <- as.integer(len)
  stopifnot(!is.na(len))

  evalue <- as.numeric(evalue)
  stopifnot(!is.na(evalue))

  identity <- as.integer(identity)
  stopifnot(!is.na(identity))

  red$query <- "FR_RU_filt"
  green_filt <- green |>
    filter_blast_p1(
      min_len = len,
      max_e_value = evalue,
      min_identity = identity
    )

  stopifnot(nrow(red) > 0)
  stopifnot(nrow(green_filt) > 0)

  blast_joined <- filter_red_green(red, green_filt)
  return(blast_joined)
}

# --- 3. Snakemake Execution Block ---
if (exists("snakemake")) {
  sink(snakemake@log[[1]])
  on.exit(sink(), add = TRUE)

  blast_red <- read_tsv(snakemake@input[[1]], show_col_types = FALSE)
  blast_green <- parse_blast(snakemake@input[[2]], "FR_green")

  output_filtered <- main(
    red = blast_red,
    green = blast_green,
    identity = snakemake@params[[1]],
    evalue = snakemake@params[[2]],
    len = snakemake@params[[3]]
  )

  write_delim(output_filtered, file = snakemake@output[[1]], delim = "\t")
  print("Finished. No errors.")
}
