#' @title Filter gene (blaSHV) BLAST hits
#' @description This script filters gene BLAST hits to keep only reads where hits
#' are located within flanking regions (FRs). It also filters based on reads
#' that passed previous steps.
#' @section Input:
#' 1. Table with gene (blaSHV) BLAST hits.
#' 2. Filtered table of RED+GREEN FR metadata.
#' 3. Snakemake parameter: maximum e-value.
#' @section Output:
#' Filtered table of gene hits.
#' @md

# --- 1. Load Libraries ---
library(readr)
library(purrr)
library(tidyr)
suppressPackageStartupMessages(library(dplyr))

# --- 2. Define Functions ---
#' Check if a hit is within coordinates
#'
#' @param bla_start Integer. Start coordinate of the gene hit.
#' @param bla_end Integer. End coordinate of the gene hit.
#' @param fr_start Integer. Start coordinate of the flanking region.
#' @param fr_end Integer. End coordinate of the flanking region.
#'
#' @return Logical. TRUE if hit is within coordinates.
#'
is_within <- function(bla_start, bla_end, fr_start, fr_end) {
  if (fr_start < fr_end) {
    between(bla_start, fr_start, fr_end) &&
      between(bla_end, fr_start, fr_end)
  } else {
    between(bla_start, fr_end, fr_start) &&
      between(bla_end, fr_end, fr_start)
  }
}

#' Check if gene hits are within flanks for a read
#'
#' @param read_id Character. Subject/Read ID.
#' @param fr_df Data frame/tibble. FR metadata.
#' @param bla_df Data frame/tibble. Gene BLAST hits.
#'
#' @return Logical. TRUE if any gene hit is OUTSIDE flanks.
#'
genes_within_flanks <- function(read_id, fr_df, bla_df) {
  fr_bla_df <- fr_df |>
    filter(subject == read_id) |>
    left_join(bla_df, by = "subject") |>
    drop_na(start.subject.y) |>
    select(
      subject, start.subject.x, end.subject.x,
      end.red, start.subject.y, end.subject.y, orientation
    )

  answer <- any(
    !map2_lgl(
      fr_bla_df$start.subject.y,
      fr_bla_df$end.subject.y,
      ~ is_within(.x, .y, fr_bla_df$end.red[1], fr_bla_df$end.subject.x[1])
    )
  )
  return(answer)
}

#' Main execution function for filtering gene hits
#'
#' @param bla Data frame. Gene BLAST hits.
#' @param fr Data frame. FR metadata.
#' @param evalue Numeric. Maximum e-value.
#'
#' @return A filtered tibble of gene hits.
#'
main <- function(bla, fr, evalue) {
  evalue <- as.numeric(evalue)
  stopifnot(!is.na(evalue))

  fr$aberrant <- map_vec(fr$subject, ~ genes_within_flanks(., fr, bla))
  fr <- filter(fr, !aberrant)

  bla_filt <- bla |>
    filter(subject %in% unique(fr$subject)) |>
    filter(e.value <= evalue)

  names(bla_filt) <- c(
    "query", "subject", "identity", "length", "mismatch",
    "gaps", "start.query", "end.query", "start.subject",
    "end.subject", "e.value", "bit.score"
  )
  return(bla_filt)
}

# --- 3. Snakemake Execution Block ---
if (exists("snakemake")) {
  sink(snakemake@log[[1]])
  on.exit(sink(), add = TRUE)

  bla_df <- read_delim(snakemake@input[[1]],
    show_col_types = FALSE,
    col_names = FALSE
  )
  names(bla_df) <- c(
    "query", "subject", "identity", "length", "mismatch",
    "gaps", "start.query", "end.query", "start.subject",
    "end.subject", "e.value", "bit.score"
  )
  fr_df <- read_delim(snakemake@input[[2]], show_col_types = FALSE)

  bla_df_filt <- main(bla = bla_df, fr = fr_df, evalue = snakemake@params[[1]])

  write_delim(bla_df_filt, file = snakemake@output[[1]], delim = "\t")
  print("Finished. No errors.")
}
