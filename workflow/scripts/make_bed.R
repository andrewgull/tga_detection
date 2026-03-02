#' @title Make BED file from BLAST outfmt 6
#' @description This script converts a BLAST output table (tabular format 6) into
#' a standard 6-column BED file. it handles both filtered (with headers) and
#' non-filtered (without headers) BLAST tables.
#' @section Input:
#' BLAST output table in format 6 (tabular).
#' @section Output:
#' BED file with columns: chrom, start, end, name, score, strand.
#' @md

# --- 1. Load Libraries ---
suppressPackageStartupMessages(library(dplyr))
library(readr)

# --- 2. Define Functions ---
#' Check if a BLAST Table is a Filtered Version
#'
#' Determines whether a BLAST table is a "filtered" version based on its filename.
#' Filtered versions are expected to contain headers, whereas non-filtered versions are not.
#'
#' @param filename Character. The path or name of the input BLAST table file.
#'
#' @return Logical. Returns `TRUE` if "filtered" is found in the filename, `FALSE` otherwise.
#'
is_filtered <- function(filename) {
  # if the blast table has headers it's been filtred
  # otherwise it's not
  filtered <- grepl("filtered", filename, fixed = TRUE)
  if (filtered) {
    # message to log file
    print("filtred version of the blast table is detected")
    return(TRUE)
  } else {
    # message to log file
    print("Non-filtred version of the blast table is detected")
    return(FALSE)
  }
}

#' Read a BLAST outfmt 6 Table
#'
#' Reads a BLAST output file in format 6 (tabular). It automatically detects if the table
#' is a filtered version (with headers) or a standard version (without headers) using [is_filtered()]
#' and ensures the resulting data frame has consistent, standard column names.
#'
#' @param input_blast Character. Path to the input BLAST table file.
#'
#' @return A tibble with standard BLAST outfmt 6 columns: `query`, `subject`, `identity`,
#'   `length`, `mismatch`, `gaps`, `start.query`, `end.query`, `start.subject`,
#'   `end.subject`, `e.value`, and `bit.score`.
#'
read_blast <- function(input_blast) {
  if (is_filtered(input_blast)) {
    blast_bla <- read_tsv(input_blast,
      show_col_types = FALSE,
      col_names = TRUE
    )
    # check if the correct variant of the blast table was read
    # incorrect version of the table (without headers)
    # will have 'X1' as the first column name
    stopifnot(grepl("query", names(blast_bla)[1], fixed = TRUE))
  } else {
    blast_bla <- read_tsv(input_blast,
      show_col_types = FALSE,
      col_names = FALSE
    )
    # check if the correct variant of the blast table was read
    # incorrect version of the table (with headers)
    # will have 'query' and 'blaSHV-1' in X1
    stopifnot(length(unique(blast_bla$X1)) == 1)
    # check dimensions
    stopifnot(ncol(blast_bla) == 12)
    # set names
    names(blast_bla) <- c(
      "query", "subject", "identity",
      "length", "mismatch", "gaps",
      "start.query", "end.query", "start.subject",
      "end.subject", "e.value", "bit.score"
    )
  }
  # check if the table is empty
  stopifnot(nrow(blast_bla) > 0)
  return(blast_bla)
}

#' Convert BLAST Results to BED Format
#'
#' Transforms a BLAST results data frame into a standard 6-column BED format.
#' The function calculates 0-indexed coordinates, determines the strand based on
#' start/end positions, and ensures that the `start` coordinate is always less than
#' the `end` coordinate, as required by the BED specification.
#'
#' @param blast Data frame. A BLAST results table, typically the output of [read_blast()].
#'
#' @return A tibble in BED format with columns: `chrom`, `start`, `end`, `name`, `score`, and `strand`.
#'
blast2bed <- function(blast) {
  # convert to BED format
  bla_bed <-
    blast |>
    select(subject, start.subject, end.subject, query, e.value) |>
    mutate(
      strand = if_else(start.subject < end.subject, "+", "-"),
      start = pmin(start.subject, end.subject) - 1,
      end = pmax(start.subject, end.subject)
    ) |>
    rename(
      "chrom" = subject,
      "name" = query,
      "score" = e.value
    ) |>
    select(chrom, start, end, name, score, strand)

  return(bla_bed)
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

  # read
  blast_tab <- read_blast(snakemake@input[[1]])
  # convert
  bed_tab <- blast2bed(blast_tab)
  # save to file
  write_tsv(bed_tab, snakemake@output[[1]], col_names = FALSE)
  print("Finished. No errors.")
}
