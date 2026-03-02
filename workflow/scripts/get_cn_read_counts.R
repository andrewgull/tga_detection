#' @title Calculate theoretical read counts per copy number (CN)
#' @description This script calculates the number of reads that are theoretically
#' long enough to contain each possible gene copy number, based on flanking
#' region (FR) and repeat unit (RU) positions and read lengths.
#' @section Input:
#' 1. Table with BLAST results filtered by FR+RU orientation and distance.
#' 2. Filtered FASTA reads (to get lengths).
#' 3. Snakemake parameters: maximum CN, length increment, and base length.
#' @section Output:
#' Table with CNs and theoretical read counts.
#' @md

# --- 1. Load Libraries ---
suppressPackageStartupMessages(library(dplyr))
library(readr)
library(purrr)
library(tibble)

# --- 2. Define Functions ---
#' Count reads capable of containing a specific CN
#'
#' @param df Data frame/tibble. BLAST results with end.red and optionally read.len.
#' @param cn Integer. Gene copy number.
#' @param b Integer. Base length (min distance for CN 0).
#' @param i Integer. Length increment per copy.
#' @param direct Logical. TRUE if hits are in direct orientation.
#'
#' @return Integer count of reads.
#'
count_reads_cn <- function(df, cn = 0, b = 1320, i = 3450, direct = TRUE) {
  stopifnot(sum(grepl("end.red", names(df), fixed = TRUE)) == 1)
  if (direct) {
    df |>
      mutate(keep = end.red > (b + i * cn)) |>
      filter(keep) |>
      nrow()
  } else {
    stopifnot(sum(grepl("read.len", names(df), fixed = TRUE)) == 1)
    df |>
      mutate(keep = (read.len - end.red) > (b + i * cn)) |>
      filter(keep) |>
      nrow()
  }
}

#' Get read lengths from FASTA
#'
#' @param fasta Character. Path to the FASTA file.
#'
#' @return A tibble with `subject` and `read.len`.
#'
read_length <- function(fasta) {
  lengths <- Biostrings::fasta.seqlengths(fasta)
  len_df <- tibble(
    "subject" = sub("^(.*?) runid=.*", "\\1", names(lengths)),
    "read.len" = lengths
  )
  return(len_df)
}

#' Separate direct from reverse hits
#'
#' @param df Data frame/tibble. Filtered hits.
#' @param orientation Character. "direct" or "reverse".
#' @param reads_len Data frame/tibble. Read lengths (required for reverse).
#'
#' @return A separated tibble.
#'
separate <- function(df, orientation, reads_len = NULL) {
  stopifnot(sum(grepl("orient", names(df), fixed = TRUE)) == 1)
  stopifnot(sum(grepl("subject", names(df), fixed = TRUE)) == 1)
  if (orientation == "direct") {
    df |> filter(orient == "direct")
  } else {
    stopifnot(!is.null(reads_len))
    stopifnot(sum(grepl("subject", names(reads_len), fixed = TRUE)) == 1)
    df |>
      filter(orient == "reverse") |>
      left_join(reads_len, by = "subject")
  }
}

#' Count reads for each CN variant
#'
#' @param cn_array Integer vector. Array of CN values.
#' @param hits_orient Data frame/tibble. Hits separated by orientation.
#' @param min_len Integer. Base minimal length.
#' @param increment Integer. Length increment per copy.
#' @param direct Logical. Orientation.
#'
#' @return Integer vector of counts.
#'
count_reads <- function(cn_array, hits_orient, min_len, increment, direct) {
  map_int(cn_array, ~ count_reads_cn(
    hits_orient,
    cn = .,
    b = min_len,
    i = increment,
    direct = direct
  ))
}

#' Main execution function for theoretical CN counts
#'
#' @param in_table Character. Path to hits table.
#' @param reads Character. Path to FASTA reads.
#' @param max_cn Integer. Maximum CN variant.
#' @param min_len Integer. Base length.
#' @param increment Integer. RU increment.
#'
#' @return A tibble with CN and counts.
#'
main <- function(in_table, reads, max_cn, min_len, increment) {
  max_cn <- as.integer(max_cn)
  stopifnot(!is.na(max_cn))

  min_len <- as.integer(min_len)
  stopifnot(!is.na(min_len))

  increment <- as.integer(increment)
  stopifnot(!is.na(increment))

  fr_repunit <- read.table(in_table, sep = "\t", header = TRUE)
  stopifnot(ncol(fr_repunit) == 7)
  stopifnot(nrow(fr_repunit) != 0)

  reads_len_df <- read_length(reads)

  hits_direct <- separate(fr_repunit, "direct")
  hits_reverse <- separate(fr_repunit, "reverse", reads_len = reads_len_df)

  cnv_variants <- seq(0, max_cn, 1)

  counts_direct <- count_reads(cnv_variants, hits_direct,
    min_len = min_len, increment = increment,
    direct = TRUE
  )

  counts_reverse <- count_reads(cnv_variants, hits_reverse,
    min_len = min_len, increment = increment,
    direct = FALSE
  )

  n_reads_cn_tot <- counts_direct + counts_reverse

  output_table <- data.frame(
    "CN" = cnv_variants,
    "n_reads_theoretical" = n_reads_cn_tot
  )
  return(output_table)
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

  cnv_theoretical <- main(
    in_table = snakemake@input[[1]],
    reads = snakemake@input[[2]],
    max_cn = snakemake@params[[1]],
    increment = snakemake@params[[2]],
    min_len = snakemake@params[[3]]
  )

  write.table(cnv_theoretical,
    file = snakemake@output[[1]],
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  print("Finished. No errors.")
}
