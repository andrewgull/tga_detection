#' @title Count merged blaSHV hits in filtered reads
#' @description This script processes merged intervals from a BED file and counts
#' hits against a filtered BLAST table. It calculates the number of gene hits based
#' on a provided gene length and joins this with the filtered results.
#' @section Input:
#' 1. BED file with merged intervals
#' 2. Filtered table of FR blast hits
#' 3. Path to snakemake params (gene length)
#' @section Output:
#' TSV file with blaSHV counts per subject
#' @md

# --- 1. Load Libraries ---
suppressPackageStartupMessages(library(dplyr))
library(readr)

# --- 2. Define Main Function ---
#' Count merged blaSHV hits in filtered reads
#'
#' @param df_merge Data frame/tibble. Table with merged intervals (typically from BED).
#' @param df_filt Data frame/tibble. Filtered table with flanking region (FR) blast hits.
#' @param len Numeric/Integer. Length of the target gene (e.g., blaSHV) used for normalization.
#'
#' @return A tibble with `subject` and `n.blaSHV.merged` columns.
#' @export
#'
main <- function(df_merge, df_filt, len) {
  # convert/check len
  len <- as.integer(len)
  stopifnot(!is.na(len))

  df_merge |>
    group_by(X1) |>
    summarise(sum.merged.hits = sum(X3 - X2 + 1)) |>
    mutate(n.blaSHV.merged = round(sum.merged.hits / len, 0)) |>
    rename("subject" = X1) |>
    # filter out the reads the 'bad reads' from previous filtering rounds
    right_join(df_filt, by = "subject") |>
    select(subject, n.blaSHV.merged)
}

# --- 3. Snakemake Execution Block ---
if (exists("snakemake")) {
  sink(snakemake@log[[1]])
  on.exit(sink(), add = TRUE)

  bla_merge <- read_tsv(snakemake@input[[1]],
    col_names = FALSE,
    show_col_types = FALSE
  )
  blast_filt <- read_tsv(snakemake@input[[2]],
    show_col_types = FALSE
  )

  bla_merged_df <- main(bla_merge, blast_filt, snakemake@params[[1]])

  # save results to file
  write_tsv(x = bla_merged_df, file = snakemake@output[[1]])
  print("Finished. No errors.")
}
