##############################################
# script to filter gene (blaSHV) blast hits
# keeping only those reads where blaSHV hits
# are located only within FRs
# and those that are located on reads filtered
# in the previous steps
# input: blast hits table
# input: filtered table RED+RU+filt+GREEN
# input: max e-value
# output: filtered table of blaSHV gene hits
##############################################

#### OPEN LOG ####
sink(snakemake@log[[1]])

#### LIBRARIES ####
library(readr)
library(purrr)
library(tidyr)
suppressPackageStartupMessages(library(dplyr))

#### FUNCTIONS ####
# Function to check whether a pair of coordinates (bla)
# is within another pair of coordinates (FRs)
# Returns TRUE if yes, FALSE if no
is_within <- function(bla_start, bla_end, fr_start, fr_end) {
  if (fr_start < fr_end) {
    # Direct strand: check if both ends of the gene
    # are inside the FR coordinates
    return(between(bla_start, fr_start, fr_end) &&
             between(bla_end, fr_start, fr_end))
  } else {
    # Reverse strand: check the same condition but with FR ends swapped
    return(between(bla_start, fr_end, fr_start) &&
             between(bla_end, fr_end, fr_start))
  }
}

# Function to check if any bla hit is located outside FRs for a given read id
# Returns FALSE if no bla hit is located outside FRs,
# TRUE if at least one hit is outside FRs
genes_within_flanks <- function(read_id, fr_df, bla_df) {
  fr_bla_df <- fr_df %>%
    filter(subject == read_id) %>%
    left_join(bla_df, by = "subject") %>%
    drop_na(start.subject.y) %>%
    select(subject, start.subject.x, end.subject.x,
           end.red, start.subject.y, end.subject.y, orientation)

  answer <- any(
    !map2_lgl(
      fr_bla_df$start.subject.y,
      fr_bla_df$end.subject.y,
      ~ is_within(.x, .y, fr_bla_df$end.red[1], fr_bla_df$end.subject.x[1])
    )
  )
  return(answer)
}

# Main function that ties everything together
main <- function(bla, fr, evalue) {
  # Find reads with bla hits outside FRs
  fr$aberrant <- map_vec(fr$subject, ~ genes_within_flanks(., fr, bla))

  # Keep non-aberrant (FALSE) reads only
  fr <- filter(fr, !aberrant)

  # Filter bla hits keeping only those present in FR table
  bla_filt <- bla %>%
    filter(subject %in% unique(fr$subject)) %>%
    filter(e.value <= evalue)

  names(bla_filt) <- c("query", "subject", "identity", "length", "mismatch",
                       "gaps", "start.query", "end.query", "start.subject",
                       "end.subject", "e.value", "bit.score")
  return(bla_filt)
}

#### RUN ####
tryCatch({
  bla_df <- read_delim(snakemake@input[[1]], show_col_types = FALSE)
  fr_df <- read_delim(snakemake@input[[2]], show_col_types = FALSE)

  bla_df_filt <- main(bla = bla_df, fr = fr_df, evalue = snakemake@params[[1]])

  # Save to file
  write_delim(bla_df_filt, file = snakemake@output[[1]], delim = "\t")
  print("Finished. No errors.")
}, error = function(e) {
  print(paste("Error: ", e$message))
})

#### CLOSE LOG ####
sink()
