##########################################
# script to filter away reads with
# at least one blaSHV hit outside the FRs
# input 1:
# input 2:
# output:
##########################################

#### OPEN LOG ####
sink(snakemake@log[[1]])

#### LIBRARIES ####
suppressPackageStartupMessages(library(dplyr))
library(readr)

#### FUNCTIONS ####
# function to show if (gene) coords are within (FR) range
is_between <- function(bla_start, bla_end, fr_start, fr_end) {
    # this is direct strand
    if (fr_start < fr_end) {
      # check if both ends of the gene are inside the FR coords
      if (between(bla_start, fr_start, fr_end) &
          between(bla_end, fr_start, fr_end)) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    # this is reverse strand
    } else {
      # check the same condition but with FR ends swapped
      if (between(bla_start, fr_end, fr_start) &
          between(bla_end, fr_end, fr_start)) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
  }

# solution that can be applied to the full table
genes_within_flanks <- function(read_id, fr_df, bla_df) {
  fr_bla_df <- fr_df %>%
    filter(subject == read_id) %>%
    # this will copy values of FR coords across rows
    left_join(bla_df, by = "subject") %>%
    # drop na in .y -> no bla hits in this read
    drop_na(start.subject.y) %>% 
    # BLA columns are prefixed by .y
    # start.subject.x and end.subject.x are the same across all rows now
    select(
      subject,
      start.subject.x,
      end.subject.x,
      end.red,
      start.subject.y,
      end.subject.y,
      orientation
    )
  
  # FR coords are always the first elements of the end.red, end.subject.x columns
  answer <- any(
    !map2_lgl(
      fr_bla_df$start.subject.y,
      fr_bla_df$end.subject.y,
      ~ is_between(.x, .y, fr_bla_df$end.red[1], fr_bla_df$end.subject.x[1])
    )
  )
  
  return(answer)
}

# get the actual count of aberrant reads
count_aberrant_reads <- function(fr_ru_df, bla_df) {
  sum(map_vec(fr_ru_df$subject, ~ genes_within_flanks(., fr_ru_df, bla_df)))
}

# function that ties everything together
# it outputs a table with percentage of aberrant reads
# but I need TO FILTER THEM AWAY
main <- function(sample, prefix = "../../results/tables/") {
  fr_ru_path <- paste0(prefix, sample, "/blast_joined.tsv")
  bla_path <- paste0(prefix, sample, "/blast_blaSHV_filtered.tsv")
  fr_df <- read_delim(fr_ru_path, show_col_types = F)
  bla_df <- read_delim(bla_path, show_col_types = F)
  n_abb <- count_aberrant_reads(fr_df, bla_df)
  n_tot <- nrow(fr_df)
  
  out_df <- tibble("sample" = sample, "n_aberrant" = n_abb, "n_total" = n_tot)
  return(out_df)
}

