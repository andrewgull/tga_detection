##############################################
# script to filter gene (blaSHV) blast hits
# keeping only those reads where blaSHV hits
# are located only within FRs
# and those that are located on reads filtered
# in the previous steps
# input: blast hits table
# input: filtered table RED+RU+filt+GREEN
# input: max e-vlaue
# output: filtered table of blaSHV gene hits
##############################################

#### OPEN LOG ####
sink(snakemake@log[[1]])

#### LIBRARIES ####
library(ggplot2)
library(readr)
library(purrr)
suppressPackageStartupMessages(library(dplyr))

#### FUNCTIONS ####
# function to check whether a pair of coords (bla)
# is within another pair of coords (FRs)
# returns TRUE if yes, FALSE if no
is_within <-
  function(bla_start,
           bla_end,
           fr_start,
           fr_end) {
    # this is direct strand
    if (fr_start < fr_end) {
      # check if both ends of the gene are inside the FR coords
      if (between(bla_start, fr_start, fr_end) &&
            between(bla_end, fr_start, fr_end)) {
        return(TRUE)
      } else {
        return(FALSE)
      }
      # this is reverse strand
    } else {
      # check the same condition but with FR ends swapped
      if (between(bla_start, fr_end, fr_start) &&
            between(bla_end, fr_end, fr_start)) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
  }

# this function can be applied to each read id in a table
# input: read id
# input: table with FR hits
# input: table with bla hits
# return: FALSE if no bla hit is located outside FRs
# TRUE if at least one hit is outside FRs
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
  # FR coords are always the first elements of end.red, end.subject.x columns
  answer <- any(
    !map2_lgl(
      fr_bla_df$start.subject.y,
      fr_bla_df$end.subject.y,
      ~ is_between(.x, .y, fr_bla_df$end.red[1], fr_bla_df$end.subject.x[1])
    )
  )
  return(answer)
}

# function that ties everything together
main <- function(bla, fr, evalue) {
  # find reads with bal hits outside FRs
  fr$aberrant <- map_vec(fr$subject, ~ genes_within_flanks(., fr, bla))

  # keep non-aberrant (FALSE) reads only
  fr <- filter(fr, !aberrant)

  # filter bla hits keeping only those present in FR table
  bla_filt <-
    bla %>%
    # leave only the reads in filtered FR blast
    filter(X2 %in% unique(fr$subject)) %>%
    # remove unreliable hits
    filter(X11 <= evalue)

  names(bla_filt) <- c("query", "subject", "identity", "length", "mismatch",
                       "gaps", "start.query", "end.query", "start.subject",
                       "end.subject", "e.value", "bit.score")
  return(bla_filt)
}

#### RUN ####
bla_df <- read_delim(snakemake@input[[1]],
                     col_names = FALSE,
                     show_col_types = FALSE)
fr_df <- read_delim(snakemake@input[[2]],
                    show_col_types = FALSE)

bla_df_filt <- main(bla = bla_df, fr = fr_df, evalue = snakemake@params[[1]])

# save to file
write_delim(bla_df_filt, file = snakemake@output[[1]], delim = "\t")
print("Finished. No erorrs.")

#### CLOSE LOG ####
sink()
