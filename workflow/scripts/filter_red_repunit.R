#############################################################
# script to filter reads containing correct combination of
# the red flanking region and the repeat unit
# the 'correct combination' that is
# 'a red FR first and a repeat unit right next to it'
# input: blast hits table for RR;
# blast hits table for RU/IS
# max evalue for blast hits
# minimal length for RR blast hits
# minimal length for RU/IS hits
# minimal identity
# max allowed distance between FR and RU/IS
# output: filtered table of hits
#############################################################

#### OPEN LOG ####
sink(snakemake@log[[1]])

#### LIBRARIES ####
suppressPackageStartupMessages(library(dplyr))
library(readr)

#### FUNCTIONS ####
parse_blast <- function(file_path, region_name) {
  # read blast table
  df <-
    read_delim(file_path, col_names = FALSE, show_col_types = FALSE)
  # check if the table is empty
  stopifnot(nrow(df) > 0)
  # check if nrow is wrong
  stopifnot(ncol(df) == 12)
  # assign names
  names(df) <- c(
    "query",
    "subject",
    "identity",
    "length",
    "mismatch",
    "gaps",
    "start.query",
    "end.query",
    "start.subject",
    "end.subject",
    "e.value",
    "bit.score"
  )
  # create orientation column
  df$orientation <-
    ifelse(df$start.subject < df$end.subject, "direct", "reverse")
  # rename query
  df$query <- region_name
  return(df)
}

filter_blast <- function(blast_df, min_len, max_e_value, min_identity) {
  # return: a blast table filtered by minimal length, evalue and hit identity
  filt_df <- filter(blast_df,
                    length >= min_len,
                    e.value <= max_e_value,
                    identity >= min_identity)
  return(filt_df)
}

filter_ru_fr <-
  function(fr_df, ru_df, max_distance) {
    # filters red blast and rep unit blast by orientation and max
    # distance btw the red region and the rep.unit
    # return: filtered tibble
    left_join(fr_df, ru_df, by = "subject") %>%
      filter(orientation.x == orientation.y) %>%
      mutate(
        distance = if_else(
          orientation.x == "reverse",
          end.subject.y - start.subject.x + 1,
          start.subject.x - end.subject.y + 1
        )
      ) %>%
      filter(abs(distance) <= max_distance)
  }

main <- function(red, rep, ru_len, fr_len, e, identity, maxd) {
  # filter rep unite table
  # NB: min length is different from FR's min length
  blast_red <- parse_blast(red, "FR_red") %>%
    filter_blast(min_len = fr_len,
                 max_e_value = e,
                 min_identity = identity)
  blast_rep <- parse_blast(rep, "Rep_unit") %>%
    filter_blast(min_len = ru_len,
                 max_e_value = e,
                 min_identity = identity)
  # filter combination of both FR and rep unit
  # 1st: same orientation
  # 2nd: close to each other
  blast_joined <-
    filter_ru_fr(blast_red, blast_rep, max_distance = maxd) %>%
    # keep some columns
    select(subject, start.subject.x, end.subject.x,
           start.subject.y, end.subject.y, distance, orientation.x)
  # give them better names
  names(blast_joined) <- c("subject", "start.red", "end.red",
                           "start.rep.unit", "end.rep.unit", "dist", "orient")
  return(blast_joined)
}

#### RUN ####
output_table <- main(red = snakemake@input[[1]],
                     rep = snakemake@input[[2]],
                     identity = snakemake@params[[1]],
                     e = snakemake@params[[2]],
                     fr_len = snakemake@params[[3]],
                     ru_len = snakemak@params[[4]],
                     maxd = snakemake@params[[5]])

# Save to file
write_delim(output_table, file = snakemake@output[[1]], delim = "\t")
print("Finished. No erorrs.")

#### CLOSE LOG ####
sink()
