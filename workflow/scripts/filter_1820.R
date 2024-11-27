#############################################################################
# script to filter reads with both FR and RU in correct orientation
# and proximity AND at least 1820 nt from the beginning of the FR
# input: table with blast results filtered by FR+RU orientation and distance
# input: min distance from the start of FR towards RU (1820)
# output: table with reads at least 1820 nt in length from the start of FR
#############################################################################

#### OPEN LOG ####
sink(snakemake@log[[1]])

#### LIBRARIES ####
suppressPackageStartupMessages(library(dplyr))
library(readr)

#### FUNCTIONS ####
filter_by_distance <- function(df, dist = 1820) {
  # df: a table with 7 columns
  # colnames must be as follows
  df %>%
    mutate(keep = if_else(
      orient == "direct",
      end.red > dist,
      start.rep.unit > dist
    )) %>%
    filter(keep) %>%
    select(-keep)
}

main <- function(df, dist) {
  # df: input table
  # dist: input distance
  # return: output table
  # read the input table
  input_table <- read_tsv(df, show_col_types = FALSE)

  # make size checks
  stopifnot(ncol(input_table) == 7)
  stopifnot(nrow(input_table) != 0)

  # filter
  filtered_df <- filter_by_distance(input_table, dist = dist)
  return(filtered_df)
}

#### RUN ####
output_table <- main(df = snakemake@input[[1]],
                     dist = snakemake@params[[1]])

# save to file
write_delim(output_table, file = snakemake@output[[1]], delim = "\t")
print("Finished. No erorrs.")

#### CLOSE LOG ####
sink()
