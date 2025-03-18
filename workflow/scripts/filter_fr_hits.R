############################################################
# script to filter FR blast hits by presence of
# both RED and GREEN in correct orientation
# input: table of RR hits after Red+RU+1820 filtering
# input: table of GR hits
# input: max e-value
# input: minimal hit length
# input: minimal hit identity
# output: filtered and joined table of GR+RR+RU+1820 hits
############################################################

#### OPEN LOG ####
sink(snakemake@log[[1]])

#### LIBRARIES ####
suppressPackageStartupMessages(library(dplyr))
library(readr)

#### FUNCTIONS ####
parse_blast <- function(file_path, region_name) {
  # read blast table
  df <- read_delim(file_path, col_names = FALSE, show_col_types = FALSE)
  # check if the table is empty
  stopifnot(nrow(df) > 0)
  # check if nrow is wrong
  stopifnot(ncol(df) == 12)
  # assign names
  names(df) <- c(
    "query", "subject", "identity", "length", "mismatch",
    "gaps", "start.query", "end.query", "start.subject",
    "end.subject", "e.value", "bit.score"
  )
  # create orientation column
  df$orientation <- ifelse(df$start.subject < df$end.subject, "direct", "reverse") # nolint: line_length_linter.
  # rename query
  df$query <- region_name
  return(df)
}

filter_blast_p1 <- function(blast_df, min_len, max_e_value, min_identity) {
  # return: a blast table filtered by minimal length, evalue and hit identity
  # apply only to GREEN hits table
  filt_df <- filter(
    blast_df,
    length >= min_len,
    e.value <= max_e_value,
    identity >= min_identity
  )
  filt_df
}

get_multihits_ids <- function(df) {
  # df: filtered blast table
  # return: subject (read ids) containing multiple query hits
  df %>%
    group_by(subject) %>%
    count() %>%
    filter(n > 1) %>%
    pull(subject)
}

filter_red_green <- function(df_red_ru, df_green) {
  # find read IDs containing multiple query hits
  reads_multiple_hits <- union(
    get_multihits_ids(df_red_ru),
    get_multihits_ids(df_green)
  )
  # filter out reads with multiple query hits
  # df_red_ru have differen column names!
  df_red_ru_filt <- df_red_ru %>%
    filter(!subject %in% reads_multiple_hits)
  df_green_filt <- df_green %>%
    filter(!subject %in% reads_multiple_hits)
  # Filter by orientation
  # Add distance between FRs
  df_joined <- full_join(df_red_ru_filt, df_green_filt, by = "subject") %>%
    filter(!is.na(query.x), !is.na(query.y), orient == orientation) %>%
    mutate(
      green.red.distance = end.red - start.subject,
      distance.btw.FR = if_else(green.red.distance < 0,
        green.red.distance * -1,
        green.red.distance * 1
      )
    )
  df_joined
}

main <- function(red, green, len, evalue, identity) {
  # add readable query name
  red$query <- "FR_RU_filt"
  green_filt <- green %>%
    # apply basic filtering of blast results
    filter_blast_p1(
      min_len = len,
      max_e_value = evalue,
      min_identity = identity
    )
  # check if they are not empty
  stopifnot(nrow(red) > 0)
  stopifnot(nrow(green_filt) > 0)
  # Filtering, part 2
  blast_joined <- filter_red_green(red, green_filt)
  blast_joined
}

#### RUN ####
blast_red <- read_tsv(snakemake@input[[1]], show_col_types = FALSE)
blast_green <- parse_blast(snakemake@input[[2]], "FR_green")
output_filtered <- main(
  red = blast_red,
  green = blast_green,
  identity = snakemake@params[[1]],
  evalue = snakemake@params[[2]],
  len = snakemake@params[[3]]
)

# write to file
write_delim(output_filtered, file = snakemake@output[[1]], delim = "\t")
print("Finished. No erorrs.")

#### CLOSE LOG ####
sink()
