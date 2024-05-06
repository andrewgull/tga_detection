# script to filter reads containng correct combiation of
# the red flanking region and the repeat unit
# the 'correct combination' that is 'a red FR first and a repeat unit right next to it'

library(optparse)

#### CLI parsing ####
option_list <- list(
  make_option(c("-r", "--blast_red"),
              type = "character",
              default = NULL,
              help = "blast hits table for FR1 (red)",
              metavar = "character"),

  make_option(c("-u", "--blast_repunit"),
              type = "character",
              default = NULL,
              help = "blast hits table for FR2 (green)",
              metavar = "character"),

  make_option(c("-o", "--output"),
              type = "character",
              default = NULL,
              help = "filtered and joined blast tables",
              metavar = "character"),

  make_option(c("-e", "--evalue"),
              type = "double",
              default = 0.00001,
              help = "max e-value",
              metavar = "double"),

  make_option(c("-l", "--min_len"),
              type = "integer",
              default = 1500,
              help = "min hit length",
              metavar = "int"),

  make_option(c("-i", "--identity"),
              type = "integer",
              default = 75,
              help = "min identity",
              metavar = "integer")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$blast_red)) {
  print_help(opt_parser)
  stop("Input file 1 must be provided", call. = FALSE)
}
if (is.null(opt$blast_green)) {
  print_help(opt_parser)
  stop("Input file 2 must be provided", call. = FALSE)
}
if (is.null(opt$output)) {
  print_help(opt_parser)
  stop("Output file must be provided", call. = FALSE)
}

#### Libraries ####
suppressPackageStartupMessages(library(dplyr))
library(readr)

#### Functions ####
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
    ifelse(df$start.subject < df$end.subject, "direct", "reverse") # nolint: line_length_linter.
  # rename query
  df$query <- region_name
  return(df)
}

filter_blast_part1 <- function(blast_df, min_len, max_e_value, min_identity) {
  # return: a blast table filtered by minimal length, evalue and hit identity
  filt_df <- filter(blast_df,
                    length >= min_len,
                    e.value <= max_e_value,
                    identity >= min_identity)
  return(filt_df)
}

filter_rep_unit_flanking_region <-
  function(fr_df, ru_df, min_distance) {
    left_join(blast_red, blast_repunit, by = 'subject') %>%
      filter(orientation.x == orientation.y) %>%
      mutate(
        distance = if_else(
          orientation.x == "reverse",
          end.subject.y - start.subject.x + 1,
          start.subject.x - end.subject.y + 1
        )
      ) %>%
      filter(abs(distance) <= min_distance)
  }

####

# filter rep unite table
# NB: min length is different from FR's min length
blast_red <- parse_blast("results/tables/76595_D_mh/blast_red.tsv", "FR_red") %>% 
  filter_blast_part1(min_len = 500, max_e_value = 0.00001, min_identity = 0.75)

blast_repunit <- parse_blast("results/tables/76595_D_mh/blast_repeat_unit.tsv", "Rep_unit") %>% 
  filter_blast_part1(min_len = 700, max_e_value = 0.00001, min_identity = 0.75)

# filter combination of both FR and rep unit
# 1st: same orientation
# 2nd: next to each other

red_repunit <- filter_rep_unit_flanking_region(blast_red, blast_repunit, min_distance = 20)

# write it down?
# what the next script expects as input?

