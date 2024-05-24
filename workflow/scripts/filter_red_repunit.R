# script to filter reads containing correct combination of
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
              help = "blast hits table for repeat unit (IS)",
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

  make_option(c("-l", "--min_len_fr"),
              type = "integer",
              default = 500,
              help = "min FR length allowed",
              metavar = "int"),

  make_option(c("-k", "--min_len_repunit"),
              type = "integer",
              default = 700,
              help = "min REP UNIT length allowed",
              metavar = "int"),

  make_option(c("-i", "--identity"),
              type = "integer",
              default = 75,
              help = "min identity",
              metavar = "integer"),

  make_option(c("-d", "--max_distance"),
              type = "integer",
              default = 30,
              help = "max allowed distance between FR and REP UNIT",
              metavar = "int")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$blast_red)) {
  print_help(opt_parser)
  stop("Input file 1 must be provided", call. = FALSE)
}
if (is.null(opt$blast_repunit)) {
  print_help(opt_parser)
  stop("Input file 2 (rep.unit) must be provided", call. = FALSE)
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

filter_blast <- function(blast_df, min_len, max_e_value, min_identity) {
  # return: a blast table filtered by minimal length, evalue and hit identity
  filt_df <- filter(blast_df,
                    length >= min_len,
                    e.value <= max_e_value,
                    identity >= min_identity)
  return(filt_df)
}

filter_rep_unit_flanking_region <-
  function(fr_df, ru_df, max_distance) {
    # filters red blast and rep unit blast by orientation and max
    # distance btw the red region and the rep.unit
    # return: filtered tibble
    left_join(fr_df, ru_df, by = 'subject') %>%
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

####

# filter rep unite table
# NB: min length is different from FR's min length
blast_red <- parse_blast(opt$blast_red, "FR_red") %>%
  filter_blast(min_len = opt$min_len_fr,
                     max_e_value = opt$evalue, 
                     min_identity = opt$identity)

blast_repunit <- parse_blast(opt$blast_repunit, "Rep_unit") %>%
  filter_blast(min_len = opt$min_len_repunit,
                     max_e_value = opt$evalue, 
                     min_identity = opt$identity)

# filter combination of both FR and rep unit
# 1st: same orientation
# 2nd: close to each other
blast_joined <-
  filter_rep_unit_flanking_region(blast_red,
                                  blast_repunit,
                                  max_distance = opt$max_distance) %>% 
  # keep some columns
  select(subject, start.subject.x, end.subject.x, 
         start.subject.y, end.subject.y, distance, orientation.x)

# give them better names
names(blast_joined) <- c("subject", "start.red", "end.red", 
                         "start.rep.unit", "end.rep.unit", "dist", "orient")

# Save results
write_delim(blast_joined, file = opt$output, delim = "\t")
print("Finished. No erorrs.")
