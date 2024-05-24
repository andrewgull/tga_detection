# script to filter reads with FR+RU in correct orientation
# and proximity AND at least 1820 nt from the beginning of the FR
# requires readr & dplyr

library(optparse)

#### CLI parsing ####
option_list <- list(
  make_option(c("-i", "--input"),
              type = "character",
              default = NULL,
              help = "a table file with blast results filtered by FR+RU orientation and distance",
              metavar = "character"),
  make_option(c("-o", "--output"),
              type = "character",
              default = NULL,
              help = "a table with reads at least 1820 nt in length from the start of FR",
              metavar = "character"),
  make_option(c("-d", "--distance"),
              type = "integer",
              default = 1820,
              help = "min distance from the start of FR towards RU",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input table file must be provided", call. = FALSE)
}
if (is.null(opt$output)) {
  print_help(opt_parser)
  stop("Output file must be provided", call. = FALSE)
}

#### Libraries ####
suppressPackageStartupMessages(library(dplyr))
library(readr)

#### Functions ####
filter_by_distance <- function(df, dist=1820) {
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

# read the iput table
input_table <- read_tsv(opt$input, show_col_types = FALSE)
# make size checks
stopifnot(ncol(input_table) == 7)
stopifnot(nrow(input_table) != 0)

# filter
output_table <- filter_by_distance(input_table, dist = opt$distance)

# save results
write_delim(output_table, file = opt$output, delim = "\t")
print("Finished. No erorrs.")
