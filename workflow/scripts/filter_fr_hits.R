# script to filter FR blast hits

library(optparse)

# CLI parsing
option_list <- list(
  make_option(c("-r", "--blast_red"),
              type = "character",
              default = NULL,
              help = "blast hits table for FR1 (red)",
              metavar = "character"),

  make_option(c("-g", "--blast_green"),
              type = "character",
              default = NULL,
              help = "blast hits table for FR2 (green)",
              metavar = "character"),

  make_option(c("-o", "--output"),
              type = "character",
              default = NULL,
              help = "filtered and joined blast tables",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$blast_red)){
  print_help(opt_parser)
  stop("Input file 1 must be provided", call. = FALSE)
}
if (is.null(opt$blast_green)){
  print_help(opt_parser)
  stop("Input file 2 must be provided", call. = FALSE)
}
if (is.null(opt$output)){
  print_help(opt_parser)
  stop("Output file must be provided", call. = FALSE)
}


#### Libraries ####
suppressPackageStartupMessages(library(dplyr))
library(readr)

#### Functions ####
parse_blast <- function(file_path, region_name) {
  # read blast table
  df <- read_delim(file_path, col_names = FALSE, show_col_types = FALSE)
  # assign names
  names(df) <- c("query", "subject", "identity", "length", "mismatch",
                 "gaps", "start.query", "end.query", "start.subject",
                 "end.subject", "e.value", "bit.score")
  # create orientation column
  df$orientation <- ifelse(df$start.subject < df$end.subject, "direct", "reverse") # nolint: line_length_linter.
  # rename query
  df$query <- region_name
  return(df)
}

filter_blast_p1 <- function(blast_df, min_len, max_e_value, min_identity) {
  filter(blast_df,
         length >= min_len,
         e.value <= max_e_value,
         identity >= min_identity)
}

get_multihits_ids <- function(df) {
  # df: filtered blast table
  df %>%
    group_by(subject) %>%
    count() %>%
    filter(n > 1) %>%
    pull(subject)
}

filter_blast_p2 <- function(df1, df2) {
  # find read IDs containig multiple query hits
  reads_multiple_hits <- union(get_multihits_ids(df1),
                               get_multihits_ids(df2))
  df1_filt <- df1 %>%
    filter(!subject %in% reads_multiple_hits)
  df2_filt <- df2 %>%
    filter(!subject %in% reads_multiple_hits)

  # Filter out reads with single hits
  # Filter by orientation
  # Add distance between FRs
  df_joined <- full_join(df1_filt, df2_filt, by = "subject") %>%
    filter(!is.na(query.x), !is.na(query.y), orientation.x == orientation.y) %>%
    mutate(green.red.distance = end.subject.x - start.subject.y,
           distance.btw.FR = if_else(green.red.distance < 0, green.red.distance * -1, green.red.distance * 1)) # nolint: line_length_linter.
  return(df_joined)
}

#### Main workflow ####
blast_red <- parse_blast(opt$blast_red, "FR_red")
blast_green <- parse_blast(opt$blast_green, "FR_green")

#### Filtering, part 1
filtered_blasts <- lapply(list(blast_red, blast_green), function(item) {
  filter_blast_p1(item, min_len = 1500, max_e_value = 10 ** -5, min_identity = 75) # nolint: line_length_linter.
})

### Filtering, part 2
blast_joined <- filter_blast_p2(filtered_blasts[[1]], filtered_blasts[[2]])

# Save results
write_delim(blast_joined, file = opt$output, delim = "\t")
