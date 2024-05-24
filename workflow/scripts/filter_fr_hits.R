# script to filter FR blast hits by
# presence of both RED and GREEN in correct orientation

library(optparse)

# CLI parsing
option_list <- list(
  make_option(c("-r", "--blast_red"),
              type = "character",
              default = NULL,
              help = "a table after Red+RU+1820 filtering",
              metavar = "character"),
  
  make_option(c("-g", "--blast_green"),
              type = "character",
              default = NULL,
              help = "blast hits table for FR2 (green)",
              metavar = "character"),
  
  make_option(c("-o", "--output"),
              type = "character",
              default = NULL,
              help = "filtered and joined table",
              metavar = "character"),
  
  make_option(c("-e", "--evalue"),
              type = "double",
              default = 0.00001,
              help = "max e-value",
              metavar = "double"),
  
  make_option(c("-l", "--min_len"),
              type = "integer",
              default = 500,
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
  df <- read_delim(file_path, col_names = FALSE, show_col_types = FALSE)
  # check if the table is empty
  stopifnot(nrow(df) > 0)
  # check if nrow is wrong
  stopifnot(ncol(df) == 12)
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
  # return: a blast table filtered by minimal length, evalue and hit identity
  # apply only to GREEN hits table
  filt_df <- filter(blast_df,
                    length >= min_len,
                    e.value <= max_e_value,
                    identity >= min_identity)
  return(filt_df)
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
  reads_multiple_hits <- union(get_multihits_ids(df_red_ru),
                               get_multihits_ids(df_green))
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
    mutate(green.red.distance = end.red - start.subject,
           distance.btw.FR = if_else(green.red.distance < 0, 
                                     green.red.distance * -1, 
                                     green.red.distance * 1))
  return(df_joined)
}

#### 

#### read and basic filtering if needed
blast_red <- read_tsv(opt$blast_red, show_col_types = FALSE)
blast_red$query <- "FR_RU_filt"
blast_green <- parse_blast(opt$blast_green, "FR_green") %>% 
  # apply basic filtering of blast results
  filter_blast_p1(min_len = opt$min_len, 
                  max_e_value = opt$evalue, 
                  min_identity = opt$identity)

# check if they are not empty
stopifnot(nrow(blast_red) > 0)
stopifnot(nrow(blast_green) > 0)

### Filtering, part 2
blast_joined <- filter_red_green(blast_red, blast_green)

# Save results
write_delim(blast_joined, file = opt$output, delim = "\t")
print("Finished. No erorrs.")