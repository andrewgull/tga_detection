# script to parse blast tables
# and make histogram of 

library(optparse)

# CLI parsing
option_list <- list(
  make_option(c("-i", "--input1"),
              type = "character",
              default = NULL,
              help = "blast hits table for FR1",
              metavar = "character"),

  make_option(c("-e", "--input2"),
              type = "character",
              default = NULL,
              help = "blast hits table for FR2",
              metavar = "character"),

  make_option(c("-s", "--output_hist"),
              type = "character",
              default = NULL,
              help = "histogram with BL copy number distribution",
              metavar = "character"),

  make_option(c("-t", "--output_table"),
              type = "character",
              default = NULL,
              help = "table with BL copy number distribution",
              metavar = "character"),

  make_option(c("-l", "--hit_len"),
              type = "integer",
              default = 250,
              help = "minimal length of a hit (default=250)",
              metavar = "integer"),

  make_option(c("-u", "--unit_len"),
              type = "integer",
              default = 3500,
              help = "length of an amplified unit (default=3500)",
              metavar = "integer")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input1)){
  print_help(opt_parser)
  stop("Input file 1 must be provided", call. = FALSE)
}
if (is.null(opt$input2)){
  print_help(opt_parser)
  stop("Input file 2 must be provided", call. = FALSE)
}
if (is.null(opt$output_hist)){
  print_help(opt_parser)
  stop("Output file (pdf) must be provided", call. = FALSE)
}
if (is.null(opt$output_table)){
  print_help(opt_parser)
  stop("Output file (csv) must be provided", call. = FALSE)
}

# READ THE DATA

# libraries
suppressPackageStartupMessages(library(dplyr))
library(readr)
library(ggplot2)

# a function to process BLAST output
parse_blast_df <-
  function(df_name,
           minimal_length = 250,
           filter_on_positive = TRUE) {
    # read blast results
    df <-
      read_delim(df_name,
                 col_names = FALSE,
                 col_select = c(X1, X2, X3, X4, X9, X10),
                 show_col_types = FALSE)
    names(df) <-
      c("query",
        "subject",
        "identity",
        "length",
        "start.subject",
        "end.subject")
    df <- filter(df, length >= minimal_length)
    # add column indicating that the hit is on positive strain or not
    df <-
      df %>% mutate(pos_strand = (end.subject - start.subject) > 0,
                    .keep = "all")
    # find such positive strand rows that have both FR, 
    # and negative strand rows that have both FR
    if (filter_on_positive) {
      df_strand <- df %>% filter(pos_strand)
    } else {
      df_strand <- df %>% filter(!pos_strand)
    }
    return(df_strand)
  }

# a function to join processed BLAST output
join_filtered_blast <-
  function(df1_path, df2_path, min_len, on_pos) {
    parsed_list <- lapply(c(df1_path, df2_path), function(x) {
      parse_blast_df(x,
                     minimal_length = min_len,
                     filter_on_positive = on_pos)
    })
    df12_joined <-
      inner_join(parsed_list[[1]], parsed_list[[2]], by = "subject")
    return(df12_joined)
  }

# COLLECT POS STRAND ROWS W BOTH FR
fr12_pos <-
  join_filtered_blast(opt$input1,
                      opt$input2,
                      min_len = opt$hit_len,
                      on_pos = TRUE)

# COLLECT NEG STRAND ROWS W BOTH FR
fr12_neg <-
  join_filtered_blast(opt$input1,
                      opt$input2,
                      min_len = opt$hit_len,
                      on_pos = FALSE)

# BIND ROWS
fr12_all <- bind_rows(fr12_pos, fr12_neg)

# FIND CNV
n_copies <- fr12_all %>%
  mutate(n.copies = round(abs(start.subject.y - end.subject.x) / opt$unit_len))
n_copies_counts <-
  n_copies %>% group_by(n.copies) %>% summarize(counts = n())

# WRITE THE TABLE W COUNTS
write.csv(n_copies_counts,
          opt$output_table,
          row.names = FALSE,
          quote = FALSE)

# PLOT
cnv_histogram <- ggplot(n_copies, aes(n.copies)) +
  geom_histogram(binwidth = 0.5, fill = "steelblue") +
  xlab("AU") +
  ggtitle("BL-genes copy number distribution",
          subtitle = paste0("N total spanning reads ", nrow(n_copies)))

ggsave(plot = cnv_histogram, filename = opt$output_hist)