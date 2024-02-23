# script to calulate and plot blaSHV counts from BED files with merged coordinates
library(optparse)

# CLI parsing
option_list <- list(
  make_option(c("-i", "--input"),
              type = "character",
              default = NULL,
              help = "BED file with merged coords",
              metavar = "input.bed"),
  make_option(c("-p", "--plot"),
              type = "character",
              default = NULL,
              help = "output plot",
              metavar = "plot.png"),
  make_option(c("-t", "--table"),
              type = "character",
              default = NULL,
              help = "output table file",
              metavar = "output.tsv")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file must be provided", call. = FALSE)
}
if (is.null(opt$plot)) {
  print_help(opt_parser)
  stop("Output plot must be provided", call. = FALSE)
}
if (is.null(opt$table)) {
  print_help(opt_parser)
  stop("Output table must be provided", call. = FALSE)
}

# other packages
suppressPackageStartupMessages(library(dplyr))
library(ggplot2)

bla_merge <- readr::read_tsv(opt$input, col_names = FALSE, show_col_types = FALSE)
# check if the input is empty
stopifnot(nrow(bla_merge) > 0)

# make counts table
counts_df <- bla_merge %>%
  group_by(X1) %>%
  count()

# make plot
counts_plot <- counts_df %>%
  ggplot(aes(n)) +
  geom_histogram(binwidth = 1) +
  geom_rug()

# save outputs
ggsave(filename = opt$plot, plot = counts_plot, height = 7, width = 10, units = "in")
readr::write_delim(x = counts_df, file = opt$table, delim = "\t")
print("Fininshed. No errors")
