# script to plot a histogram of distnces between red and green flanking regions 

library(optparse)

# CLI parsing
option_list <- list(
  make_option(c("-i", "--input"),
              type = "character",
              default = NULL,
              help = "joined and filtered table of blast hits of FRs",
              metavar = "file.tsv"),

  make_option(c("-o", "--output"),
              type = "character",
              default = NULL,
              help = "output file name (*.png)",
              metavar = "histogram.png")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file must be provided", call. = FALSE)
}
if (is.null(opt$output)) {
  print_help(opt_parser)
  stop("Output file must be provided", call. = FALSE)
}

library(ggplot2)

blast_df <- readr::read_delim(opt$input, show_col_types = FALSE)

hist_plot <- ggplot(blast_df, aes(distance.btw.FR)) +
  geom_histogram(bins = 200, fill = "steelblue") +
  geom_rug() +
  xlab("distance between FRs")

ggsave(filename = opt$output, plot = hist_plot, height = 7, width = 10, units = "in") # nolint: line_length_linter.
print("Histogram was saved to a file. No errors.")
