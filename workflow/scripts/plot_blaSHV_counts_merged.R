# plot merged bla hits

library(optparse)

# CLI parsing
option_list <- list(
  make_option(c("-i", "--input"),
              type = "character",
              default = NULL,
              help = "bed file with merged intervals",
              metavar = "file.bed"),
  make_option(c("-b", "--blast"),
              type = "character",
              default = NULL,
              help = "filtered table of blaSHV blast hits",
              metavar = "file.tsv"),
  make_option(c("-l", "--length"),
              type = "integer",
              default = 860,
              help = "blaSHV gene length, nt",
              metavar = "int"),
  make_option(c("-p", "--output_plot"),
              type = "character",
              default = NULL,
              help = "plot (*.png)",
              metavar = "histogram.png"),
  make_option(c("-a", "--output_table"),
              type = "character",
              default = NULL,
              help = "table file TSV",
              metavar = "table.tsv")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file must be provided", call. = FALSE)
}
if (is.null(opt$output_plot)) {
  print_help(opt_parser)
  stop("Output file for plot must be provided", call. = FALSE)
}
if (is.null(opt$output_table)) {
  print_help(opt_parser)
  stop("Output file for table must be provided", call. = FALSE)
}
if (is.null(opt$blast)) {
  print_help(opt_parser)
  stop("Filtered blast file must be provided", call. = FALSE)
}
suppressPackageStartupMessages(library(dplyr))
library(ggplot2)
bla_merge <- readr::read_tsv(opt$input, col_names = FALSE, show_col_types = FALSE) # nolint: line_length_linter.
blast_filt <- readr::read_tsv(opt$blast, show_col_types = FALSE) # nolint: line_length_linter.

bla_merged_df <- bla_merge %>%
  group_by(X1) %>%
  summarise(sum.merged.hits = sum(X3 - X2 + 1)) %>%
  mutate(n.blaSHV.merged = round(sum.merged.hits / opt$length, 0)) %>%
  rename("subject" = X1) %>%
  # filter out the reads the 'bad reads' from previous filtering rounds
  right_join(blast_filt, by = "subject") %>%
  select(subject, n.blaSHV.merged, distance.btw.FR)

hist_plot <- bla_merged_df %>%
  ggplot(aes(n.blaSHV.merged)) +
  geom_histogram() +
  geom_rug() +
  xlab("n blaSHV merged hits")

# save outputs
ggsave(filename = opt$output_plot, plot = hist_plot, height = 7, width = 10, units = "in") # nolint: line_length_linter.
readr::write_tsv(x = bla_merged_df, file = opt$output_table)
print("Histogram was saved to a file. No errors.")