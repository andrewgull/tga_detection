# plot merged bla hits

library(optparse)

# CLI parsing
option_list <- list(
  make_option(c("-i", "--input"),
              type = "character",
              default = NULL,
              help = "bed file with merged intervals",
              metavar = "file.bed"),
  make_option(c("-l", "--length"),
              type = "integer",
              default = 860,
              help = "blaSHV gene length, nt",
              metavar = "int"),
  make_option(c("-o", "--output"),
              type = "character",
              default = NULL,
              help = "plot (*.png)",
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

suppressPackageStartupMessages(library(dplyr))
library(ggplot2)
bla_merge <- read_tsv(opt$input, col_names = FALSE)

bla_merged_df <- bla_merge %>%
  group_by(X1) %>%
  summarise(sum.merged.hits = sum(X3 - X2 + 1)) %>%
  mutate(n.blaSHV.merged = round(sum.merged.hits / opt$length, 0)) %>%
  rename("subject" = X1) %>%
  # filter out the reads the 'bad reads' from previous filtering rounds
  right_join(blast_joined_filt, by = "subject") %>%
  select(subject, n.blaSHV.merged, n.blaSHV, green.red.distance.pos)

hist_plot <- bla_merged_df %>%
  ggplot(aes(n.blaSHV.merged)) +
  geom_histogram() +
  geom_rug() +
  xlab("n blaSHV merged hits")

ggsave(filename = opt$output, plot = hist_plot, height = 7, width = 10, units = "in") # nolint: line_length_linter.
print("Histogram was saved to a file. No errors.")