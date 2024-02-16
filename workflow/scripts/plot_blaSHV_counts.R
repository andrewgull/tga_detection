# script to plot blaSHV counts
library(optparse)

# CLI parsing
option_list <- list(
  make_option(c("-i", "--input"),
              type = "character",
              default = NULL,
              help = "filtered table of blaSHV blast hits",
              metavar = "file.tsv"),
  make_option(c("-l", "--length"),
              type = "integer",
              default = 860,
              help = "blaSHV gene length, nt",
              metavar = "int"),
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
suppressPackageStartupMessages(library(dplyr))

# gene length
bla_len <- opt$length
# read the filtered blast table
blast_blaSHV <- readr::read_tsv(opt$input, show_col_types = FALSE)
# check if it is empty
stopifnot(nrow(blast_blaSHV) > 0)

# n bla = (sum of bla hits/bla length)
hist_plot <- blast_blaSHV %>%
  group_by(subject) %>%
  summarise(sum.bla.hit.len = sum(length)) %>%
  mutate(n.blaSHV = round(sum.bla.hit.len / bla_len)) %>%
  ggplot(aes(n.blaSHV)) +
  geom_histogram(bins = 100, fill = "darkslategray4") +
  geom_rug() +
  xlab("n blaSHV (sum hits)")

ggsave(filename = opt$output, plot = hist_plot, height = 7, width = 10, units = "in") # nolint: line_length_linter.
print("Histogram was saved to a file. No errors.")
