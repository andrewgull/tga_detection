# script to plot read length and distance between FRs

library(optparse)

# CLI parsing
option_list <- list(
  make_option(c("-b", "--blast"),
              type = "character",
              default = NULL,
              help = "joined and filtered table of blast hits of FRs",
              metavar = "file.tsv"),
  make_option(c("-r", "--reads"),
              type = "character",
              default = NULL,
              help = "Nanopore reads in FASTA format",
              metavar = "reads.fasta"),
  make_option(c("-o", "--output"),
              type = "character",
              default = NULL,
              help = "output file name (*.png)",
              metavar = "histogram.png")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$blast)) {
  print_help(opt_parser)
  stop("Input blast table must be provided", call. = FALSE)
}
if (is.null(opt$reads)) {
  print_help(opt_parser)
  stop("Input reads file must be provided", call. = FALSE)
}
if (is.null(opt$output)) {
  print_help(opt_parser)
  stop("Output file must be provided", call. = FALSE)
}

library(ggplot2)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Biostrings))

# read filtered and joined blast
blast_df <- read.table(file = opt$blast, header = TRUE, sep = "\t")
# read reads
reads <- readDNAStringSet(opt$reads)
read_lengths <- width(reads)
reads_len_df <- tibble("subject" = sub("^(.*?) runid=.*", "\\1", names(reads)),
                       "read.len" = read_lengths)

# join the two data sets
blast_df <- left_join(blast_df,
                      reads_len_df,
                      by = "subject")
# make plot
len_dist_plot <-
  ggplot(blast_df, aes(distance.btw.FR, read.len)) +
  geom_point(alpha = 0.3) +
  xlab("distance between FRs") +
  ylab("read length")

ggsave(filename = opt$output, plot = len_dist_plot, height = 7, width = 10, units = "in") # nolint: line_length_linter.
print("Plot was saved to a file. No errors.")
